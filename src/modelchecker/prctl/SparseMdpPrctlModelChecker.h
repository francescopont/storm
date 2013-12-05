/*
 * SparseMdpPrctlModelChecker.h
 *
 *  Created on: 15.02.2013
 *      Author: Christian Dehnert
 */

#ifndef STORM_MODELCHECKER_PRCTL_SPARSEMDPPRCTLMODELCHECKER_H_
#define STORM_MODELCHECKER_PRCTL_SPARSEMDPPRCTLMODELCHECKER_H_

#include <vector>
#include <stack>
#include <fstream>

#include "src/modelchecker/prctl/AbstractModelChecker.h"
#include "src/solver/AbstractNondeterministicLinearEquationSolver.h"
#include "src/solver/GmmxxLinearEquationSolver.h"
#include "src/models/Mdp.h"
#include "src/utility/vector.h"
#include "src/utility/graph.h"
#include "src/utility/solver.h"
#include "src/settings/Settings.h"
#include "src/storage/TotalScheduler.h"

namespace storm {
    namespace modelchecker {
        namespace prctl {
            
            /*!
             * This class represents the base class for all PRCTL model checkers for MDPs.
             */
            template<class Type>
            class SparseMdpPrctlModelChecker : public AbstractModelChecker<Type> {
                
            public:
                /*!
                 * Constructs a SparseMdpPrctlModelChecker with the given model.
                 *
                 * @param model The MDP to be checked.
                 */
                explicit SparseMdpPrctlModelChecker(storm::models::Mdp<Type> const& model, std::shared_ptr<storm::solver::AbstractNondeterministicLinearEquationSolver<Type>> nondeterministicLinearEquationSolver = storm::utility::solver::getNondeterministicLinearEquationSolver<Type>()) : AbstractModelChecker<Type>(model), minimumOperatorStack(), nondeterministicLinearEquationSolver(nondeterministicLinearEquationSolver) {
                    // Intentionally left empty.
                }
                
                /*!
                 * Copy constructs a SparseMdpPrctlModelChecker from the given model checker. In particular, this means that the newly
                 * constructed model checker will have the model of the given model checker as its associated model.
                 */
                explicit SparseMdpPrctlModelChecker(storm::modelchecker::prctl::SparseMdpPrctlModelChecker<Type> const& modelchecker)
                : AbstractModelChecker<Type>(modelchecker),  minimumOperatorStack(), nondeterministicLinearEquationSolver(storm::utility::solver::getNondeterministicLinearEquationSolver<Type>()) {
                    // Intentionally left empty.
                }
                
                /*!
                 * Returns a constant reference to the MDP associated with this model checker.
                 * @returns A constant reference to the MDP associated with this model checker.
                 */
                storm::models::Mdp<Type> const& getModel() const {
                    return AbstractModelChecker<Type>::template getModel<storm::models::Mdp<Type>>();
                }
                
                /*!
                 * Checks the given formula that is a P/R operator without a bound.
                 *
                 * @param formula The formula to check.
                 * @returns The set of states satisfying the formula represented by a bit vector.
                 */
                virtual std::vector<Type> checkNoBoundOperator(const storm::property::prctl::AbstractNoBoundOperator<Type>& formula) const {
                    // Check if the operator was an non-optimality operator and report an error in that case.
                    if (!formula.isOptimalityOperator()) {
                        LOG4CPLUS_ERROR(logger, "Formula does not specify neither min nor max optimality, which is not meaningful over nondeterministic models.");
                        throw storm::exceptions::InvalidArgumentException() << "Formula does not specify neither min nor max optimality, which is not meaningful over nondeterministic models.";
                    }
                    minimumOperatorStack.push(formula.isMinimumOperator());
                    std::vector<Type> result = formula.check(*this, false);
                    minimumOperatorStack.pop();
                    return result;
                }
                
                /*!
                 * Computes the probability to satisfy phi until psi within a limited number of steps for each state.
                 *
                 * @param phiStates A bit vector indicating which states satisfy phi.
                 * @param psiStates A bit vector indicating which states satisfy psi.
                 * @param stepBound The upper bound for the number of steps.
                 * @param qualitative A flag indicating whether the check only needs to be done qualitatively, i.e. if the
                 * results are only compared against the bounds 0 and 1. If set to true, this will most likely results that are only
                 * qualitatively correct, i.e. do not represent the correct value, but only the correct relation with respect to the
                 * bounds 0 and 1.
                 * @return The probabilities for satisfying phi until psi within a limited number of steps for each state.
                 * If the qualitative flag is set, exact probabilities might not be computed.
                 */
                std::vector<Type> checkBoundedUntil(storm::storage::BitVector const& phiStates, storm::storage::BitVector const& psiStates, uint_fast64_t stepBound, bool qualitative) const {
                    std::vector<Type> result(this->getModel().getNumberOfStates());

                    // Determine the states that have 0 probability of reaching the target states.
                    storm::storage::BitVector statesWithProbabilityGreater0;
                    if (this->minimumOperatorStack.top()) {
                        statesWithProbabilityGreater0 = storm::utility::graph::performProbGreater0A(this->getModel().getTransitionMatrix(), this->getModel().getNondeterministicChoiceIndices(), this->getModel().getBackwardTransitions(), phiStates, psiStates, true, stepBound);
                    } else {
                        statesWithProbabilityGreater0 = storm::utility::graph::performProbGreater0E(this->getModel().getTransitionMatrix(), this->getModel().getNondeterministicChoiceIndices(), this->getModel().getBackwardTransitions(), phiStates, psiStates, true, stepBound);
                    }
                    
                    // Check if we already know the result (i.e. probability 0) for all initial states and
                    // don't compute anything in this case.
                    if (this->getModel().getInitialStates().isDisjointFrom(statesWithProbabilityGreater0)) {
                        LOG4CPLUS_INFO(logger, "The probabilities for the initial states were determined in a preprocessing step."
                                       << " No exact probabilities were computed.");
                        // Set the values for all maybe-states to 0.5 to indicate that their probability values are not 0 (and
                        // not necessarily 1).
                        storm::utility::vector::setVectorValues(result, statesWithProbabilityGreater0, Type(0.5));
                    } else {
                        // In this case we have have to compute the probabilities.
                        
                        // We can eliminate the rows and columns from the original transition probability matrix that have probability 0.
                        storm::storage::SparseMatrix<Type> submatrix = this->getModel().getTransitionMatrix().getSubmatrix(statesWithProbabilityGreater0, this->getModel().getNondeterministicChoiceIndices());
                        
                        // Get the "new" nondeterministic choice indices for the submatrix.
                        std::vector<uint_fast64_t> subNondeterministicChoiceIndices = storm::utility::vector::getConstrainedOffsetVector(this->getModel().getNondeterministicChoiceIndices(), statesWithProbabilityGreater0);
                        
                        // Compute the new set of target states in the reduced system.
                        storm::storage::BitVector rightStatesInReducedSystem = psiStates % statesWithProbabilityGreater0;
                        
                        // Make all rows absorbing that satisfy the second sub-formula.
                        submatrix.makeRowsAbsorbing(rightStatesInReducedSystem, subNondeterministicChoiceIndices);
                        
                        // Create the vector with which to multiply.
                        std::vector<Type> subresult(statesWithProbabilityGreater0.getNumberOfSetBits());
                        storm::utility::vector::setVectorValues(subresult, rightStatesInReducedSystem, storm::utility::constGetOne<Type>());
                        
                        this->nondeterministicLinearEquationSolver->performMatrixVectorMultiplication(this->minimumOperatorStack.top(), submatrix, subresult, subNondeterministicChoiceIndices, nullptr, stepBound);
                        
                        // Set the values of the resulting vector accordingly.
                        storm::utility::vector::setVectorValues(result, statesWithProbabilityGreater0, subresult);
                        storm::utility::vector::setVectorValues(result, ~statesWithProbabilityGreater0, storm::utility::constGetZero<Type>());
                    }
                    
                    return result;
                }
                
                /*!
                 * Checks the given formula that is a bounded-until formula.
                 *
                 * @param formula The formula to check.
                 * @param qualitative A flag indicating whether the formula only needs to be evaluated qualitatively, i.e. if the
                 * results are only compared against the bounds 0 and 1. If set to true, this will most likely results that are only
                 * qualitatively correct, i.e. do not represent the correct value, but only the correct relation with respect to the
                 * bounds 0 and 1.
                 * @return The probabilities for the given formula to hold on every state of the model associated with this model
                 * checker. If the qualitative flag is set, exact probabilities might not be computed.
                 */
                virtual std::vector<Type> checkBoundedUntil(storm::property::prctl::BoundedUntil<Type> const& formula, bool qualitative) const {
                    return checkBoundedUntil(formula.getLeft().check(*this), formula.getRight().check(*this), formula.getBound(), qualitative);
                }
                
                /*!
                 * Computes the probability to reach the given set of states in the next step for each state.
                 *
                 * @param nextStates A bit vector defining the states to reach in the next state.
                 * @param qualitative A flag indicating whether the formula only needs to be evaluated qualitatively, i.e. if the
                 * results are only compared against the bounds 0 and 1. If set to true, this will most likely results that are only
                 * qualitatively correct, i.e. do not represent the correct value, but only the correct relation with respect to the
                 * bounds 0 and 1.
                 * @return The probabilities to reach the gien set of states in the next step for each state. If the
                 * qualitative flag is set, exact probabilities might not be computed.
                 */
                virtual std::vector<Type> checkNext(storm::storage::BitVector const& nextStates, bool qualitative) const {
                    // Create the vector with which to multiply and initialize it correctly.
                    std::vector<Type> result(this->getModel().getNumberOfStates());
                    storm::utility::vector::setVectorValues(result, nextStates, storm::utility::constGetOne<Type>());
                    
                    this->nondeterministicLinearEquationSolver->performMatrixVectorMultiplication(this->minimumOperatorStack.top(), this->getModel().getTransitionMatrix(), result, this->getModel().getNondeterministicChoiceIndices());
                    
                    return result;
                }
                
                /*!
                 * Checks the given formula that is a next formula.
                 *
                 * @param formula The formula to check.
                 * @param qualitative A flag indicating whether the formula only needs to be evaluated qualitatively, i.e. if the
                 * results are only compared against the bounds 0 and 1. If set to true, this will most likely results that are only
                 * qualitatively correct, i.e. do not represent the correct value, but only the correct relation with respect to the
                 * bounds 0 and 1.
                 * @return The probabilities for the given formula to hold on every state of the model associated with this model
                 * checker. If the qualitative flag is set, exact probabilities might not be computed.
                 */
                virtual std::vector<Type> checkNext(const storm::property::prctl::Next<Type>& formula, bool qualitative) const {
                    return checkNext(formula.getChild().check(*this), qualitative);
                }
                
                /*!
                 * Checks the given formula that is a bounded-eventually formula.
                 *
                 * @param formula The formula to check.
                 * @param qualitative A flag indicating whether the formula only needs to be evaluated qualitatively, i.e. if the
                 * results are only compared against the bounds 0 and 1. If set to true, this will most likely results that are only
                 * qualitatively correct, i.e. do not represent the correct value, but only the correct relation with respect to the
                 * bounds 0 and 1.
                 * @returns The probabilities for the given formula to hold on every state of the model associated with this model
                 * checker. If the qualitative flag is set, exact probabilities might not be computed.
                 */
                virtual std::vector<Type> checkBoundedEventually(const storm::property::prctl::BoundedEventually<Type>& formula, bool qualitative) const {
                    // Create equivalent temporary bounded until formula and check it.
                    storm::property::prctl::BoundedUntil<Type> temporaryBoundedUntilFormula(new storm::property::prctl::Ap<Type>("true"), formula.getChild().clone(), formula.getBound());
                    return this->checkBoundedUntil(temporaryBoundedUntilFormula, qualitative);
                }
                
                /*!
                 * Checks the given formula that is an eventually formula.
                 *
                 * @param formula The formula to check.
                 * @param qualitative A flag indicating whether the formula only needs to be evaluated qualitatively, i.e. if the
                 * results are only compared against the bounds 0 and 1. If set to true, this will most likely results that are only
                 * qualitatively correct, i.e. do not represent the correct value, but only the correct relation with respect to the
                 * bounds 0 and 1.
                 * @returns The probabilities for the given formula to hold on every state of the model associated with this model
                 * checker. If the qualitative flag is set, exact probabilities might not be computed.
                 */
                virtual std::vector<Type> checkEventually(const storm::property::prctl::Eventually<Type>& formula, bool qualitative) const {
                    // Create equivalent temporary until formula and check it.
                    storm::property::prctl::Until<Type> temporaryUntilFormula(new storm::property::prctl::Ap<Type>("true"), formula.getChild().clone());
                    return this->checkUntil(temporaryUntilFormula, qualitative);
                }
                
                /*!
                 * Checks the given formula that is a globally formula.
                 *
                 * @param formula The formula to check.
                 * @param qualitative A flag indicating whether the formula only needs to be evaluated qualitatively, i.e. if the
                 * results are only compared against the bounds 0 and 1. If set to true, this will most likely results that are only
                 * qualitatively correct, i.e. do not represent the correct value, but only the correct relation with respect to the
                 * bounds 0 and 1.
                 * @returns The probabilities for the given formula to hold on every state of the model associated with this model
                 * checker. If the qualitative flag is set, exact probabilities might not be computed.
                 */
                virtual std::vector<Type> checkGlobally(const storm::property::prctl::Globally<Type>& formula, bool qualitative) const {
                    // Create "equivalent" temporary eventually formula and check it.
                    storm::property::prctl::Eventually<Type> temporaryEventuallyFormula(new storm::property::prctl::Not<Type>(formula.getChild().clone()));
                    std::vector<Type> result = this->checkEventually(temporaryEventuallyFormula, qualitative);
                    
                    // Now subtract the resulting vector from the constant one vector to obtain final result.
                    storm::utility::vector::subtractFromConstantOneVector(result);
                    return result;
                }
                
                /*!
                 * Check the given formula that is an until formula.
                 *
                 * @param formula The formula to check.
                 * @param qualitative A flag indicating whether the formula only needs to be evaluated qualitatively, i.e. if the
                 * results are only compared against the bounds 0 and 1. If set to true, this will most likely results that are only
                 * qualitatively correct, i.e. do not represent the correct value, but only the correct relation with respect to the
                 * bounds 0 and 1.
                 * @param scheduler If <code>qualitative</code> is false and this vector is non-null and has as many elements as
                 * there are states in the MDP, this vector will represent a scheduler for the model that achieves the probability
                 * returned by model checking. To this end, the vector will hold the nondeterministic choice made for each state.
                 * @return The probabilities for the given formula to hold on every state of the model associated with this model
                 * checker. If the qualitative flag is set, exact probabilities might not be computed.
                 */
                virtual std::vector<Type> checkUntil(const storm::property::prctl::Until<Type>& formula, bool qualitative) const {
                    return this->checkUntil(this->minimumOperatorStack.top(), formula.getLeft().check(*this), formula.getRight().check(*this), qualitative).first;
                }
                
                /*!
                 * Computes the extremal probability to satisfy phi until psi for each state in the model.
                 *
                 * @param minimize If set, the probability is minimized and maximized otherwise.
                 * @param phiStates A bit vector indicating which states satisfy phi.
                 * @param psiStates A bit vector indicating which states satisfy psi.
                 * @param qualitative A flag indicating whether the formula only needs to be evaluated qualitatively, i.e. if the
                 * results are only compared against the bounds 0 and 1. If set to true, this will most likely results that are only
                 * qualitatively correct, i.e. do not represent the correct value, but only the correct relation with respect to the
                 * bounds 0 and 1.
                 * @param scheduler If <code>qualitative</code> is false and this vector is non-null and has as many elements as
                 * there are states in the MDP, this vector will represent a scheduler for the model that achieves the probability
                 * returned by model checking. To this end, the vector will hold the nondeterministic choice made for each state.
                 * @return The probabilities for the satisfying phi until psi for each state of the model. If the
                 * qualitative flag is set, exact probabilities might not be computed.
                 */
                static std::pair<std::vector<Type>, storm::storage::TotalScheduler> computeUnboundedUntilProbabilities(bool minimize, storm::storage::SparseMatrix<Type> const& transitionMatrix, std::vector<uint_fast64_t> nondeterministicChoiceIndices, storm::storage::SparseMatrix<Type> const& backwardTransitions, storm::storage::BitVector const& initialStates, storm::storage::BitVector const& phiStates, storm::storage::BitVector const& psiStates, std::shared_ptr<storm::solver::AbstractNondeterministicLinearEquationSolver<Type>> nondeterministicLinearEquationSolver, bool qualitative) {
                    size_t numberOfStates = phiStates.size();
                    
                    // We need to identify the states which have to be taken out of the matrix, i.e.
                    // all states that have probability 0 and 1 of satisfying the until-formula.
                    std::pair<storm::storage::BitVector, storm::storage::BitVector> statesWithProbability01;
                    if (minimize) {
                        statesWithProbability01 = storm::utility::graph::performProb01Min(transitionMatrix, nondeterministicChoiceIndices, backwardTransitions, phiStates, psiStates);
                    } else {
                        statesWithProbability01 = storm::utility::graph::performProb01Max(transitionMatrix, nondeterministicChoiceIndices, backwardTransitions, phiStates, psiStates);
                    }
                    storm::storage::BitVector statesWithProbability0 = std::move(statesWithProbability01.first);
                    storm::storage::BitVector statesWithProbability1 = std::move(statesWithProbability01.second);
                    
                    storm::storage::BitVector maybeStates = ~(statesWithProbability0 | statesWithProbability1);
                    LOG4CPLUS_INFO(logger, "Found " << statesWithProbability0.getNumberOfSetBits() << " 'no' states.");
                    LOG4CPLUS_INFO(logger, "Found " << statesWithProbability1.getNumberOfSetBits() << " 'yes' states.");
                    LOG4CPLUS_INFO(logger, "Found " << maybeStates.getNumberOfSetBits() << " 'maybe' states.");
                    
                    // Create resulting vector.
                    std::vector<Type> result(numberOfStates);
                    
                    // Check whether we need to compute exact probabilities for some states.
                    if (initialStates.isDisjointFrom(maybeStates) || qualitative) {
                        if (qualitative) {
                            LOG4CPLUS_INFO(logger, "The formula was checked qualitatively. No exact probabilities were computed.");
                        } else {
                            LOG4CPLUS_INFO(logger, "The probabilities for the initial states were determined in a preprocessing step."
                                           << " No exact probabilities were computed.");
                        }
                        // Set the values for all maybe-states to 0.5 to indicate that their probability values
                        // are neither 0 nor 1.
                        storm::utility::vector::setVectorValues<Type>(result, maybeStates, Type(0.5));
                    } else {
                        // In this case we have have to compute the probabilities.
                        
                        // First, we can eliminate the rows and columns from the original transition probability matrix for states
                        // whose probabilities are already known.
                        storm::storage::SparseMatrix<Type> submatrix = transitionMatrix.getSubmatrix(maybeStates, nondeterministicChoiceIndices);
                        
                        // Get the "new" nondeterministic choice indices for the submatrix.
                        std::vector<uint_fast64_t> subNondeterministicChoiceIndices = storm::utility::vector::getConstrainedOffsetVector(nondeterministicChoiceIndices, maybeStates);
                        
                        // Prepare the right-hand side of the equation system. For entry i this corresponds to
                        // the accumulated probability of going from state i to some 'yes' state.
                        std::vector<Type> b = transitionMatrix.getConstrainedRowSumVector(maybeStates, nondeterministicChoiceIndices, statesWithProbability1, submatrix.getRowCount());
                        
                        // Create vector for results for maybe states.
                        std::vector<Type> x(maybeStates.getNumberOfSetBits());
                        
                        // Solve the corresponding system of equations.
                        nondeterministicLinearEquationSolver->solveEquationSystem(minimize, submatrix, x, b, subNondeterministicChoiceIndices);
                        
                        // Set values of resulting vector according to result.
                        storm::utility::vector::setVectorValues<Type>(result, maybeStates, x);
                    }
                    
                    // Set values of resulting vector that are known exactly.
                    storm::utility::vector::setVectorValues<Type>(result, statesWithProbability0, storm::utility::constGetZero<Type>());
                    storm::utility::vector::setVectorValues<Type>(result, statesWithProbability1, storm::utility::constGetOne<Type>());
                    
                    // Finally, compute a scheduler that achieves the extramal value.
                    storm::storage::TotalScheduler scheduler = computeExtremalScheduler(minimize, transitionMatrix, nondeterministicChoiceIndices, result);

                    return std::make_pair(result, scheduler);
                }
                
                std::pair<std::vector<Type>, storm::storage::TotalScheduler> checkUntil(bool minimize, storm::storage::BitVector const& phiStates, storm::storage::BitVector const& psiStates, bool qualitative) const {
                    return computeUnboundedUntilProbabilities(minimize, this->getModel().getTransitionMatrix(), this->getModel().getNondeterministicChoiceIndices(), this->getModel().getBackwardTransitions(), this->getModel().getInitialStates(), phiStates, psiStates, this->nondeterministicLinearEquationSolver, qualitative);
                }
                
                /*!
                 * Checks the given formula that is an instantaneous reward formula.
                 *
                 * @param formula The formula to check.
                 * @param qualitative A flag indicating whether the formula only needs to be evaluated qualitatively, i.e. if the
                 * results are only compared against the bound 0. If set to true, this will most likely results that are only
                 * qualitatively correct, i.e. do not represent the correct value, but only the correct relation with respect to the
                 * bound 0.
                 * @return The reward values for the given formula for every state of the model associated with this model
                 * checker. If the qualitative flag is set, exact values might not be computed.
                 */
                virtual std::vector<Type> checkInstantaneousReward(const storm::property::prctl::InstantaneousReward<Type>& formula, bool qualitative) const {
                    // Only compute the result if the model has a state-based reward model.
                    if (!this->getModel().hasStateRewards()) {
                        LOG4CPLUS_ERROR(logger, "Missing (state-based) reward model for formula.");
                        throw storm::exceptions::InvalidPropertyException() << "Missing (state-based) reward model for formula.";
                    }
                    
                    // Initialize result to state rewards of the model.
                    std::vector<Type> result(this->getModel().getStateRewardVector());
                    
                    this->nondeterministicLinearEquationSolver->performMatrixVectorMultiplication(this->minimumOperatorStack.top(), this->getModel().getTransitionMatrix(), result, this->getModel().getNondeterministicChoiceIndices(), nullptr, formula.getBound());
                    
                    return result;
                }
                
                /*!
                 * Check the given formula that is a cumulative reward formula.
                 *
                 * @param formula The formula to check.
                 * @param qualitative A flag indicating whether the formula only needs to be evaluated qualitatively, i.e. if the
                 * results are only compared against the bound 0. If set to true, this will most likely results that are only
                 * qualitatively correct, i.e. do not represent the correct value, but only the correct relation with respect to the
                 * bound 0.
                 * @return The reward values for the given formula for every state of the model associated with this model
                 * checker. If the qualitative flag is set, exact values might not be computed.
                 */
                virtual std::vector<Type> checkCumulativeReward(const storm::property::prctl::CumulativeReward<Type>& formula, bool qualitative) const {
                    // Only compute the result if the model has at least one reward model.
                    if (!this->getModel().hasStateRewards() && !this->getModel().hasTransitionRewards()) {
                        LOG4CPLUS_ERROR(logger, "Missing reward model for formula.");
                        throw storm::exceptions::InvalidPropertyException() << "Missing reward model for formula.";
                    }
                    
                    // Compute the reward vector to add in each step based on the available reward models.
                    std::vector<Type> totalRewardVector;
                    if (this->getModel().hasTransitionRewards()) {
                        totalRewardVector = this->getModel().getTransitionMatrix().getPointwiseProductRowSumVector(this->getModel().getTransitionRewardMatrix());
                        if (this->getModel().hasStateRewards()) {
                            storm::utility::vector::addVectorsInPlace(totalRewardVector, this->getModel().getStateRewardVector());
                        }
                    } else {
                        totalRewardVector = std::vector<Type>(this->getModel().getStateRewardVector());
                    }
                    
                    // Initialize result to either the state rewards of the model or the null vector.
                    std::vector<Type> result;
                    if (this->getModel().hasStateRewards()) {
                        result = std::vector<Type>(this->getModel().getStateRewardVector());
                    } else {
                        result.resize(this->getModel().getNumberOfStates());
                    }
                    
                    this->nondeterministicLinearEquationSolver->performMatrixVectorMultiplication(this->minimumOperatorStack.top(), this->getModel().getTransitionMatrix(), result, this->getModel().getNondeterministicChoiceIndices(), &totalRewardVector, formula.getBound());
                    
                    return result;
                }
                
                /*!
                 * Checks the given formula that is a reachability reward formula.
                 *
                 * @param formula The formula to check.
                 * @param qualitative A flag indicating whether the formula only needs to be evaluated qualitatively, i.e. if the
                 * results are only compared against the bound 0. If set to true, this will most likely results that are only
                 * qualitatively correct, i.e. do not represent the correct value, but only the correct relation with respect to the
                 * bound 0.
                 * @return The reward values for the given formula for every state of the model associated with this model
                 * checker. If the qualitative flag is set, exact values might not be computed.
                 */
                virtual std::vector<Type> checkReachabilityReward(const storm::property::prctl::ReachabilityReward<Type>& formula, bool qualitative) const {
                    return this->checkReachabilityReward(this->minimumOperatorStack.top(), formula.getChild().check(*this), qualitative).first;
                }
                
                /*!
                 * Computes the expected reachability reward that is gained before a target state is reached for each state.
                 *
                 * @param minimize If set, the reward is to be minimized and maximized otherwise.
                 * @param targetStates The target states before which rewards can be gained.
                 * @param qualitative A flag indicating whether the formula only needs to be evaluated qualitatively, i.e. if the
                 * results are only compared against the bound 0. If set to true, this will most likely results that are only
                 * qualitatively correct, i.e. do not represent the correct value, but only the correct relation with respect to the
                 * bound 0.
                 * @param scheduler If <code>qualitative</code> is false and this vector is non-null and has as many elements as
                 * there are states in the MDP, this vector will represent a scheduler for the model that achieves the probability
                 * returned by model checking. To this end, the vector will hold the nondeterministic choice made for each state.
                 * @return The expected reward values gained before a target state is reached for each state. If the
                 * qualitative flag is set, exact values might not be computed.
                 */
                virtual std::pair<std::vector<Type>, storm::storage::TotalScheduler> checkReachabilityReward(bool minimize, storm::storage::BitVector const& targetStates, bool qualitative) const {
                    // Only compute the result if the model has at least one reward model.
                    if (!(this->getModel().hasStateRewards() || this->getModel().hasTransitionRewards())) {
                        LOG4CPLUS_ERROR(logger, "Missing reward model for formula.");
                        throw storm::exceptions::InvalidPropertyException() << "Missing reward model for formula.";
                    }
                    
                    // Determine which states have a reward of infinity by definition.
                    storm::storage::BitVector infinityStates;
                    storm::storage::BitVector trueStates(this->getModel().getNumberOfStates(), true);
                    if (minimize) {
                        infinityStates = std::move(storm::utility::graph::performProb1A(this->getModel().getTransitionMatrix(), this->getModel().getNondeterministicChoiceIndices(), this->getModel().getBackwardTransitions(), trueStates, targetStates));
                    } else {
                        infinityStates = std::move(storm::utility::graph::performProb1E(this->getModel().getTransitionMatrix(), this->getModel().getNondeterministicChoiceIndices(), this->getModel().getBackwardTransitions(), trueStates, targetStates));
                    }
                    infinityStates.complement();
                    storm::storage::BitVector maybeStates = ~targetStates & ~infinityStates;
                    LOG4CPLUS_INFO(logger, "Found " << infinityStates.getNumberOfSetBits() << " 'infinity' states.");
                    LOG4CPLUS_INFO(logger, "Found " << targetStates.getNumberOfSetBits() << " 'target' states.");
                    LOG4CPLUS_INFO(logger, "Found " << maybeStates.getNumberOfSetBits() << " 'maybe' states.");
                    
                    // Create resulting vector.
                    std::vector<Type> result(this->getModel().getNumberOfStates());
                    
                    // Check whether we need to compute exact rewards for some states.
                    if (this->getModel().getInitialStates().isDisjointFrom(maybeStates)) {
                        LOG4CPLUS_INFO(logger, "The rewards for the initial states were determined in a preprocessing step."
                                       << " No exact rewards were computed.");
                        // Set the values for all maybe-states to 1 to indicate that their reward values
                        // are neither 0 nor infinity.
                        storm::utility::vector::setVectorValues<Type>(result, maybeStates, storm::utility::constGetOne<Type>());
                    } else {
                        // In this case we have to compute the reward values for the remaining states.
                        
                        // We can eliminate the rows and columns from the original transition probability matrix for states
                        // whose reward values are already known.
                        storm::storage::SparseMatrix<Type> submatrix = this->getModel().getTransitionMatrix().getSubmatrix(maybeStates, this->getModel().getNondeterministicChoiceIndices());
                        
                        // Get the "new" nondeterministic choice indices for the submatrix.
                        std::vector<uint_fast64_t> subNondeterministicChoiceIndices = storm::utility::vector::getConstrainedOffsetVector(this->getModel().getNondeterministicChoiceIndices(), maybeStates);
                        
                        // Prepare the right-hand side of the equation system. For entry i this corresponds to
                        // the accumulated probability of going from state i to some 'yes' state.
                        std::vector<Type> b(submatrix.getRowCount());
                        
                        if (this->getModel().hasTransitionRewards()) {
                            // If a transition-based reward model is available, we initialize the right-hand
                            // side to the vector resulting from summing the rows of the pointwise product
                            // of the transition probability matrix and the transition reward matrix.
                            std::vector<Type> pointwiseProductRowSumVector = this->getModel().getTransitionMatrix().getPointwiseProductRowSumVector(this->getModel().getTransitionRewardMatrix());
                            storm::utility::vector::selectVectorValues(b, maybeStates, this->getModel().getNondeterministicChoiceIndices(), pointwiseProductRowSumVector);
                            
                            if (this->getModel().hasStateRewards()) {
                                // If a state-based reward model is also available, we need to add this vector
                                // as well. As the state reward vector contains entries not just for the states
                                // that we still consider (i.e. maybeStates), we need to extract these values
                                // first.
                                std::vector<Type> subStateRewards(b.size());
                                storm::utility::vector::selectVectorValuesRepeatedly(subStateRewards, maybeStates, this->getModel().getNondeterministicChoiceIndices(), this->getModel().getStateRewardVector());
                                storm::utility::vector::addVectorsInPlace(b, subStateRewards);
                            }
                        } else {
                            // If only a state-based reward model is  available, we take this vector as the
                            // right-hand side. As the state reward vector contains entries not just for the
                            // states that we still consider (i.e. maybeStates), we need to extract these values
                            // first.
                            storm::utility::vector::selectVectorValuesRepeatedly(b, maybeStates, this->getModel().getNondeterministicChoiceIndices(), this->getModel().getStateRewardVector());
                        }
                        
                        // Create vector for results for maybe states.
                        std::vector<Type> x(maybeStates.getNumberOfSetBits());
                        
                        // Solve the corresponding system of equations.
                        this->nondeterministicLinearEquationSolver->solveEquationSystem(minimize, submatrix, x, b, subNondeterministicChoiceIndices);
                        
                        // Set values of resulting vector according to result.
                        storm::utility::vector::setVectorValues<Type>(result, maybeStates, x);
                    }
                    
                    // Set values of resulting vector that are known exactly.
                    storm::utility::vector::setVectorValues(result, targetStates, storm::utility::constGetZero<Type>());
                    storm::utility::vector::setVectorValues(result, infinityStates, storm::utility::constGetInfinity<Type>());
                    
                    // Finally, compute a scheduler that achieves the extramal value.
                    storm::storage::TotalScheduler scheduler = computeExtremalScheduler(this->minimumOperatorStack.top(), this->getModel().getTransitionMatrix(), this->getModel().getNondeterministicChoiceIndices(), result, this->getModel().hasStateRewards() ? &this->getModel().getStateRewardVector() : nullptr, this->getModel().hasTransitionRewards() ? &this->getModel().getTransitionRewardMatrix() : nullptr);

                    return std::make_pair(result, scheduler);
                }
                
            protected:
                /*!
                 * Computes the vector of choices that need to be made to minimize/maximize the model checking result for each state.
                 *
                 * @param minimize If set, all choices are resolved such that the solution value becomes minimal and maximal otherwise.
                 * @param nondeterministicResult The model checking result for nondeterministic choices of all states.
                 * @param takenChoices The output vector that is to store the taken choices.
                 * @param nondeterministicChoiceIndices The assignment of states to their nondeterministic choices in the matrix.
                 */
                static storm::storage::TotalScheduler computeExtremalScheduler(bool minimize, storm::storage::SparseMatrix<Type> const& transitionMatrix, std::vector<uint_fast64_t> const& nondeterministicChoiceIndices, std::vector<Type> const& result, std::vector<Type> const* stateRewardVector = nullptr, storm::storage::SparseMatrix<Type> const* transitionRewardMatrix = nullptr) {
                    std::vector<Type> temporaryResult(nondeterministicChoiceIndices.size() - 1);
                    std::vector<Type> nondeterministicResult(result);
                    storm::solver::GmmxxLinearEquationSolver<Type> solver;
                    solver.performMatrixVectorMultiplication(transitionMatrix, nondeterministicResult, nullptr, 1);
                    if (stateRewardVector != nullptr || transitionRewardMatrix != nullptr) {
                        std::vector<Type> totalRewardVector;
                        if (transitionRewardMatrix != nullptr) {
                            totalRewardVector = transitionMatrix.getPointwiseProductRowSumVector(*transitionRewardMatrix);
                            if (stateRewardVector != nullptr) {
                                std::vector<Type> stateRewards(totalRewardVector.size());
                                storm::utility::vector::selectVectorValuesRepeatedly(stateRewards, storm::storage::BitVector(stateRewardVector->size(), true), nondeterministicChoiceIndices, *stateRewardVector);
                                storm::utility::vector::addVectorsInPlace(totalRewardVector, stateRewards);
                            }
                        } else {
                            totalRewardVector.resize(nondeterministicResult.size());
                            storm::utility::vector::selectVectorValuesRepeatedly(totalRewardVector, storm::storage::BitVector(stateRewardVector->size(), true), nondeterministicChoiceIndices, *stateRewardVector);
                        }
                        storm::utility::vector::addVectorsInPlace(nondeterministicResult, totalRewardVector);
                    }
                    
                    std::vector<uint_fast64_t> choices(result.size());
                    
                    if (minimize) {
                        storm::utility::vector::reduceVectorMin(nondeterministicResult, temporaryResult, nondeterministicChoiceIndices, &choices);
                    } else {
                        storm::utility::vector::reduceVectorMax(nondeterministicResult, temporaryResult, nondeterministicChoiceIndices, &choices);
                    }

                    return storm::storage::TotalScheduler(choices);
                }
                
                /*!
                 * A stack used for storing whether we are currently computing min or max probabilities or rewards, respectively.
                 * The topmost element is true if and only if we are currently computing minimum probabilities or rewards.
                 */
                mutable std::stack<bool> minimumOperatorStack;
                
                /*!
                 * A solver that is used for solving systems of linear equations that are the result of nondeterministic choices.
                 */
                std::shared_ptr<storm::solver::AbstractNondeterministicLinearEquationSolver<Type>> nondeterministicLinearEquationSolver;
            };
        } // namespace prctl
    } // namespace modelchecker
} // namespace storm

#endif /* STORM_MODELCHECKER_PRCTL_SPARSEMDPPRCTLMODELCHECKER_H_ */
