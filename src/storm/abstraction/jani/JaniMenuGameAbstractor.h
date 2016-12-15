#pragma once

#include "storm/storage/dd/DdType.h"

#include "storm/abstraction/MenuGameAbstractor.h"
#include "storm/abstraction/MenuGame.h"
#include "storm/abstraction/RefinementCommand.h"
#include "storm/abstraction/ValidBlockAbstractor.h"
#include "storm/abstraction/jani/JaniAbstractionInformation.h"
#include "storm/abstraction/jani/AutomatonAbstractor.h"

#include "storm/storage/dd/Add.h"

#include "storm/storage/expressions/Expression.h"

#include "storm/settings/modules/AbstractionSettings.h"

namespace storm {
    namespace utility {
        namespace solver {
            class SmtSolverFactory;
        }
    }
    
    namespace models {
        namespace symbolic {
            template<storm::dd::DdType Type, typename ValueType>
            class StochasticTwoPlayerGame;
        }
    }
    
    namespace jani {
        // Forward-declare concrete Model class.
        class Model;
    }
    
    namespace abstraction {
        namespace jani {
            
            template <storm::dd::DdType DdType, typename ValueType>
            class JaniMenuGameAbstractor : public MenuGameAbstractor<DdType, ValueType> {
            public:
                /*!
                 * Constructs an abstractor for the given model.
                 *
                 * @param model The concrete model for which to build the abstraction.
                 * @param smtSolverFactory A factory that is to be used for creating new SMT solvers.
                 */
                JaniMenuGameAbstractor(storm::jani::Model const& model, std::shared_ptr<storm::utility::solver::SmtSolverFactory> const& smtSolverFactory);
                
                JaniMenuGameAbstractor(JaniMenuGameAbstractor const&) = default;
                JaniMenuGameAbstractor& operator=(JaniMenuGameAbstractor const&) = default;
                JaniMenuGameAbstractor(JaniMenuGameAbstractor&&) = default;
                JaniMenuGameAbstractor& operator=(JaniMenuGameAbstractor&&) = default;
                
                /*!
                 * Uses the current set of predicates to derive the abstract menu game in the form of an ADD.
                 *
                 * @return The abstract stochastic two player game.
                 */
                MenuGame<DdType, ValueType> abstract() override;
                
                /*!
                 * Retrieves information about the abstraction.
                 *
                 * @return The abstraction information object.
                 */
                AbstractionInformation<DdType> const& getAbstractionInformation() const override;
                
                /*!
                 * Retrieves the guard predicate of the given player 1 choice.
                 *
                 * @param player1Choice The choice for which to retrieve the guard.
                 * @return The guard of the player 1 choice.
                 */
                storm::expressions::Expression const& getGuard(uint64_t player1Choice) const override;
                
                /*!
                 * Retrieves a mapping from variables to expressions that define their updates wrt. to the given player
                 * 1 choice and auxiliary choice.
                 */
                std::map<storm::expressions::Variable, storm::expressions::Expression> getVariableUpdates(uint64_t player1Choice, uint64_t auxiliaryChoice) const override;
                
                /*!
                 * Retrieves the range of player 1 choices.
                 */
                std::pair<uint64_t, uint64_t> getPlayer1ChoiceRange() const override;
                
                /*!
                 * Retrieves the expression that characterizes the initial states.
                 */
                storm::expressions::Expression getInitialExpression() const override;
                
                /*!
                 * Retrieves the set of states (represented by a BDD) satisfying the given predicate, assuming that it
                 * was either given as an initial predicate or used as a refining predicate later.
                 *
                 * @param predicate The predicate for which to retrieve the states.
                 * @return The BDD representing the set of states.
                 */
                storm::dd::Bdd<DdType> getStates(storm::expressions::Expression const& predicate);
                
                /*!
                 * Performs the given refinement command.
                 *
                 * @param command The command to perform.
                 */
                virtual void refine(RefinementCommand const& command) override;

                /*!
                 * Exports the current state of the abstraction in the dot format to the given file.
                 *
                 * @param filename The name of the file to which to write the dot output.
                 * @param highlightStates A BDD characterizing states that will be highlighted.
                 * @param filter A filter that is applied to select which part of the game to export.
                 */
                void exportToDot(std::string const& filename, storm::dd::Bdd<DdType> const& highlightStates, storm::dd::Bdd<DdType> const& filter) const override;
                
            private:                
                /*!
                 * Builds the stochastic game representing the abstraction of the program.
                 *
                 * @return The stochastic game.
                 */
                std::unique_ptr<MenuGame<DdType, ValueType>> buildGame();
                
                /*!
                 * Decodes the given choice over the auxiliary and successor variables to a mapping from update indices
                 * to bit vectors representing the successors under these updates.
                 *
                 * @param choice The choice to decode.
                 */
                std::map<uint_fast64_t, storm::storage::BitVector> decodeChoiceToUpdateSuccessorMapping(storm::dd::Bdd<DdType> const& choice) const;
                
                // The concrete model this abstractor refers to.
                std::reference_wrapper<storm::jani::Model const> model;

                // A factory that can be used to create new SMT solvers.
                std::shared_ptr<storm::utility::solver::SmtSolverFactory> smtSolverFactory;
                
                // An object containing all information about the abstraction like predicates and the corresponding DDs.
                JaniAbstractionInformation<DdType> abstractionInformation;
                
                // The abstract modules of the abstract program.
                std::vector<AutomatonAbstractor<DdType, ValueType>> automata;
                
                // A state-set abstractor used to determine the initial states of the abstraction.
                StateSetAbstractor<DdType, ValueType> initialStateAbstractor;
                
                // An object that is used to compute the valid blocks.
                ValidBlockAbstractor<DdType> validBlockAbstractor;
                                
                // An ADD characterizing the probabilities and source/target locations of edges and their updates.
                storm::dd::Add<DdType, ValueType> edgeDecoratorAdd;
                
                // The current game-based abstraction.
                std::unique_ptr<MenuGame<DdType, ValueType>> currentGame;
                
                // A flag storing whether a refinement was performed.
                bool refinementPerformed;
            };
        }
    }
}