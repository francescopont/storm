#include "src/modelchecker/AbstractModelChecker.h"

#include "src/modelchecker/results/QualitativeCheckResult.h"
#include "src/modelchecker/results/QuantitativeCheckResult.h"
#include "src/utility/constants.h"
#include "src/utility/macros.h"
#include "src/exceptions/NotImplementedException.h"
#include "src/exceptions/InvalidOperationException.h"
#include "src/exceptions/InvalidArgumentException.h"
#include "src/exceptions/InternalTypeErrorException.h"

namespace storm {
    namespace modelchecker {
        std::unique_ptr<CheckResult> AbstractModelChecker::check(CheckTask<storm::logic::Formula> const& checkTask) {
            storm::logic::Formula const& formula = checkTask.getFormula();
            STORM_LOG_THROW(this->canHandle(formula), storm::exceptions::InvalidArgumentException, "The model checker is not able to check the formula '" << formula << "'.");
            if (formula.isStateFormula()) {
                return this->checkStateFormula(checkTask.substituteFormula(formula.asStateFormula()));
            } else if (formula.isPathFormula()) {
                if (checkTask.computeProbabilities()) {
                    return this->computeProbabilities(checkTask.substituteFormula(formula.asPathFormula()));
                } else if (checkTask.computeRewards()) {
                    return this->computeRewards(checkTask.substituteFormula(formula.asPathFormula()));
                }
            }
            STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "The given formula '" << formula << "' is invalid.");
        }
        
        std::unique_ptr<CheckResult> AbstractModelChecker::computeProbabilities(CheckTask<storm::logic::PathFormula> const& checkTask) {
            storm::logic::PathFormula const& pathFormula = checkTask.getFormula();
            if (pathFormula.isBoundedUntilFormula()) {
                return this->computeBoundedUntilProbabilities(checkTask.substituteFormula(pathFormula.asBoundedUntilFormula()));
            } else if (pathFormula.isConditionalPathFormula()) {
                return this->computeConditionalProbabilities(checkTask.substituteFormula(pathFormula.asConditionalPathFormula()));
            } else if (pathFormula.isEventuallyFormula()) {
                return this->computeEventuallyProbabilities(checkTask.substituteFormula(pathFormula.asEventuallyFormula()));
            } else if (pathFormula.isGloballyFormula()) {
                return this->computeGloballyProbabilities(checkTask.substituteFormula(pathFormula.asGloballyFormula()));
            } else if (pathFormula.isUntilFormula()) {
                return this->computeUntilProbabilities(checkTask.substituteFormula(pathFormula.asUntilFormula()));
            } else if (pathFormula.isNextFormula()) {
                return this->computeNextProbabilities(checkTask.substituteFormula(pathFormula.asNextFormula()));
            }
            STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "The given formula '" << pathFormula << "' is invalid.");
        }
        
        std::unique_ptr<CheckResult> AbstractModelChecker::computeBoundedUntilProbabilities(CheckTask<storm::logic::BoundedUntilFormula> const& checkTask) {
            STORM_LOG_THROW(false, storm::exceptions::NotImplementedException, "This model checker does not support the formula: " << checkTask.getFormula() << ".");
        }
        
        std::unique_ptr<CheckResult> AbstractModelChecker::computeConditionalProbabilities(CheckTask<storm::logic::ConditionalPathFormula> const& checkTask) {
            STORM_LOG_THROW(false, storm::exceptions::NotImplementedException, "This model checker does not support the formula: " << checkTask.getFormula() << ".");
        }
        
        std::unique_ptr<CheckResult> AbstractModelChecker::computeEventuallyProbabilities(CheckTask<storm::logic::EventuallyFormula> const& checkTask) {
            storm::logic::EventuallyFormula const& pathFormula = checkTask.getFormula();
            storm::logic::UntilFormula newFormula(storm::logic::Formula::getTrueFormula(), pathFormula.getSubformula().asSharedPointer());
            return this->computeUntilProbabilities(checkTask.substituteFormula(newFormula));
        }
        
        std::unique_ptr<CheckResult> AbstractModelChecker::computeGloballyProbabilities(CheckTask<storm::logic::GloballyFormula> const& checkTask) {
            STORM_LOG_THROW(false, storm::exceptions::NotImplementedException, "This model checker does not support the formula: " << checkTask.getFormula() << ".");
        }
        
        std::unique_ptr<CheckResult> AbstractModelChecker::computeNextProbabilities(CheckTask<storm::logic::NextFormula> const& checkTask) {
            STORM_LOG_THROW(false, storm::exceptions::NotImplementedException, "This model checker does not support the formula: " << checkTask.getFormula() << ".");
        }
        
        std::unique_ptr<CheckResult> AbstractModelChecker::computeUntilProbabilities(CheckTask<storm::logic::UntilFormula> const& checkTask) {
            STORM_LOG_THROW(false, storm::exceptions::NotImplementedException, "This model checker does not support the formula: " << checkTask.getFormula() << ".");
        }
        
        std::unique_ptr<CheckResult> AbstractModelChecker::computeRewards(CheckTask<storm::logic::PathFormula> const& checkTask) {
            storm::logic::PathFormula const& rewardPathFormula = checkTask.getFormula();
            if (rewardPathFormula.isCumulativeRewardFormula()) {
                return this->computeCumulativeRewards(checkTask.substituteFormula(rewardPathFormula.asCumulativeRewardFormula()));
            } else if (rewardPathFormula.isInstantaneousRewardFormula()) {
                return this->computeInstantaneousRewards(checkTask.substituteFormula(rewardPathFormula.asInstantaneousRewardFormula()));
            } else if (rewardPathFormula.isEventuallyFormula()) {
                return this->computeReachabilityRewards(checkTask.substituteFormula(rewardPathFormula.asEventuallyFormula()));
            } else if (rewardPathFormula.isLongRunAverageRewardFormula()) {
                return this->computeLongRunAverageRewards(checkTask.substituteFormula(rewardPathFormula.asLongRunAverageRewardFormula()));
            }
            STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "The given formula '" << rewardPathFormula << "' is invalid.");
        }
        
        std::unique_ptr<CheckResult> AbstractModelChecker::computeCumulativeRewards(CheckTask<storm::logic::CumulativeRewardFormula> const& checkTask) {
            STORM_LOG_THROW(false, storm::exceptions::NotImplementedException, "This model checker does not support the formula: " << checkTask.getFormula() << ".");
        }
        
        std::unique_ptr<CheckResult> AbstractModelChecker::computeInstantaneousRewards(CheckTask<storm::logic::InstantaneousRewardFormula> const& checkTask) {
            STORM_LOG_THROW(false, storm::exceptions::NotImplementedException, "This model checker does not support the formula: " << checkTask.getFormula() << ".");
        }
        
        std::unique_ptr<CheckResult> AbstractModelChecker::computeReachabilityRewards(CheckTask<storm::logic::EventuallyFormula> const& checkTask) {
            STORM_LOG_THROW(false, storm::exceptions::NotImplementedException, "This model checker does not support the formula: " << checkTask.getFormula() << ".");
        }
        
        std::unique_ptr<CheckResult> AbstractModelChecker::computeLongRunAverageRewards(CheckTask<storm::logic::LongRunAverageRewardFormula> const& checkTask) {
            STORM_LOG_THROW(false, storm::exceptions::NotImplementedException, "This model checker does not support the formula: " << checkTask.getFormula() << ".");
        }
        
        std::unique_ptr<CheckResult> AbstractModelChecker::computeLongRunAverageProbabilities(CheckTask<storm::logic::StateFormula> const& checkTask) {
            STORM_LOG_THROW(false, storm::exceptions::NotImplementedException, "This model checker does not support the formula: " << checkTask.getFormula() << ".");
        }
        
        std::unique_ptr<CheckResult> AbstractModelChecker::computeExpectedTimes(CheckTask<storm::logic::EventuallyFormula> const& checkTask) {
            STORM_LOG_THROW(false, storm::exceptions::NotImplementedException, "This model checker does not support the formula: " << checkTask.getFormula() << ".");
        }
        
        std::unique_ptr<CheckResult> AbstractModelChecker::checkStateFormula(CheckTask<storm::logic::StateFormula> const& checkTask) {
            storm::logic::StateFormula const& stateFormula = checkTask.getFormula();
            if (stateFormula.isBinaryBooleanStateFormula()) {
                return this->checkBinaryBooleanStateFormula(checkTask.substituteFormula(stateFormula.asBinaryBooleanStateFormula()));
            } else if (stateFormula.isUnaryBooleanStateFormula()) {
                return this->checkUnaryBooleanStateFormula(checkTask.substituteFormula(stateFormula.asUnaryBooleanStateFormula()));
            } else if (stateFormula.isBooleanLiteralFormula()) {
                return this->checkBooleanLiteralFormula(checkTask.substituteFormula(stateFormula.asBooleanLiteralFormula()));
            } else if (stateFormula.isProbabilityOperatorFormula()) {
                return this->checkProbabilityOperatorFormula(checkTask.substituteFormula(stateFormula.asProbabilityOperatorFormula()));
            } else if (stateFormula.isRewardOperatorFormula()) {
                return this->checkRewardOperatorFormula(checkTask.substituteFormula(stateFormula.asRewardOperatorFormula()));
            } else if (stateFormula.isExpectedTimeOperatorFormula()) {
                return this->checkExpectedTimeOperatorFormula(checkTask.substituteFormula(stateFormula.asExpectedTimeOperatorFormula()));
            } else if (stateFormula.isLongRunAverageOperatorFormula()) {
                return this->checkLongRunAverageOperatorFormula(checkTask.substituteFormula(stateFormula.asLongRunAverageOperatorFormula()));
            } else if (stateFormula.isAtomicExpressionFormula()) {
                return this->checkAtomicExpressionFormula(checkTask.substituteFormula(stateFormula.asAtomicExpressionFormula()));
            } else if (stateFormula.isAtomicLabelFormula()) {
                return this->checkAtomicLabelFormula(checkTask.substituteFormula(stateFormula.asAtomicLabelFormula()));
            } else if (stateFormula.isBooleanLiteralFormula()) {
                return this->checkBooleanLiteralFormula(checkTask.substituteFormula(stateFormula.asBooleanLiteralFormula()));
            }
            STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "The given formula '" << stateFormula << "' is invalid.");
        }
        
        std::unique_ptr<CheckResult> AbstractModelChecker::checkAtomicExpressionFormula(CheckTask<storm::logic::AtomicExpressionFormula> const& checkTask) {
            storm::logic::AtomicExpressionFormula const& stateFormula = checkTask.getFormula();
            std::stringstream stream;
            stream << stateFormula.getExpression();
            return this->checkAtomicLabelFormula(checkTask.substituteFormula(storm::logic::AtomicLabelFormula(stream.str())));
        }
        
        std::unique_ptr<CheckResult> AbstractModelChecker::checkAtomicLabelFormula(CheckTask<storm::logic::AtomicLabelFormula> const& checkTask) {
            STORM_LOG_THROW(false, storm::exceptions::NotImplementedException, "This model checker does not support the formula: " << checkTask.getFormula() << ".");
        }
        
        std::unique_ptr<CheckResult> AbstractModelChecker::checkBinaryBooleanStateFormula(CheckTask<storm::logic::BinaryBooleanStateFormula> const& checkTask) {
            storm::logic::BinaryBooleanStateFormula const& stateFormula = checkTask.getFormula();
            STORM_LOG_THROW(stateFormula.getLeftSubformula().isStateFormula() && stateFormula.getRightSubformula().isStateFormula(), storm::exceptions::InvalidArgumentException, "The given formula is invalid.");
            
            std::unique_ptr<CheckResult> leftResult = this->check(checkTask.substituteFormula<storm::logic::Formula>(stateFormula.getLeftSubformula().asStateFormula()));
            std::unique_ptr<CheckResult> rightResult = this->check(checkTask.substituteFormula<storm::logic::Formula>(stateFormula.getRightSubformula().asStateFormula()));
            
            STORM_LOG_THROW(leftResult->isQualitative() && rightResult->isQualitative(), storm::exceptions::InternalTypeErrorException, "Expected qualitative results.");
            
            if (stateFormula.isAnd()) {
                leftResult->asQualitativeCheckResult() &= rightResult->asQualitativeCheckResult();
            } else if (stateFormula.isOr()) {
                leftResult->asQualitativeCheckResult() |= rightResult->asQualitativeCheckResult();
            } else {
                STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "The given formula '" << stateFormula << "' is invalid.");
            }
            
            return leftResult;
        }
        
        std::unique_ptr<CheckResult> AbstractModelChecker::checkBooleanLiteralFormula(CheckTask<storm::logic::BooleanLiteralFormula> const& checkTask) {
            STORM_LOG_THROW(false, storm::exceptions::NotImplementedException, "This model checker does not support the formula: " << checkTask.getFormula() << ".");
        }
        
        std::unique_ptr<CheckResult> AbstractModelChecker::checkProbabilityOperatorFormula(CheckTask<storm::logic::ProbabilityOperatorFormula> const& checkTask) {
            storm::logic::ProbabilityOperatorFormula const& stateFormula = checkTask.getFormula();
            STORM_LOG_THROW(stateFormula.getSubformula().isValidProbabilityPathFormula(), storm::exceptions::InvalidArgumentException, "The given formula is invalid.");
            
            std::unique_ptr<CheckResult> result = this->computeProbabilities(checkTask.substituteFormula(stateFormula.getSubformula().asPathFormula()));
            
            if (stateFormula.hasBound()) {
                STORM_LOG_THROW(result->isQuantitative(), storm::exceptions::InvalidOperationException, "Unable to perform comparison operation on non-quantitative result.");
                return result->asQuantitativeCheckResult().compareAgainstBound(stateFormula.getComparisonType(), stateFormula.getThreshold());
            } else {
                return result;
            }
        }
        
        std::unique_ptr<CheckResult> AbstractModelChecker::checkRewardOperatorFormula(CheckTask<storm::logic::RewardOperatorFormula> const& checkTask) {
            storm::logic::RewardOperatorFormula const& stateFormula = checkTask.getFormula();
            STORM_LOG_THROW(stateFormula.getSubformula().isValidRewardPathFormula(), storm::exceptions::InvalidArgumentException, "The given formula is invalid.");
            
            std::unique_ptr<CheckResult> result = this->computeRewards(checkTask.substituteFormula(stateFormula.getSubformula().asPathFormula()));
            
            if (checkTask.isBoundSet()) {
                STORM_LOG_THROW(result->isQuantitative(), storm::exceptions::InvalidOperationException, "Unable to perform comparison operation on non-quantitative result.");
                return result->asQuantitativeCheckResult().compareAgainstBound(checkTask.getBoundComparisonType(), checkTask.getBoundThreshold());
            } else {
                return result;
            }
        }
        
		std::unique_ptr<CheckResult> AbstractModelChecker::checkExpectedTimeOperatorFormula(CheckTask<storm::logic::ExpectedTimeOperatorFormula> const& checkTask) {
            storm::logic::ExpectedTimeOperatorFormula const& stateFormula = checkTask.getFormula();
			STORM_LOG_THROW(stateFormula.getSubformula().isEventuallyFormula(), storm::exceptions::InvalidArgumentException, "The given formula is invalid.");
            
            std::unique_ptr<CheckResult> result = this->computeExpectedTimes(checkTask.substituteFormula(stateFormula.getSubformula().asEventuallyFormula()));
            
            if (checkTask.isBoundSet()) {
                STORM_LOG_THROW(result->isQuantitative(), storm::exceptions::InvalidOperationException, "Unable to perform comparison operation on non-quantitative result.");
                return result->asQuantitativeCheckResult().compareAgainstBound(checkTask.getBoundComparisonType(), checkTask.getBoundThreshold());
            } else {
                return result;
            }
        }
        
		std::unique_ptr<CheckResult> AbstractModelChecker::checkLongRunAverageOperatorFormula(CheckTask<storm::logic::LongRunAverageOperatorFormula> const& checkTask) {
            storm::logic::LongRunAverageOperatorFormula const& stateFormula = checkTask.getFormula();
			STORM_LOG_THROW(stateFormula.getSubformula().isStateFormula(), storm::exceptions::InvalidArgumentException, "The given formula is invalid.");
            
            std::unique_ptr<CheckResult> result = this->computeLongRunAverageProbabilities(checkTask.substituteFormula(stateFormula.getSubformula().asStateFormula()));
            
            if (checkTask.isBoundSet()) {
                STORM_LOG_THROW(result->isQuantitative(), storm::exceptions::InvalidOperationException, "Unable to perform comparison operation on non-quantitative result.");
                return result->asQuantitativeCheckResult().compareAgainstBound(checkTask.getBoundComparisonType(), checkTask.getBoundThreshold());
            } else {
                return result;
            }
        }
        
        std::unique_ptr<CheckResult> AbstractModelChecker::checkUnaryBooleanStateFormula(CheckTask<storm::logic::UnaryBooleanStateFormula> const& checkTask) {
            storm::logic::UnaryBooleanStateFormula const& stateFormula = checkTask.getFormula();
            std::unique_ptr<CheckResult> subResult = this->check(checkTask.substituteFormula<storm::logic::Formula>(stateFormula.getSubformula()));
            STORM_LOG_THROW(subResult->isQualitative(), storm::exceptions::InternalTypeErrorException, "Expected qualitative result.");
            if (stateFormula.isNot()) {
                subResult->asQualitativeCheckResult().complement();
            } else {
                STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "The given formula '" << stateFormula << "' is invalid.");
            }
            return subResult;
        }
    }
}