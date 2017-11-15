#include "gtest/gtest.h"
#include "storm-config.h"

#include "test/storm_gtest.h"

#include "storm/parser/FormulaParser.h"
#include "storm/logic/Formulas.h"
#include "storm/models/sparse/Mdp.h"
#include "storm/models/symbolic/Mdp.h"
#include "storm/models/sparse/StandardRewardModel.h"
#include "storm/modelchecker/prctl/SparseMdpPrctlModelChecker.h"
#include "storm/modelchecker/prctl/HybridMdpPrctlModelChecker.h"
#include "storm/modelchecker/prctl/SymbolicMdpPrctlModelChecker.h"
#include "storm/modelchecker/results/QuantitativeCheckResult.h"
#include "storm/modelchecker/results/ExplicitQualitativeCheckResult.h"
#include "storm/modelchecker/results/SymbolicQualitativeCheckResult.h"
#include "storm/modelchecker/results/QualitativeCheckResult.h"
#include "storm/settings/SettingsManager.h"
#include "storm/api/builder.h"
#include "storm/api/model_descriptions.h"
#include "storm/api/properties.h"
#include "storm/environment/solver/MinMaxSolverEnvironment.h"

#include "storm/settings/modules/GeneralSettings.h"
#include "storm/settings/modules/CoreSettings.h"
#include "storm/parser/AutoParser.h"

namespace {
    class SparseDoubleValueIterationEnvironment {
    public:
        static const storm::dd::DdType ddType = storm::dd::DdType::Sylvan; // Unused for sparse models
        static const bool isExact = false;
        typedef double ValueType;
        typedef storm::models::sparse::Mdp<ValueType> ModelType;
        static storm::Environment createEnvironment() {
            storm::Environment env;
            env.solver().minMax().setMethod(storm::solver::MinMaxMethod::ValueIteration);
            env.solver().minMax().setPrecision(storm::utility::convertNumber<storm::RationalNumber>(1e-8));
            return env;
        }
        // TODO: this should be part of the environment
        static storm::settings::modules::CoreSettings::Engine engine() { return storm::settings::modules::CoreSettings::Engine::Sparse; }
    };
    class SparseDoubleSoundValueIterationEnvironment {
    public:
        static const storm::dd::DdType ddType = storm::dd::DdType::Sylvan; // Unused for sparse models
        static const bool isExact = false;
        typedef double ValueType;
        typedef storm::models::sparse::Mdp<ValueType> ModelType;
        static storm::Environment createEnvironment() {
            storm::Environment env;
            env.solver().setForceSoundness(true);
            env.solver().minMax().setMethod(storm::solver::MinMaxMethod::ValueIteration);
            env.solver().minMax().setPrecision(storm::utility::convertNumber<storm::RationalNumber>(1e-6));
            return env;
        }
        // TODO: this should be part of the environment
        static storm::settings::modules::CoreSettings::Engine engine() { return storm::settings::modules::CoreSettings::Engine::Sparse; }
    };
    class SparseRationalPolicyIterationEnvironment {
    public:
        static const storm::dd::DdType ddType = storm::dd::DdType::Sylvan; // Unused for sparse models
        static const bool isExact = true;
        typedef storm::RationalNumber ValueType;
        typedef storm::models::sparse::Mdp<ValueType> ModelType;
        static storm::Environment createEnvironment() {
            storm::Environment env;
            env.solver().minMax().setMethod(storm::solver::MinMaxMethod::PolicyIteration);
            return env;
        }
        // TODO: this should be part of the environment
        static storm::settings::modules::CoreSettings::Engine engine() { return storm::settings::modules::CoreSettings::Engine::Sparse; }
    };
    class SparseRationalRationalSearchEnvironment {
    public:
        static const storm::dd::DdType ddType = storm::dd::DdType::Sylvan; // Unused for sparse models
        static const bool isExact = true;
        typedef storm::RationalNumber ValueType;
        typedef storm::models::sparse::Mdp<ValueType> ModelType;
        static storm::Environment createEnvironment() {
            storm::Environment env;
            env.solver().minMax().setMethod(storm::solver::MinMaxMethod::RationalSearch);
            return env;
        }
        // TODO: this should be part of the environment
        static storm::settings::modules::CoreSettings::Engine engine() { return storm::settings::modules::CoreSettings::Engine::Sparse; }
    };
    
    template<typename TestType>
    class MdpPrctlModelCheckerTest : public ::testing::Test {
    public:
        typedef typename TestType::ValueType ValueType;
        typedef typename storm::models::sparse::Mdp<ValueType> SparseModelType;
        typedef typename storm::models::symbolic::Mdp<TestType::ddType, ValueType> SymbolicModelType;
        
        MdpPrctlModelCheckerTest() : _environment(TestType::createEnvironment()) {}
        storm::Environment const& env() const { return _environment; }
        ValueType parseNumber(std::string const& input) const { return storm::utility::convertNumber<ValueType>(input);}
        ValueType precision() const { return TestType::isExact ? parseNumber("0") : parseNumber("1e-6");}
        bool isSparseModel() const { return std::is_same<typename TestType::ModelType, SparseModelType>::value; }
        bool isSymbolicModel() const { return std::is_same<typename TestType::ModelType, SymbolicModelType>::value; }
        
        template <typename MT = typename TestType::ModelType>
        typename std::enable_if<std::is_same<MT, SparseModelType>::value, std::pair<std::shared_ptr<MT>, std::vector<std::shared_ptr<storm::logic::Formula const>>>>::type
        buildModelFormulas(std::string const& pathToPrismFile, std::string const& formulasAsString, std::string const& constantDefinitionString = "") const {
            std::pair<std::shared_ptr<MT>, std::vector<std::shared_ptr<storm::logic::Formula const>>> result;
            storm::prism::Program program = storm::api::parseProgram(pathToPrismFile);
            program = storm::utility::prism::preprocess(program, constantDefinitionString);
            result.second = storm::api::extractFormulasFromProperties(storm::api::parsePropertiesForPrismProgram(formulasAsString, program));
            result.first = storm::api::buildSparseModel<ValueType>(program, result.second)->template as<MT>();
            return result;
        }
        
        template <typename MT = typename TestType::ModelType>
        typename std::enable_if<std::is_same<MT, SymbolicModelType>::value, std::pair<std::shared_ptr<MT>, std::vector<std::shared_ptr<storm::logic::Formula const>>>>::type
        buildModelFormulas(std::string const& pathToPrismFile, std::string const& formulasAsString, std::string const& constantDefinitionString = "") const {
            std::pair<std::shared_ptr<MT>, std::vector<std::shared_ptr<storm::logic::Formula const>>> result;
            storm::prism::Program program = storm::api::parseProgram(pathToPrismFile);
            program = storm::utility::prism::preprocess(program, constantDefinitionString);
            result.second = storm::api::extractFormulasFromProperties(storm::api::parsePropertiesForPrismProgram(formulasAsString, program));
            result.first = storm::api::buildSymbolicModel<TestType::ddType, ValueType>(program, result.second)->template as<MT>();
            return result;
        }
        
        std::vector<storm::modelchecker::CheckTask<storm::logic::Formula, ValueType>> getTasks(std::vector<std::shared_ptr<storm::logic::Formula const>> const& formulas) const {
            std::vector<storm::modelchecker::CheckTask<storm::logic::Formula, ValueType>> result;
            for (auto const& f : formulas) {
                result.emplace_back(*f);
            }
            return result;
        }
        
        template <typename MT = typename TestType::ModelType>
        typename std::enable_if<std::is_same<MT, SparseModelType>::value, std::shared_ptr<storm::modelchecker::AbstractModelChecker<MT>>>::type
        createModelChecker(std::shared_ptr<MT> const& model) const {
            return std::make_shared<storm::modelchecker::SparseMdpPrctlModelChecker<SparseModelType>>(*model);
        }
        
        template <typename MT = typename TestType::ModelType>
        typename std::enable_if<std::is_same<MT, SymbolicModelType>::value, std::shared_ptr<storm::modelchecker::AbstractModelChecker<MT>>>::type
        createModelChecker(std::shared_ptr<MT> const& model) const {
            if (TestType::engine() == storm::settings::modules::CoreSettings::Engine::Hybrid) {
                return std::make_shared<storm::modelchecker::HybridMdpPrctlModelChecker<SymbolicModelType>>(*model);
            } else if (TestType::engine() == storm::settings::modules::CoreSettings::Engine::Dd) {
                return std::make_shared<storm::modelchecker::SymbolicMdpPrctlModelChecker<SymbolicModelType>>(*model);
            }
        }
        
        bool getQualitativeResultAtInitialState(std::shared_ptr<storm::models::Model<ValueType>> const& model, std::unique_ptr<storm::modelchecker::CheckResult>& result) {
            auto filter = getInitialStateFilter(model);
            result->filter(*filter);
            return result->asQualitativeCheckResult().forallTrue();
        }
        
        ValueType getQuantitativeResultAtInitialState(std::shared_ptr<storm::models::Model<ValueType>> const& model, std::unique_ptr<storm::modelchecker::CheckResult>& result) {
            auto filter = getInitialStateFilter(model);
            result->filter(*filter);
            return result->asQuantitativeCheckResult<ValueType>().getMin();
        }
    
    private:
        storm::Environment _environment;
        
        std::unique_ptr<storm::modelchecker::QualitativeCheckResult> getInitialStateFilter(std::shared_ptr<storm::models::Model<ValueType>> const& model) const {
            if (isSparseModel()) {
                return std::make_unique<storm::modelchecker::ExplicitQualitativeCheckResult>(model->template as<SparseModelType>()->getInitialStates());
            } else {
                return std::make_unique<storm::modelchecker::SymbolicQualitativeCheckResult<TestType::ddType>>(model->template as<SymbolicModelType>()->getReachableStates(), model->template as<SymbolicModelType>()->getInitialStates());
            }
        }
    };
  
    typedef ::testing::Types<
            SparseDoubleValueIterationEnvironment,
            SparseDoubleSoundValueIterationEnvironment,
            SparseRationalPolicyIterationEnvironment,
            SparseRationalRationalSearchEnvironment
    > TestingTypes;
    
    TYPED_TEST_CASE(MdpPrctlModelCheckerTest, TestingTypes);
    

    TYPED_TEST(MdpPrctlModelCheckerTest, Dice) {
        std::string formulasString = "Pmin=? [F \"two\"]";
                 formulasString += "; Pmax=? [F \"two\"]";
                 formulasString += "; Pmin=? [F \"three\"]";
                 formulasString += "; Pmax=? [F \"three\"]";
                 formulasString += "; Pmin=? [F \"four\"]";
                 formulasString += "; Pmax=? [F \"four\"]";
                 formulasString += "; Rmin=? [F \"done\"]";
                 formulasString += "; Rmax=? [F \"done\"]";
        
        auto modelFormulas = this->buildModelFormulas(STORM_TEST_RESOURCES_DIR "/mdp/two_dice.nm", formulasString);
        auto model = std::move(modelFormulas.first);
        auto tasks = this->getTasks(modelFormulas.second);
        EXPECT_EQ(169ul, model->getNumberOfStates());
        EXPECT_EQ(436ul, model->getNumberOfTransitions());
        ASSERT_EQ(model->getType(), storm::models::ModelType::Mdp);
        auto checker = this->createModelChecker(model);
        std::unique_ptr<storm::modelchecker::CheckResult> result;
        
        result = checker->check(this->env(), tasks[0]);
        EXPECT_NEAR(this->parseNumber("1/36"), this->getQuantitativeResultAtInitialState(model, result), this->precision());
        
        result = checker->check(this->env(), tasks[1]);
        EXPECT_NEAR(this->parseNumber("1/36"), this->getQuantitativeResultAtInitialState(model, result), this->precision());
        
        result = checker->check(this->env(), tasks[2]);
        EXPECT_NEAR(this->parseNumber("2/36"), this->getQuantitativeResultAtInitialState(model, result), this->precision());
        
        result = checker->check(this->env(), tasks[3]);
        EXPECT_NEAR(this->parseNumber("2/36"), this->getQuantitativeResultAtInitialState(model, result), this->precision());
        
        result = checker->check(this->env(), tasks[4]);
        EXPECT_NEAR(this->parseNumber("3/36"), this->getQuantitativeResultAtInitialState(model, result), this->precision());
        
        result = checker->check(this->env(), tasks[5]);
        EXPECT_NEAR(this->parseNumber("3/36"), this->getQuantitativeResultAtInitialState(model, result), this->precision());
        
        result = checker->check(this->env(), tasks[6]);
        EXPECT_NEAR(this->parseNumber("22/3"), this->getQuantitativeResultAtInitialState(model, result), this->precision());
        
        result = checker->check(this->env(), tasks[7]);
        EXPECT_NEAR(this->parseNumber("22/3"), this->getQuantitativeResultAtInitialState(model, result), this->precision());
    }
    
    TYPED_TEST(MdpPrctlModelCheckerTest, AsynchronousLeader) {
        std::string formulasString = "Pmin=? [F \"elected\"]";
                 formulasString += "; Pmax=? [F \"elected\"]";
                 formulasString += "; Pmin=? [F<=25 \"elected\"]";
                 formulasString += "; Pmax=? [F<=25 \"elected\"]";
                 formulasString += "; Rmin=? [F \"elected\"]";
                 formulasString += "; Rmax=? [F \"elected\"]";

        auto modelFormulas = this->buildModelFormulas(STORM_TEST_RESOURCES_DIR "/mdp/leader4.nm", formulasString);
        auto model = std::move(modelFormulas.first);
        auto tasks = this->getTasks(modelFormulas.second);
        EXPECT_EQ(3172ul, model->getNumberOfStates());
        EXPECT_EQ(7144ul, model->getNumberOfTransitions());
        ASSERT_EQ(model->getType(), storm::models::ModelType::Mdp);
        auto checker = this->createModelChecker(model);
        std::unique_ptr<storm::modelchecker::CheckResult> result;
   
        
        result = checker->check(this->env(), tasks[0]);
        EXPECT_NEAR(this->parseNumber("1"), this->getQuantitativeResultAtInitialState(model, result), this->precision());
        
        result = checker->check(this->env(), tasks[1]);
        EXPECT_NEAR(this->parseNumber("1"), this->getQuantitativeResultAtInitialState(model, result), this->precision());
        
        result = checker->check(this->env(), tasks[2]);
        EXPECT_NEAR(this->parseNumber("1/16"), this->getQuantitativeResultAtInitialState(model, result), this->precision());
        
        result = checker->check(this->env(), tasks[3]);
        EXPECT_NEAR(this->parseNumber("1/16"), this->getQuantitativeResultAtInitialState(model, result), this->precision());
        
        result = checker->check(this->env(), tasks[4]);
        EXPECT_NEAR(this->parseNumber("30/7"), this->getQuantitativeResultAtInitialState(model, result), this->precision());
        
        result = checker->check(this->env(), tasks[5]);
        EXPECT_NEAR(this->parseNumber("30/7"), this->getQuantitativeResultAtInitialState(model, result), this->precision());
    }
    
    TYPED_TEST(MdpPrctlModelCheckerTest, TinyRewards) {
        std::string formulasString = "Rmin=? [F \"target\"]";
        auto modelFormulas = this->buildModelFormulas(STORM_TEST_RESOURCES_DIR "/mdp/tiny_rewards.nm", formulasString);
        auto model = std::move(modelFormulas.first);
        auto tasks = this->getTasks(modelFormulas.second);
        EXPECT_EQ(3ul, model->getNumberOfStates());
        EXPECT_EQ(4ul, model->getNumberOfTransitions());
        ASSERT_EQ(model->getType(), storm::models::ModelType::Mdp);
        auto checker = this->createModelChecker(model);
        std::unique_ptr<storm::modelchecker::CheckResult> result;

        result = checker->check(this->env(), tasks[0]);
        EXPECT_NEAR(this->parseNumber("1"), this->getQuantitativeResultAtInitialState(model, result), this->precision());
    }
}