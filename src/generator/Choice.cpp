#include "src/generator/Choice.h"

#include "src/utility/constants.h"

namespace storm {
    namespace generator {
        
        template<typename ValueType, typename StateType>
        Choice<ValueType, StateType>::Choice(uint_fast64_t actionIndex, bool markovian) : markovian(markovian), actionIndex(actionIndex), distribution(), totalMass(storm::utility::zero<ValueType>()), choiceReward(storm::utility::zero<ValueType>()) {
            // Intentionally left empty.
        }
        
        template<typename ValueType, typename StateType>
        typename storm::storage::Distribution<ValueType, StateType>::iterator Choice<ValueType, StateType>::begin() {
            return distribution.begin();
        }
        
        template<typename ValueType, typename StateType>
        typename storm::storage::Distribution<ValueType, StateType>::const_iterator Choice<ValueType, StateType>::begin() const {
            return distribution.cbegin();
        }
        
        template<typename ValueType, typename StateType>
        typename storm::storage::Distribution<ValueType, StateType>::iterator Choice<ValueType, StateType>::end() {
            return distribution.end();
        }
        
        template<typename ValueType, typename StateType>
        typename storm::storage::Distribution<ValueType, StateType>::const_iterator Choice<ValueType, StateType>::end() const {
            return distribution.cend();
        }
        
        template<typename ValueType, typename StateType>
        void Choice<ValueType, StateType>::addChoiceLabel(uint_fast64_t label) {
            if (!choiceLabels) {
                choiceLabels = LabelSet();
            }
            choiceLabels->insert(label);
        }
        
        template<typename ValueType, typename StateType>
        void Choice<ValueType, StateType>::addChoiceLabels(LabelSet const& labelSet) {
            if (!choiceLabels) {
                choiceLabels = LabelSet();
            }
            choiceLabels->insert(labelSet.begin(), labelSet.end());
        }
        
        template<typename ValueType, typename StateType>
        boost::container::flat_set<uint_fast64_t> const& Choice<ValueType, StateType>::getChoiceLabels() const {
            return *choiceLabels;
        }
        
        template<typename ValueType, typename StateType>
        uint_fast64_t Choice<ValueType, StateType>::getActionIndex() const {
            return actionIndex;
        }
        
        template<typename ValueType, typename StateType>
        ValueType Choice<ValueType, StateType>::getTotalMass() const {
            return totalMass;
        }
        
        template<typename ValueType, typename StateType>
        ValueType& Choice<ValueType, StateType>::getOrAddEntry(StateType const& state) {
            auto stateProbabilityPair = distribution.find(state);
            
            if (stateProbabilityPair == distribution.end()) {
                distribution[state] = ValueType();
            }
            return distribution.at(state);
        }
        
        template<typename ValueType, typename StateType>
        ValueType const& Choice<ValueType, StateType>::getOrAddEntry(StateType const& state) const {
            auto stateProbabilityPair = distribution.find(state);
            
            if (stateProbabilityPair == distribution.end()) {
                distribution[state] = ValueType();
            }
            return distribution.at(state);
        }
        
        template<typename ValueType, typename StateType>
        void Choice<ValueType, StateType>::addProbability(StateType const& state, ValueType const& value) {
            totalMass += value;
            distribution[state] += value;
        }
        
        template<typename ValueType, typename StateType>
        void Choice<ValueType, StateType>::addChoiceReward(ValueType const& value) {
            choiceReward += value;
        }
        
        template<typename ValueType, typename StateType>
        std::size_t Choice<ValueType, StateType>::size() const {
            return distribution.size();
        }
        
        template<typename ValueType, typename StateType>
        std::ostream& operator<<(std::ostream& out, Choice<ValueType, StateType> const& choice) {
            out << "<";
            for (auto const& stateProbabilityPair : choice) {
                out << stateProbabilityPair.first << " : " << stateProbabilityPair.second << ", ";
            }
            out << ">";
            return out;
        }
        
    }
}