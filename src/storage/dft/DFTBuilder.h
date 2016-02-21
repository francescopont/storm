
#ifndef DFTBUILDER_H
#define	DFTBUILDER_H

#include "DFTElements.h"
#include "elements/DFTRestriction.h"
#include <iostream>
#include <unordered_map>
#include <map>
#include "src/utility/macros.h"

namespace storm {
    namespace storage {
        template<typename ValueType>
        class DFT;

        template<typename ValueType>
        class DFTBuilder {

            using DFTElementPointer = std::shared_ptr<DFTElement<ValueType>>;
            using DFTElementVector = std::vector<DFTElementPointer>;
            using DFTGatePointer = std::shared_ptr<DFTGate<ValueType>>;
            using DFTGateVector = std::vector<DFTGatePointer>;
            using DFTDependencyPointer = std::shared_ptr<DFTDependency<ValueType>>;
            using DFTRestrictionPointer = std::shared_ptr<DFTRestriction<ValueType>>;

        private:
            std::size_t mNextId = 0;
            std::string topLevelIdentifier;
            std::unordered_map<std::string, DFTElementPointer> mElements;
            std::unordered_map<DFTElementPointer, std::vector<std::string>> mChildNames;
            std::unordered_map<DFTRestrictionPointer, std::vector<std::string>> mRestrictionChildNames;
            std::vector<DFTDependencyPointer> mDependencies;
            std::vector<DFTRestrictionPointer> mRestrictions;
            
        public:
            DFTBuilder() {
                
            }
            
            bool addAndElement(std::string const& name, std::vector<std::string> const& children) {
                return addStandardGate(name, children, DFTElementType::AND);
            }
            
            bool addOrElement(std::string const& name, std::vector<std::string> const& children) {
                return addStandardGate(name, children, DFTElementType::OR);
            }
            
            bool addPandElement(std::string const& name, std::vector<std::string> const& children) {
                return addStandardGate(name, children, DFTElementType::PAND);
            }
            
            bool addPorElement(std::string const& name, std::vector<std::string> const& children) {
                return addStandardGate(name, children, DFTElementType::POR);
            }
            
            bool addSpareElement(std::string const& name, std::vector<std::string> const& children) {
                return addStandardGate(name, children, DFTElementType::SPARE);
            }

            bool addSequenceEnforcer(std::string const& name, std::vector<std::string> const& children) {
                return addRestriction(name, children, DFTElementType::SEQ);
            }

            bool addMutex(std::string const& name, std::vector<std::string> const& children) {
                return addRestriction(name, children, DFTElementType::MUTEX);
            }
            
            bool addDepElement(std::string const& name, std::vector<std::string> const& children, ValueType probability) {
                if(children.size() <= 1) {
                    STORM_LOG_ERROR("Dependencies require at least two children");
                }
                if(mElements.count(name) != 0) {
                    // Element with that name already exists.
                    return false;
                }

                if (storm::utility::isZero(probability)) {
                    // Element is superfluous
                    return true;
                }
                std::string trigger = children[0];

                //TODO Matthias: collect constraints for SMT solving
                //0 <= probability <= 1
                if (!storm::utility::isOne(probability) && children.size() > 2) {
                    //TODO Matthias: introduce additional element for probability and then add pdeps with probability 1 to children
                    STORM_LOG_ERROR("Probability != 1 for more than one child currently not supported.");
                    return false;
                }

                for (size_t i = 1; i < children.size(); ++i) {
                    std::string nameDep = name + "_" + std::to_string(i);
                    if(mElements.count(nameDep) != 0) {
                        // Element with that name already exists.
                        STORM_LOG_ERROR("Element with name: " << nameDep << " already exists.");
                        return false;
                    }
                    assert(storm::utility::isOne(probability) || children.size() == 2);
                    DFTDependencyPointer element = std::make_shared<DFTDependency<ValueType>>(mNextId++, nameDep, trigger, children[i], probability);
                    mElements[element->name()] = element;
                    mDependencies.push_back(element);
                }
                return true;
            }

            bool addVotElement(std::string const& name, unsigned threshold, std::vector<std::string> const& children) {
                assert(children.size() > 0);
                if(mElements.count(name) != 0) {
                    STORM_LOG_ERROR("Element with name: " << name << " already exists.");
                    return false;
                }
                // It is an and-gate
                if(children.size() == threshold) {
                    return addAndElement(name, children);
                }
                // It is an or-gate
                if(threshold == 1) {
                    return addOrElement(name, children);
                }
                
                if(threshold > children.size()) {
                    STORM_LOG_ERROR("Voting gates with threshold higher than the number of children is not supported.");
                    return false;
                }
                DFTElementPointer element = std::make_shared<DFTVot<ValueType>>(mNextId++, name, threshold);
                
                mElements[name] = element;
                mChildNames[element] = children;
                return true;
            }
            
            bool addBasicElement(std::string const& name, ValueType failureRate, ValueType dormancyFactor) {
                //TODO Matthias: collect constraints for SMT solving
                //failureRate > 0
                //0 <= dormancyFactor <= 1

                mElements[name] = std::make_shared<DFTBE<ValueType>>(mNextId++, name, failureRate, dormancyFactor);
                return true;
            }
            
            bool setTopLevel(std::string const& tle) {
                topLevelIdentifier = tle;
                return mElements.count(tle) > 0;
            }
            
            DFT<ValueType> build();

         
        private:
            
            unsigned computeRank(DFTElementPointer const& elem);
            
            bool addStandardGate(std::string const& name, std::vector<std::string> const& children, DFTElementType tp);
            bool addRestriction(std::string const& name, std::vector<std::string> const& children, DFTElementType tp);

            enum class topoSortColour {WHITE, BLACK, GREY}; 
            
            void topoVisit(DFTElementPointer const& n, std::map<DFTElementPointer, topoSortColour, OrderElementsById<ValueType>>& visited, DFTElementVector& L);

            DFTElementVector topoSort();
            
        };
    }
}



#endif	/* DFTBUILDER_H */

