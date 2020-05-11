#pragma once

#include <vector>
#include "storm/storage/BitVector.h"

namespace storm {
    namespace pomdp {
        class WinningRegion {
        public:
            WinningRegion(std::vector<uint64_t> const& observationSizes = {});

            bool update(uint64_t observation, storm::storage::BitVector const& winning);
            bool query(uint64_t observation, storm::storage::BitVector const& currently) const;
            bool isWinning(uint64_t observation, uint64_t offset) const {
                storm::storage::BitVector currently(observationSizes[observation]);
                currently.set(offset);
                return query(observation,currently);
            }

            std::vector<storm::storage::BitVector> const& getWinningSetsPerObservation(uint64_t observation) const;

            void setObservationIsWinning(uint64_t observation);

            bool observationIsWinning(uint64_t observation) const;
            storm::expressions::Expression extensionExpression(uint64_t observation, std::vector<storm::expressions::Expression>& varsForStates) const;


                uint64_t getStorageSize() const;
            uint64_t getNumberOfObservations() const;
            bool empty() const;
            void print() const;

            void storeToFile(std::string const& path) const;
            static WinningRegion loadFromFile(std::string const& path);


        private:
            std::vector<std::vector<storm::storage::BitVector>> winningRegion;
            std::vector<uint64_t> observationSizes;
        };
    }
}