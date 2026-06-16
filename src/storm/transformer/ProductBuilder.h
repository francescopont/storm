#pragma once

#include "storm/models/sparse/StateLabeling.h"
#include "storm/storage/BitVector.h"
#include "storm/storage/SparseMatrix.h"
#include "storm/transformer/ProductModel.h"

#include <deque>
#include <map>
#include <vector>

namespace storm {
namespace transformer {

template<typename Model>
class ProductBuilder {
   public:
    typedef storm::storage::SparseMatrix<typename Model::ValueType> matrix_type;

    template<typename ProductOperator>
    static typename Product<Model>::ptr buildProduct(const matrix_type& originalMatrix, ProductOperator& prodOp,
                                                     const storm::storage::BitVector& statesOfInterest) {
        bool deterministic = originalMatrix.hasTrivialRowGrouping();

        typedef storm::storage::sparse::state_type state_type;
        typedef std::pair<state_type, state_type> product_state_type;

        state_type nextState = 0;
        std::map<product_state_type, state_type> productStateToProductIndex;
        std::vector<product_state_type> productIndexToProductState;
        std::vector<state_type> prodInitial;

        // use deque for todo so that the states are handled in the order
        // of their index in the product model, which is required due to the
        // use of the SparseMatrixBuilder that can only handle linear addNextValue
        // calls
        std::deque<state_type> todo;
        for (state_type s_0 : statesOfInterest) {
            state_type q_0 = prodOp.getInitialState(s_0);

            // std::cout << "Initial: " << s_0 << ", " << q_0 << " = " << nextState << "\n";

            product_state_type s_q(s_0, q_0);
            state_type index = nextState++;
            productStateToProductIndex[s_q] = index;
            productIndexToProductState.push_back(s_q);
            prodInitial.push_back(index);
            todo.push_back(index);
        }

        storm::storage::SparseMatrixBuilder<typename Model::ValueType> builder(0, 0, 0, false, deterministic ? false : true, 0);
        std::size_t curRow = 0;
        while (!todo.empty()) {
            state_type prodIndexFrom = todo.front();
            todo.pop_front();

            product_state_type from = productIndexToProductState.at(prodIndexFrom);
            // std::cout << "Handle " << from.first << "," << from.second << " (prodIndexFrom = " << prodIndexFrom << "):\n";
            if (deterministic) {
                typename matrix_type::const_rows row = originalMatrix.getRow(from.first);
                for (auto const& entry : row) {
                    state_type t = entry.getColumn();
                    state_type p = prodOp.getSuccessor(from.second, t);
                    // std::cout << " p = " << p << "\n";
                    product_state_type t_p(t, p);
                    state_type prodIndexTo;
                    auto it = productStateToProductIndex.find(t_p);
                    if (it == productStateToProductIndex.end()) {
                        prodIndexTo = nextState++;
                        todo.push_back(prodIndexTo);
                        productIndexToProductState.push_back(t_p);
                        productStateToProductIndex[t_p] = prodIndexTo;
                        // std::cout << " Adding " << t_p.first << "," << t_p.second << " as " << prodIndexTo << "\n";
                    } else {
                        prodIndexTo = it->second;
                    }
                    // std::cout << " " << t_p.first << "," << t_p.second << ": to = " << prodIndexTo << "\n";

                    // std::cout << " addNextValue(" << prodIndexFrom << "," << prodIndexTo << "," << entry.getValue() << ")\n";
                    builder.addNextValue(prodIndexFrom, prodIndexTo, entry.getValue());
                }
            } else {
                std::size_t numRows = originalMatrix.getRowGroupSize(from.first);
                builder.newRowGroup(curRow);
                for (std::size_t i = 0; i < numRows; i++) {
                    auto const& row = originalMatrix.getRow(from.first, i);
                    for (auto const& entry : row) {
                        state_type t = entry.getColumn();
                        state_type p = prodOp.getSuccessor(from.second, t);
                        // std::cout << " p = " << p << "\n";
                        product_state_type t_p(t, p);
                        state_type prodIndexTo;
                        auto it = productStateToProductIndex.find(t_p);
                        if (it == productStateToProductIndex.end()) {
                            prodIndexTo = nextState++;
                            todo.push_back(prodIndexTo);
                            productIndexToProductState.push_back(t_p);
                            productStateToProductIndex[t_p] = prodIndexTo;
                            // std::cout << " Adding " << t_p.first << "," << t_p.second << " as " << prodIndexTo << "\n";
                        } else {
                            prodIndexTo = it->second;
                        }
                        // std::cout << " " << t_p.first << "," << t_p.second << ": to = " << prodIndexTo << "\n";

                        // std::cout << " addNextValue(" << prodIndexFrom << "," << prodIndexTo << "," << entry.getValue() << ")\n";
                        builder.addNextValue(curRow, prodIndexTo, entry.getValue());
                    }
                    curRow++;
                }
            }
        }

        state_type numberOfProductStates = nextState;

        Model product(builder.build(), storm::models::sparse::StateLabeling(numberOfProductStates));
        storm::storage::BitVector productStatesOfInterest(product.getNumberOfStates());
        for (auto& s : prodInitial) {
            productStatesOfInterest.set(s);
        }
        std::string prodSoiLabel = product.getStateLabeling().addUniqueLabel("soi", productStatesOfInterest);

        // const storm::models::sparse::StateLabeling& orignalLabels = dtmc->getStateLabeling();
        // for (originalLabels.)

        return typename Product<Model>::ptr(
            new Product<Model>(std::move(product), std::move(prodSoiLabel), std::move(productStateToProductIndex), std::move(productIndexToProductState)));
    }

    typedef storm::storage::sparse::state_type state_type;
    typedef std::pair<state_type, state_type> product_state_type;

    static typename ProductModel<Model>::ptr shrinkAndExportProductModel(const matrix_type& originalMatrix,
                                                                         const storm::storage::BitVector& statesOfInterest,
                                                                         std::vector<product_state_type> const& oldProductIndexToProductState,
                                                                         const storm::storage::BitVector& acceptingStates) {

        bool deterministic = originalMatrix.hasTrivialRowGrouping();
        state_type nextState = 0;
        std::vector<product_state_type> productIndexToProductState;
        std::map<state_type, state_type> stateRemapping;
        std::vector<state_type> prodInitial;

        // use deque for todo so that the states are handled in the order
        // of their index in the product model, which is required due to the
        // use of the SparseMatrixBuilder that can only handle linear addNextValue
        // calls
        std::deque<product_state_type> todo;
        for (state_type oldIndex : statesOfInterest) {
            STORM_LOG_ASSERT(!acceptingStates.get(oldIndex), "An initial state cannot be accepting.");
            product_state_type s_q = oldProductIndexToProductState.at(oldIndex);
            productIndexToProductState.push_back(s_q);
            state_type newIndex = nextState++;
            stateRemapping[oldIndex] = newIndex;
            prodInitial.push_back(newIndex);
            product_state_type t_p(newIndex, oldIndex);
            todo.push_back(t_p);
        }

        // handling the accepting sink state
        state_type oldIndex = acceptingStates.getNextSetIndex(0); // just as a convention
        product_state_type s_q = oldProductIndexToProductState.at(oldIndex);
        productIndexToProductState.push_back(s_q);
        state_type sink_accepting_state = nextState++;
        stateRemapping[oldIndex] = sink_accepting_state;
        product_state_type t_err(sink_accepting_state, oldIndex);
        todo.push_back(t_err);

        storm::storage::SparseMatrixBuilder<typename Model::ValueType> builder(0, 0, 0, false, deterministic ? false : true, 0);
        std::size_t curRow = 0;
        while (!todo.empty()) {
            product_state_type fromTuple = todo.front();
            todo.pop_front();
            state_type from = fromTuple.first;
            state_type oldFrom = fromTuple.second;
            STORM_LOG_ASSERT(!acceptingStates.get(oldFrom) || oldFrom == oldIndex, "This state is accepting.");

            if (deterministic) {
                typename matrix_type::const_rows row = originalMatrix.getRow(oldFrom);
                for (auto const& entry : row) {
                    state_type oldTo = entry.getColumn();

                    auto it = stateRemapping.find(oldTo);
                    if (it == stateRemapping.end() && (!acceptingStates.get(oldTo))) {
                        product_state_type t_p = oldProductIndexToProductState.at(oldTo);
                        productIndexToProductState.push_back(t_p);
                        state_type to = nextState++;
                        stateRemapping[oldTo] = to;
                        product_state_type t_tomc(to, oldTo);
                        todo.push_back(t_tomc);
                        builder.addNextValue(from, to, entry.getValue());

                    } else if (!acceptingStates.get(oldTo)) {
                        state_type to = it->second;
                        builder.addNextValue(from, to,entry.getValue());
                    } else {
                        builder.addNextValue(from, sink_accepting_state,entry.getValue());
                    }
                }
            } else {
                std::size_t numRows = originalMatrix.getRowGroupSize(oldFrom);
                builder.newRowGroup(curRow);
                for (std::size_t i = 0; i < numRows; i++) {
                    auto const& row = originalMatrix.getRow(oldFrom, i);
                    for (auto const& entry : row) {
                        state_type oldTo = entry.getColumn();
                        auto it = stateRemapping.find(oldTo);
                        if (it == stateRemapping.end() && (!acceptingStates.get(oldTo))) {
                            //std::cout << "Pushing back state: ";
                            //std::cout << oldTo << '\n';
                            product_state_type t_p = oldProductIndexToProductState.at(oldTo);
                            productIndexToProductState.push_back(t_p);
                            state_type to = nextState++;
                            stateRemapping[oldTo] = to;
                            product_state_type t_tomdp(to, oldTo);
                            todo.push_back(t_tomdp);
                            builder.addNextValue(curRow, to, entry.getValue());
                        } else if (!acceptingStates.get(oldTo)) {
                            state_type to = it->second;
                            builder.addNextValue(curRow, to, entry.getValue());
                        } else {
                            builder.addNextValue(curRow, sink_accepting_state, entry.getValue());
                        }
                    }
                    curRow++;
                }
            }
        }
        state_type numberOfProductStates = nextState;
        
        Model shrinked_product(builder.build(), storm::models::sparse::StateLabeling(numberOfProductStates));
        storm::storage::BitVector productStatesOfInterest(shrinked_product.getNumberOfStates());
        for (auto& s : prodInitial) {
            productStatesOfInterest.set(s);
        }
        std::string prodSoiLabel = shrinked_product.getStateLabeling().addUniqueLabel("soi", productStatesOfInterest);

        return typename ProductModel<Model>::ptr(
            new ProductModel<Model>(std::move(shrinked_product), std::move(productIndexToProductState), sink_accepting_state));
    }
};
}  // namespace transformer
}  // namespace storm
