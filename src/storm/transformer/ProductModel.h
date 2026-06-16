#pragma once

#include <memory>

#include "storm/storage/sparse/StateType.h"
#include "storm/storage/BitVector.h"
#include "storm/transformer/Product.h"

namespace storm {
namespace transformer{
template<typename Model>
class ProductModel {
   public:

    typedef std::shared_ptr<ProductModel<Model>> ptr;

    typedef storm::storage::sparse::state_type state_type;
    typedef std::pair<state_type, state_type> product_state_type;
    typedef std::vector<product_state_type> product_index_to_product_state_vector;

    ProductModel(Model&& productModel,
                 product_index_to_product_state_vector&& productIndexToProductState, state_type acceptingState)
        : productModel(productModel),
          productIndexToProductState(productIndexToProductState),
          acceptingState(acceptingState) {}

    ProductModel(ProductModel<Model>&& product) = default;
    ProductModel& operator=(ProductModel<Model>&& product) = default;

    Model& getModel(){
        return productModel;
    }

    product_index_to_product_state_vector& getProductIndexToProductState(){
        return productIndexToProductState;
    }

    state_type& getAcceptingState(){
        return acceptingState;
    }

    void printMapping(std::ostream& out) const {
        out << "Mapping index -> product state\n";
        for (std::size_t i = 0; i < productIndexToProductState.size(); i++) {
            out << " " << i << ": " << productIndexToProductState.at(i).first << "," << productIndexToProductState.at(i).second << "\n";
        }
    }

   private:
    Model productModel;
    product_index_to_product_state_vector productIndexToProductState;
    state_type acceptingState;

};
}  // namespace transformer
}  // namespace storm
