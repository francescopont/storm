#pragma once

#include "JaniType.h"

namespace storm {
    namespace jani {
        class BoundedType : public JaniType {

        public:
         
            enum class BaseType {
                Int,
                Real
            };
            
            BoundedType (BaseType const& type, storm::expressions::Expression const& lowerBound, storm::expressions::Expression const& upperBound);
            virtual bool isBoundedType() const override;
            
            BaseType const& getBaseType() const;
            bool isIntegerType() const;
            bool isRealType() const;
            
            virtual std::string getStringRepresentation() const override;
            virtual void substitute(std::map<storm::expressions::Variable, storm::expressions::Expression> const& substitution) override;
            virtual std::unique_ptr<JaniType> clone() const override;

            void setLowerBound(storm::expressions::Expression const& expression);
            void setUpperBound(storm::expressions::Expression const& expression);
            bool hasLowerBound() const;
            bool hasUpperBound() const;
            
            storm::expressions::Expression const& getLowerBound() const;
            storm::expressions::Expression const& getUpperBound() const;



        private:
            BaseType type;
            storm::expressions::Expression lowerBound;
            storm::expressions::Expression upperBound;
        };
    }
}
