/*
 * BoundedUntil.h
 *
 *  Created on: 27.11.2012
 *      Author: Christian Dehnert
 */

#ifndef STORM_FORMULA_LTL_BOUNDEDEVENTUALLY_H_
#define STORM_FORMULA_LTL_BOUNDEDEVENTUALLY_H_

#include "AbstractLtlFormula.h"
#include "src/formula/AbstractFormulaChecker.h"
#include <cstdint>
#include <string>
#include "src/modelchecker/ltl/ForwardDeclarations.h"

namespace storm {
namespace property {
namespace ltl {

template <class T> class BoundedEventually;

/*!
 *  @brief Interface class for model checkers that support BoundedEventually.
 *   
 *  All model checkers that support the formula class BoundedEventually must inherit
 *  this pure virtual class.
 */
template <class T>
class IBoundedEventuallyModelChecker {
    public:
		/*!
         *  @brief Evaluates BoundedEventually formula within a model checker.
         *
         *  @param obj Formula object with subformulas.
         *  @return Result of the formula for every node.
         */
        virtual std::vector<T> checkBoundedEventually(const BoundedEventually<T>& obj) const = 0;
};

/*!
 *	@brief Interface class for visitors that support BoundedEventually.
 *
 *	All visitors that support the formula class BoundedEventually must inherit
 *	this pure virtual class.
 */
template <class T>
class IBoundedEventuallyVisitor {
	public:
		/*!
		 *	@brief Evaluates BoundedEventually formula within a model checker.
		 *
		 *	@param obj Formula object with subformulas.
		 *	@return Result of the formula for every node.
		 */
		virtual void visitBoundedEventually(const BoundedEventually<T>& obj) = 0;
};

/*!
 * @brief
 * Class for an abstract (path) formula tree with a BoundedEventually node as root.
 *
 * Has one Abstract LTL formulas as sub formula/tree.
 *
 * @par Semantics
 * The formula holds iff in at most \e bound steps, formula \e child holds.
 *
 * The subtrees are seen as part of the object and deleted with the object
 * (this behavior can be prevented by setting them to NULL before deletion)
 *
 * @see AbstractLtlFormula
 */
template <class T>
class BoundedEventually : public AbstractLtlFormula<T> {

public:

	/*!
	 * Empty constructor
	 */
	BoundedEventually() {
		this->child = nullptr;
		bound = 0;
	}

	/*!
	 * Constructor
	 *
	 * @param child The child formula subtree
	 * @param bound The maximal number of steps
	 */
	BoundedEventually(AbstractLtlFormula<T>* child, uint_fast64_t bound) {
		this->child = child;
		this->bound = bound;
	}

	/*!
	 * Destructor.
	 *
	 * Also deletes the subtrees.
	 * (this behaviour can be prevented by setting the subtrees to NULL before deletion)
	 */
	virtual ~BoundedEventually() {
	  if (child != nullptr) {
		  delete child;
	  }
	}


	/*!
	 * Clones the called object.
	 *
	 * Performs a "deep copy", i.e. the subtrees of the new object are clones of the original ones
	 *
	 * @returns a new BoundedUntil-object that is identical the called object.
	 */
	virtual AbstractLtlFormula<T>* clone() const override {
		BoundedEventually<T>* result = new BoundedEventually<T>();
		result->setBound(this->getBound());
		if (this->childIsSet()) {
			result->setChild(this->getChild().clone());
		}
		return result;
	}


	/*!
	 * Calls the model checker to check this formula.
	 * Needed to infer the correct type of formula class.
	 *
	 * @note This function should only be called in a generic check function of a model checker class. For other uses,
	 *       the methods of the model checker should be used.
	 *
	 * @returns A vector indicating the probability that the formula holds for each state.
	 */
	virtual std::vector<T> check(const storm::modelchecker::ltl::AbstractModelChecker<T>& modelChecker) const override {
		return modelChecker.template as<IBoundedEventuallyModelChecker>()->checkBoundedEventually(*this);
	}

	virtual void visit(visitor::AbstractLtlFormulaVisitor<T>& visitor) const override {
		visitor.template as<IBoundedEventuallyVisitor>()->visitBoundedEventually(*this);
	}

	/*!
	 * @returns a string representation of the formula
	 */
	virtual std::string toString() const override {
		std::string result = "F<=";
		result += std::to_string(bound);
		result += " ";
		result += child->toString();
		return result;
	}

	/*!
     *  @brief Checks if the subtree conforms to some logic.
     *
     *  @param checker Formula checker object.
     *  @return true iff the subtree conforms to some logic.
     */
	virtual bool validate(const AbstractFormulaChecker<T>& checker) const override {
		return checker.validate(this->child);
	}

	/*!
	 * @returns the child node
	 */
	const AbstractLtlFormula<T>& getChild() const {
		return *child;
	}

	/*!
	 * Sets the subtree
	 * @param child the new child node
	 */
	void setChild(AbstractLtlFormula<T>* child) {
		this->child = child;
	}

	/*!
	 *
	 * @return True if the child is set, i.e. it does not point to nullptr; false otherwise
	 */
	bool childIsSet() const {
		return child != nullptr;
	}

	/*!
	 * @returns the maximally allowed number of steps for the bounded until operator
	 */
	uint_fast64_t getBound() const {
		return bound;
	}

	/*!
	 * Sets the maximally allowed number of steps for the bounded until operator
	 *
	 * @param bound the new bound.
	 */
	void setBound(uint_fast64_t bound) {
		this->bound = bound;
	}


private:
	AbstractLtlFormula<T>* child;
	uint_fast64_t bound;

};

} //namespace ltl
} //namespace property
} //namespace storm

#endif /* STORM_FORMULA_LTL_BOUNDEDUNTIL_H_ */
