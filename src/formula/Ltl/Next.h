/*
 * Next.h
 *
 *  Created on: 19.10.2012
 *      Author: Thomas Heinemann
 */

#ifndef STORM_FORMULA_LTL_NEXT_H_
#define STORM_FORMULA_LTL_NEXT_H_

#include "AbstractLtlFormula.h"
#include "src/formula/abstract/Next.h"
#include "src/formula/AbstractFormulaChecker.h"

namespace storm {
namespace property {
namespace ltl {

template <class T> class Next;

/*!
 *  @brief Interface class for model checkers that support Next.
 *   
 *  All model checkers that support the formula class Next must inherit
 *  this pure virtual class.
 */
template <class T>
class INextModelChecker {
    public:
		/*!
         *  @brief Evaluates Next formula within a model checker.
         *
         *  @param obj Formula object with subformulas.
         *  @return Result of the formula for every node.
         */
        virtual std::vector<T>* checkNext(const Next<T>& obj) const = 0;
};

/*!
 *	@brief Interface class for visitors that support Next.
 *
 *	All visitors that support the formula class Next must inherit
 *	this pure virtual class.
 */
template <class T>
class INextVisitor {
	public:
		/*!
		 *	@brief Visits Next formula.
		 *
		 *	@param obj Formula object with subformulas.
		 *	@return Result of the formula for every node.
		 */
		virtual void visitNext(const Next<T>& obj) = 0;
};

/*!
 * @brief
 * Class for an abstract (path) formula tree with a Next node as root.
 *
 * Has two Abstract LTL formulas as sub formulas/trees.
 *
 * @par Semantics
 * The formula holds iff in the next step, \e child holds
 *
 * The subtree is seen as part of the object and deleted with the object
 * (this behavior can be prevented by setting them to NULL before deletion)
 *
 * @see AbstractLtlFormula
 */
template <class T>
class Next : public storm::property::abstract::Next<T, AbstractLtlFormula<T>>,
				 public AbstractLtlFormula<T> {

public:
	/*!
	 * Empty constructor
	 */
	Next() {
		//intentionally left empty
	}

	/*!
	 * Constructor
	 *
	 * @param child The child node
	 */
	Next(AbstractLtlFormula<T>* child)
		: storm::property::abstract::Next<T, AbstractLtlFormula<T>>(child) {
		//intentionally left empty
	}

	/*!
	 * Constructor.
	 *
	 * Also deletes the subtree.
	 * (this behaviour can be prevented by setting the subtrees to NULL before deletion)
	 */
	virtual ~Next() {
	  //intentionally left empty
	}

	/*!
	 * Clones the called object.
	 *
	 * Performs a "deep copy", i.e. the subtrees of the new object are clones of the original ones
	 *
	 * @returns a new BoundedUntil-object that is identical the called object.
	 */
	virtual AbstractLtlFormula<T>* clone() const override {
		Next<T>* result = new Next<T>();
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
	virtual std::vector<T> *check(const storm::modelchecker::ltl::AbstractModelChecker<T>& modelChecker) const {
		return modelChecker.template as<INextModelChecker>()->checkNext(*this);
	}

	virtual void visit(visitor::AbstractLtlFormulaVisitor<T>& visitor) const {
		visitor.template as<INextVisitor>()->visitNext(*this);
	}
};

} //namespace ltl
} //namespace property
} //namespace storm

#endif /* STORM_FORMULA_LTL_NEXT_H_ */
