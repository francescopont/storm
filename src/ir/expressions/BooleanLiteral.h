/*
 * BooleanLiteral.h
 *
 *  Created on: 03.01.2013
 *      Author: chris
 */

#ifndef STORM_IR_EXPRESSIONS_BOOLEANLITERAL_H_
#define STORM_IR_EXPRESSIONS_BOOLEANLITERAL_H_

#include "src/ir/expressions/BaseExpression.h"

namespace storm {
namespace ir {
namespace expressions {

class BooleanLiteral : public BaseExpression {
public:
	bool value;

	BooleanLiteral(bool value) : BaseExpression(bool_), value(value) {

	}

	virtual ~BooleanLiteral() {

	}

	virtual std::shared_ptr<BaseExpression> clone(const std::map<std::string, std::string>& renaming, const std::map<std::string, uint_fast64_t>& bools, const std::map<std::string, uint_fast64_t>& ints) override {
		return std::shared_ptr<BaseExpression>(new BooleanLiteral(this->value));
	}

	virtual bool getValueAsBool(std::pair<std::vector<bool>, std::vector<int_fast64_t>> const* variableValues) const override {
		return value;
	}

	virtual void accept(ExpressionVisitor* visitor) override {
		visitor->visit(this);
	}

	virtual std::string toString() const override {
		if (value) {
			return std::string("true");
		} else {
			return std::string("false");
		}
	}
	
	virtual std::string dump(std::string prefix) const override {
		std::stringstream result;
		result << prefix << "BooleanLiteral " << this->toString() << std::endl;
		return result.str();
	}
};

}

}

}

#endif /* STORM_IR_EXPRESSIONS_BOOLEANLITERAL_H_ */
