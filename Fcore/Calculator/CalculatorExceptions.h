#ifndef CALCULATOR_EXCEPTION_H

#define CALCULATOR_EXCEPTION_H


#include <string>


using std::string;


class CalculatorException
{
public:
	
	const string& ToString() const;

protected:

	string _error;

};


class CEParamNotFound
	:
		public CalculatorException
{
public:

	CEParamNotFound
		(
			const string& param,
			int pos = 0
		);

};


class CEDivisionByZero
	:
		public CalculatorException
{
public:

	CEDivisionByZero
		(
			int pos = 0
		);

};


class CENoBracket
	:
		public CalculatorException
{
public:

	CENoBracket
		(
			int pos = 0
		);

};


class CENotANumber
	:
		public CalculatorException
{
public:

	CENotANumber
		(
			const string& param,
			int pos = 0
		);

};


class CEInvalidToken
	:
		public CalculatorException
{
public:

	CEInvalidToken
		(
			int param,
			int pos = 0
		);

};


class CENotAName
	:
		public CalculatorException
{
public:

	CENotAName
		(
			const string& param,
			int pos = 0
		);

};


class CEInvalidFunction
	:
		public CalculatorException
{
public:

	CEInvalidFunction
		(
			const string& param,
			int pos = 0
		);

};


class CEInvalidValue
	:
		public CalculatorException
{
public:

	CEInvalidValue
		(
			const string& param,
			double value,
			int pos = 0
		);

};


#endif // CALCULATOR_EXCEPTION_H