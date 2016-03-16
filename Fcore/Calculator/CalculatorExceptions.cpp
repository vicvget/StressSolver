#include "CalculatorExceptions.h"

#include <sstream>


using std::ostringstream;


// class CalculatorException

const string& CalculatorException::ToString() const
{
	return _error;
}


// class CEParamNotFound

CEParamNotFound::CEParamNotFound
	(
		const string& param,
		int pos
	)
{
	ostringstream stringStream;

	stringStream << "Error at " << pos << ": Parameter [" << param << "] not found";
	_error = stringStream.str();
}


// class CEDivisionByZero

CEDivisionByZero::CEDivisionByZero
	(
		int pos
	)
{
	ostringstream stringStream;

	stringStream << "Error at " << pos << ": Division by 0";
	_error = stringStream.str();
}


// class CENoBracket

CENoBracket::CENoBracket
	(
		int pos
	)
{
	ostringstream stringStream;

	stringStream << "Error at " << pos << ": No bracket";
	_error = stringStream.str();
}


// class CENotANumber

CENotANumber::CENotANumber
	(
		const string& param,
		int pos
	)
{
	ostringstream stringStream;

	stringStream << "Error at " << pos << ": Not a number: " << param;
	_error = stringStream.str();
}


// class CEInvalidToken

CEInvalidToken::CEInvalidToken
	(
		int param,
		int pos
	)
{
	ostringstream stringStream;

	stringStream << "Error at " << pos << ": Invalid token: " << param;
	_error = stringStream.str();
}


// class CENotAName

CENotAName::CENotAName
	(
		const string& param,
		int pos
	)
{
	ostringstream stringStream;

	stringStream << "Error at " << pos << ": Not a valid name: " << param;
	_error = stringStream.str();
}


// class CEInvalidFunction

CEInvalidFunction::CEInvalidFunction
	(
		const string& param,
		int pos
	)
{
	ostringstream stringStream;

	stringStream << "Error at " << pos << ": Invalid function: " << param;
	_error = stringStream.str();
}


// class CEInvalidValue

CEInvalidValue::CEInvalidValue
	(
		const string& param,
		double value,
		int pos
	)
{
	ostringstream stringStream;

	stringStream << "Error at " << pos << " " << param << " [" << value << "]";
	_error = stringStream.str();
}