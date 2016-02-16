#include "ComparerTest.h"

#include "../RealComparing.h"

#include <iostream>
#include <typeinfo>


template
	<
		typename Real
	>
void ComparerTest
	(
		const Real& number1,
		const Real& number2
	)
{
	std::cout
		<< typeid(Real).name()
		<< " comparer test"
		<< std::endl;
	std::cout
		<< "Equal("
		<< RealComparing<Real>::Equal(number1, number2)
		<< "): "
		<< number1
		<< " == "
		<< number2
		<< std::endl;
	std::cout
		<< "Unequal("
		<< RealComparing<Real>::Unequal(number1, number2)
		<< "): "
		<< number1
		<< " != "
		<< number2
		<< std::endl;
	std::cout
		<< "Greater("
		<< RealComparing<Real>::Greater(number1, number2)
		<< "): "
		<< number1
		<< " > "
		<< number2
		<< std::endl;
	std::cout
		<< "Less("
		<< RealComparing<Real>::Less(number1, number2)
		<< "): "
		<< number1
		<< " < "
		<< number2
		<< std::endl;
	std::cout
		<< "GreaterOrEqual("
		<< RealComparing<Real>::GreaterOrEqual(number1, number2)
		<< "): "
		<< number1
		<< " >= "
		<< number2
		<< std::endl;
	std::cout
		<< "LessOrEqual("
		<< RealComparing<Real>::LessOrEqual(number1, number2)
		<< "): "
		<< number1
		<< " <= "
		<< number2
		<< std::endl;
	std::cout << std::endl;
}


#define COMPARER_TEST_INSTANTIATION(Real) \
 \
template \
void ComparerTest<Real> \
	( \
		const Real& number1, \
		const Real& number2 \
	);

COMPARER_TEST_INSTANTIATION(float)

COMPARER_TEST_INSTANTIATION(double)

COMPARER_TEST_INSTANTIATION(long double)