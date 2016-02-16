#include "ComparerTest.h"

#include <iostream>


using Real = double;


struct Iterate
{

	Iterate(...)
	{
	}

};


template
	<
		typename... Reals,
		typename Real
	>
void ComparerTestForTypes
	(
		const Real& number1,
		const Real& number2
	)
{
	Iterate
		{
			(
				ComparerTest<Reals>
					(
						static_cast<Reals>(number1),
						static_cast<Reals>(number2)
					),
				0
			)...
		};
}

int main()
{
	const Real a = 1.2345678901234567890123456789e-20;
	const Real realValue = 3e10;
	const Real result = a / (a / realValue);

	ComparerTestForTypes<float, double, long double>(result, realValue);
}