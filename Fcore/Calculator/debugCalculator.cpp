#include "../BuildType.h"
#ifdef DEBUG_CALCULATOR

#include "../Calculator/Calculator.h"

#include <iostream>


int main()
{
	Calc::Calculator calculator;
	std::cout << calculator.Eval(" 1.0") << '\n' <<
		calculator.Eval(".0") << '\n' <<
		calculator.Eval("0.") << '\n' <<
		calculator.Eval("123456.0") << '\n';
	return 0;
}

#endif