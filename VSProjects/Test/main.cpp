#include <iostream>
#include "ComparativeTest/ComparativeTest.h"
#include "ImportedCode/StressSolverTest.h"

int main()
{
	std::cout << (ComparativeTest() ? "PASSED!\n" : "NOT PASSED\n");

	using namespace SpecialSolversTest;
	StressStrainStuff::Test3x3x10(-1);

	return 0;
}
