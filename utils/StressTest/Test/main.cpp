#include <iostream>
#include "ComparativeTest.h"
#include "StressSolverTest.h"

int main()
{
	// RHS calculation based on random data test
	//std::cout << (ComparativeTest() ? "PASSED!\n" : "NOT PASSED\n");


	// Solver  test
	using namespace SpecialSolversTest;
	StressStrainStuff::Test2x1x10(1);
	//StressStrainStuff::Test1x3x10(0);
	//StressStrainStuff::Test1x1x3(4);
	//StressStrainStuff::Test1x1x3(3);
	//StressStrainStuff::Test3x3x10(0);
	//StressStrainStuff::Test1x1x3(4);

	return 0;
}
