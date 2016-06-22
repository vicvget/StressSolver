#include <iostream>
#include "ComparativeTest.h"
#include "StressSolverTest.h"

int main()
{
	// RHS calculation based on random data test
	std::cout << (ComparativeTest() ? "PASSED!\n" : "NOT PASSED\n");


	// Solver  test
	//using namespace SpecialSolversTest;
	//StressStrainStuff::Test3x1x1(3);

	return 0;
}
