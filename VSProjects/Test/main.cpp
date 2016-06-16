#include <iostream>
#include "ComparativeTest/ComparativeTest.h"
  #include "ImportedCode/StressSolverTest.h"

int main()
{
	// RHS calculation based on random data test
	std::cout << (ComparativeTest() ? "PASSED!\n" : "NOT PASSED\n");


	// Solver  test
	using namespace SpecialSolversTest;
	StressStrainStuff::Test1x3x10(0);

	return 0;
}
