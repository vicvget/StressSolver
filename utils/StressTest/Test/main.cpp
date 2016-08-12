#include <iostream>
#include "ComparativeTest.h"
#include "StressSolverTest.h"
#include "TestFactory.h"

int main()
{
	// RHS calculation based on random data test
	//std::cout << (ComparativeTest() ? "PASSED!\n" : "NOT PASSED\n");


	// Solver  test
	using namespace SpecialSolversTest;
	//StressStrainStuff::Test2x2x5(1);
	//StressStrainStuff::Test2x1x10(1);
	//StressStrainStuff::Test10x1x3(0);
	//StressStrainStuff::Test1x1x3(4);
	//StressStrainStuff::Test1x1x3(3);
	//StressStrainStuff::Test3x3x10(0);
	//StressStrainStuff::Test1x1x3(4);
	//SpecialSolversTest::StressStrainStuff::Test1x1x3(0, SpecialSolversTest::StressStrainStuff::ECode::xlr);
	SpecialSolversTest::StressStrainStuff::Test1x11x11(0, SpecialSolversTest::StressStrainStuff::ECode::xlr);
	//SpecialSolversTest::StressStrainStuff::Test1x11x11(1, SpecialSolversTest::StressStrainStuff::ECode::xlr);
	//SpecialSolversTest::StressStrainStuff::Test10x3x3(0, SpecialSolversTest::StressStrainStuff::ECode::xlr);

	return 0;
}
