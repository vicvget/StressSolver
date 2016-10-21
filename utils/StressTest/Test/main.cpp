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

	//SpecialSolversTest::StressStrainStuff::Test1x1x3(1, SpecialSolversTest::StressStrainStuff::ECode::xlr);
	//SpecialSolversTest::StressStrainStuff::Test1x51x51(1, SpecialSolversTest::StressStrainStuff::ECode::xlr);
	//SpecialSolversTest::StressStrainStuff::Test1x11x11(1, SpecialSolversTest::StressStrainStuff::ECode::xlr);
	//SpecialSolversTest::StressStrainStuff::Test1x22x22a(0, SpecialSolversTest::StressStrainStuff::ECode::xlr);
	SpecialSolversTest::StressStrainStuff::Test1x44x44a(0, SpecialSolversTest::StressStrainStuff::ECode::xlr);
	//SpecialSolversTest::StressStrainStuff::Test1x100x100a(0, SpecialSolversTest::StressStrainStuff::ECode::xlr);
	//SpecialSolversTest::StressStrainStuff::Test10x3x3(1, SpecialSolversTest::StressStrainStuff::ECode::xlr);
	//SpecialSolversTest::StressStrainStuff::Test50x5x5(1, SpecialSolversTest::StressStrainStuff::ECode::xlr);
	//SpecialSolversTest::StressStrainStuff::Test10x5x5(1, SpecialSolversTest::StressStrainStuff::ECode::xlr);
	//SpecialSolversTest::StressStrainStuff::Test10x7x7(1, SpecialSolversTest::StressStrainStuff::ECode::xlr);

	return 0;
}
