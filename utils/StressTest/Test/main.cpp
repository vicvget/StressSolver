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
	using namespace SpecialSolversTest::StressStrainStuff;
	//Test1x1x3(0, xlr);
	Test1x51x51(1, xlr);
	//Test1x11x11(1, xlr);
	//Test1x22x22a(0, xlr);
	//Test1x44x44a(0, xlr);
	//Test1x100x100a(0, xlr);
	//Test10x3x3(1, xlr);
	//Test50x5x5(1, xlr);
	//Test10x5x5(1, xlr);
	//Test10x7x7(1, xlr);
	//Test1xXxXa(0, 0.8, 0.01, 10);
	//Test1xXxXa(0, 0.8, 0.01, 20);
	//Test1xXxXa(0, 0.8, 0.01, 30);
	//Test1xXxXa(0, 0.8, 0.01, 100);
	//Test1xXxXa(0, 0.8, 0.01, 80);
	//Test1xXxXa(0, 0.8, 0.01, 1000);
	//Test1xXxXa(0, 0.8, 0.01, 50);
	//Test1xXxXa(0, 0.8, 0.01, 20);
	//Test1xXxXa(0, 0.8, 0.01, 100);
	return 0;
}
