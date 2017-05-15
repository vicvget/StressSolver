#include <iostream>
#include "ComparativeTest.h"
#include "StressSolverTest.h"
#include "TestFactory.h"
#include <sstream>

int main(int argc, char *argv[])
{
	// RHS calculation based on random data test
	//std::cout << (ComparativeTest() ? "PASSED!\n" : "NOT PASSED\n");

	// Solver  test
	using namespace SpecialSolversTest;
	using namespace SpecialSolversTest::StressStrainStuff;
	//Test1x1x3(0, xlr);
	//Test1x51x51(1, xlr);
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
	//Test1xXxXa(1, 0.8, 0.01, 50);

	std::stringstream stringstr;
	int numberOfElements = 50;
	int numberOfSubiterations = 10;
	int solverType = 1;
	string inputIcoFileName = "in.ico";
	string outputIcoFileName = "out.ico";

	if (argc == 2)
	{	// set number of elements
		stringstr << argv[1];
		stringstr >> numberOfElements;
	}
	else
	if (argc == 3)
	{	// set number of subiterations (number of iterations is 100*numberOfSubiterations)
		stringstr << argv[1] << ' ' << argv[2];
		stringstr >> numberOfElements >> numberOfSubiterations;
	}
	else
	if (argc == 4)
	{	// set solverType
		stringstr << argv[1] << ' ' << argv[2] << ' ' << argv[3]; 
		stringstr >> numberOfElements >> numberOfSubiterations >> solverType;
	}

	if (argc == 5)
	{	// set input initial conditions file
		stringstr << argv[1] << ' ' << argv[2] << ' ' << argv[3] << ' ' << argv[4];
		stringstr >> numberOfElements >> numberOfSubiterations >> solverType >> inputIcoFileName;
	}

	if (argc == 6)
	{	// set output initial conditions file
		stringstr << argv[1] << ' ' << argv[2] << ' ' << argv[3] << ' ' << argv[4] << ' ' << argv[5];
		stringstr >> numberOfElements >> numberOfSubiterations >> solverType >> inputIcoFileName >> outputIcoFileName;
	}

	std::cout << "------------------------------" << std::endl
		<< "    STRESS TEST" << std::endl
		<< "------------------------------" << std::endl;
#ifdef INTEL_AVX
	std::cout << "INTEL AVX\n";
#endif
	std::cout << "Elements: " << numberOfElements << 'x' << numberOfElements
		<< "\nIterations: " << numberOfSubiterations << 'x' << 100
		<< "\nSolver Type: " << solverType 
		<< "\nInput ico: " << inputIcoFileName
		<< "\nOutput ico: " << outputIcoFileName
		<< std::endl;

	Test1xXxXa(solverType, 0.8, 0.01, numberOfElements, numberOfSubiterations, inputIcoFileName, outputIcoFileName);
	return 0;
}
