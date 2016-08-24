#include "../Solvers/Stress/StressStrainCppIterativeSolver.h"
#ifndef USE_KNC
#include "../Solvers/Stress/StressStrainCppIterativeSolverAVX.h"
#include "../Solvers/Stress/StressStrainCppIterativeSolverFMA.h"
#endif // !USE_KNC
#include "../AdditionalModules/fmath/Matrix3x4.h"
#include "../AdditionalModules/fmath/Matrix3x3.h"

#include <memory>

#define M_PI 3.1415926535897932384626433832795
using std::setw;
using MathHelpers::Mat3x4;
using MathHelpers::Mat3;

// X1,Y2,Z3 rotation (airplane angles)
void MakeRotationMatrix(double* ue, double* a, int stride, bool transposed = false)
{
	if (stride == 4)
	{
		Mat3x4::MakeXYZRotationMtx01(ue).Export(a);
	}
	else if (transposed)
	{
		Mat3::MakeXYZRotationMtx10(ue).Export(a);
	}
	else
	{
		Mat3::MakeXYZRotationMtx01(ue).Export(a);
	}
}

int CompareArrays(double* array1, double* array2, int count)
{
	for (size_t i = 0; i < count; i++)
	{
		if (std::abs(array1[i] - array2[i]) > 1e-8)
			return (int)i;
	}
	return 0;
}

bool ComparativeTest()
{
	double params[4] = { 1e10, 0.01, 1, 1e8 };
	double nodes[300];
	double gridStep = 1;
	double timeStep = 1e-4;

	int links[1];
	int nLinks = 0;
	int nElements = 100;
	int nThreads = 1;
	int for_stride = 3;
	int cpp_stride = 4;

	// random coordinates
	for (int i = 0; i < nElements * 3; i++)
	{
		nodes[i] = rand();

	}

	// create tested solver
	std::shared_ptr<Stress::StressStrainCppIterativeSolver> cppSolver =
		std::make_shared < Stress::StressStrainCppIterativeSolver >
		(
		params,
		links,
		nLinks,
		nodes,
		nElements,
		gridStep,
		timeStep,
		nThreads,
		cpp_stride
		);

	// create tested solver
	std::shared_ptr<Stress::StressStrainCppIterativeSolver> cppSolverUa =
		std::make_shared < Stress::StressStrainCppIterativeSolver >
		(
		params,
		links,
		nLinks,
		nodes,
		nElements,
		gridStep,
		timeStep,
		nThreads,
		3
		);

	// make random matrices
	double UE[3];
	for (size_t i = 0; i < nElements; i++)
	{
		for (size_t j = 0; j < 3; j++)
		{
			UE[j] = ((double)rand()) / RAND_MAX * 2 * M_PI;
		}
		MakeRotationMatrix(UE, cppSolver->GetRotationMatrix(i), cpp_stride);
		MakeRotationMatrix(UE, cppSolverUa->GetRotationMatrix(i), 3);
	}

	// make random velocities and speeds
	for (int i = 0; i < nElements; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			double u = (double)rand() / RAND_MAX * 100;
			double v = (double)rand() / RAND_MAX * 100;
			double w = (double)rand() / RAND_MAX * 100;
			*(cppSolver->GetElementShiftAngular(i) + j) = u;
			*(cppSolver->GetElementVelocity(i) + j) = v;
			*(cppSolver->GetElementVelocityAngular(i) + j) = w;
			*(cppSolverUa->GetElementShiftAngular(i) + j) = u;
			*(cppSolverUa->GetElementVelocity(i) + j) = v;
			*(cppSolverUa->GetElementVelocityAngular(i) + j) = w;
		}

	}

	// output arrays
	__declspec(align(64)) double SL[8] = { 0 }, VL[8] = { 0 };
	__declspec(align(64)) double SL2[8] = { 0 }, VL2[8] = { 0 };
	__declspec(align(64)) double SL3[8] = { 0 }, VL3[8] = { 0 };
	__declspec(align(64)) double SL4[8] = { 0 }, VL4[8] = { 0 };
	__declspec(align(64)) double SL5[8] = { 0 }, VL5[8] = { 0 };
	__declspec(align(64)) double SL6[8] = { 0 }, VL6[8] = { 0 };

	double cVec1[] = { -gridStep * 0.5, 0, 0 };
	double cVec2[] = { gridStep * 0.5, 0, 0 };
//	double A[36], C[36];

//	forSolver->linksh(cVec1, cVec2, SL3, VL3, A, C, 0, 1, nElements);
//	forSolver->linksh2(cVec1, cVec2, SL4, VL4, A, C, 0, 1, nElements);
//	forSolver->linksh3(cVec1, cVec2, SL5, VL5, A, C, 0, 1, nElements);

	cppSolverUa->CalculateStrainsUa
		(
		0,	// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
		&SL6[0],		// выход деформаций
		&VL6[0],		// выход изм. скоростей
		1,	// номер узла 1
		2	// номер узла 2
		);


	//std::cout << "CalculateStrains" << std::endl;
	//cppSolver->CalculateStrains
	//	(
	//	0,	// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
	//	&SL[0],		// выход деформаций
	//	&VL[0],		// выход изм. скоростей
	//	1,	// номер узла 1
	//	2	// номер узла 2
	//	);
	//cppSolver->CalculateStrains
	//	(
	//	0,	// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
	//	&SL[0],		// выход деформаций
	//	&VL[0],		// выход изм. скоростей
	//	1,	// номер узла 1
	//	3	// номер узла 2
	//	);
	//cppSolver->CalculateStrains
	//	(
	//	0,	// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
	//	&SL[0],		// выход деформаций
	//	&VL[0],		// выход изм. скоростей
	//	0,	// номер узла 1
	//	1	// номер узла 2
	//	);

#ifndef USE_KNC
	cppSolver->CalculateStrainsAVX
		(
		0,	// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
		&SL2[0],		// выход деформаций 
		&VL2[0],		// выход деформаций
		1,	// номер узла 1
		2	// номер узла 2
		);

	//cppSolver->CalculateStrainsSSE
	//	(
	//	0,	// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
	//	&SL3[0],		// выход деформаций 
	//	&VL3[0],		// выход деформаций
	//	0,	// номер узла 1
	//	1	// номер узла 2
	//	);

	//cppSolver->CalculateStrainsFMA
	//	(
	//	0,	// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
	//	&SL4[0],		// выход деформаций 
	//	&VL4[0],		// выход деформаций
	//	0,	// номер узла 1
	//	1	// номер узла 2
	//	);
#else
	std::cout << "CalculateStrainsKNC" << std::endl;
	cppSolver->CalculateStrainsKNC
		(
		0,	// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
		&SL5[0],		// выход деформаций 
		&VL5[0],		// выход деформаций
		1,	// номер узла 1
		2	// номер узла 2
		);
	cppSolver->CalculateStrainsKNC
		(
		0,	// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
		&SL5[0],		// выход деформаций 
		&VL5[0],		// выход деформаций
		1,	// номер узла 1
		3	// номер узла 2
		);
	cppSolver->CalculateStrainsKNC
	(
		0,	// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
		&SL5[0],		// выход деформаций 
		&VL5[0],		// выход деформаций
		0,	// номер узла 1
		1	// номер узла 2
	);

#endif

	const int width = 9;
	std::cout  << "SL " << 
		"cpp " << setw(width) << 
		"avx " << setw(width) << 
		"sse"  << setw(width) << 
		"fma"  << setw(width) << 
		"knc"  << setw(width) << 
		"ua"   << std::endl;
	
	for (int i = 0; i < 8; i++)
	{
		std::cout << setw(width) << SL[i] << " "
			<< setw(width) << SL2[i] << " "
			<< setw(width) << SL3[i] << " "
			<< setw(width) << SL4[i] << " "
			<< setw(width) << SL5[i] << " "
			<< setw(width) << SL6[i] << " "
			<< std::endl;
	}
	std::cout << "VL\n";
	for (int i = 0; i < 8; i++)
	{
		std::cout << setw(width) << VL[i] << " "
			<< setw(width) << VL2[i] << " "
			<< setw(width) << VL3[i] << " "
			<< setw(width) << VL4[i] << " "
			<< setw(width) << VL5[i] << " "
			<< setw(width) << VL6[i] << " "
			<< std::endl;
	}

	bool res = true;
	int cmp = 0;

	cmp = CompareArrays(SL, SL2, 8);
	if (cmp)
	{
		std::cout << "DIFF: SL,SL2 in element " << cmp << std::endl;
		res = false;
	}
	cmp = CompareArrays(SL, SL3, 8);
	if (cmp)
	{
		std::cout << "DIFF: SL,SL3 in element " << cmp << std::endl;
		res = false;
	}
	
	cmp = CompareArrays(SL, SL4, 8);
	if (cmp)
	{
		std::cout << "DIFF: SL,SL4 in element " << cmp << std::endl;
		res = false;
	}

	cmp = CompareArrays(SL, SL5, 8);
	if (cmp)
	{
		std::cout << "DIFF: SL,SL5 in element " << cmp << std::endl;
		res = false;
	}

	cppSolver->Solve(10);
	return res;
}