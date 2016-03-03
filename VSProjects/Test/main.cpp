#include "StressStrainSolverExports.h"
#include "StressStrainFortranIterativeSolver.h"
#include "StressStrainCppIterativeSolver.h"
#include <cmath>
#include "../../AdditionalModules/fmath/Vector3.h"
#include "../../AdditionalModules/fmath/Matrix3x4.h"
#include "../../AdditionalModules/fmath/Matrix3x3.h"

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

int main()
{
	double params[4] = {1e10, 0.01, 1, 1e8};
	int links[1];
	double nLinks = 0;
	double nodes[300];
	int nElements = 100;
	double gridStep = 1;
	double timeStep = 1e-4;
	double nThreads = 1;
	double for_stride = 3;
	double cpp_stride = 4;

	for (int i = 0; i < nElements*3; i++)
	{
		nodes[i] = rand();

	}

	Stress::StressStrainCppIterativeSolver* cppSolver = new Stress::StressStrainCppIterativeSolver
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

	StressStrainFortranIterativeSolver* forSolver = new StressStrainFortranIterativeSolver
		(
		params,
		links,
		nLinks,
		nodes,
		nElements,
		gridStep,
		timeStep,
		nThreads,
		for_stride
		);

	double UE[3];
	for (int i = 0; i < nElements; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			UE[j] = ((double)rand())/RAND_MAX * 2*M_PI;
		}
		MakeRotationMatrix(UE, cppSolver->GetRotationMatrix(i), cpp_stride);
		MakeRotationMatrix(UE, forSolver->GetRotationMatrix(i), for_stride, true);
	}

	for (int i = 0; i < nElements; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			double v = (double)rand() / RAND_MAX * 100;
			double w = (double)rand() / RAND_MAX * 100;
			*(cppSolver->GetElementVelocity(i) + j) = v;
			*(cppSolver->GetElementVelocityAngular(i)+ j) = w;
			*(forSolver->GetElementVelocity(i) + j) = v;
			*(forSolver->GetElementVelocityAngular(i) + j) = w;
		}

	}

	__declspec(align(32)) double SL[8] = { 0 }, VL[8] = { 0 };
	__declspec(align(32)) double SL2[8] = { 0 }, VL2[8] = { 0 };
	__declspec(align(32)) double SL3[8] = { 0 }, VL3[8] = { 0 };
	__declspec(align(32)) double SL4[8] = { 0 }, VL4[8] = { 0 };
	__declspec(align(32)) double SL5[8] = { 0 }, VL5[8] = { 0 };
	double rx, ry, rz;

	double cVec1[] = { -gridStep * 0.5, 0, 0 };
	double cVec2[] = { gridStep * 0.5, 0, 0 };
	double A[36], C[36];

	forSolver->linksh(cVec1, cVec2, SL3, VL3, A, C, 0, 1, nElements);
	forSolver->linksh2(cVec1, cVec2, SL4, VL4, A, C, 0, 1, nElements);
	forSolver->linksh3(cVec1, cVec2, SL5, VL5, A, C, 0, 1, nElements);

	cppSolver->CalculateStrains
		(
		0,	// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
		&SL[0],		// выход деформаций
		&VL[0],		// выход изм. скоростей
		0,	// номер узла 1
		1	// номер узла 2
		);


	cppSolver->CalculateStrainsAVX
		(
		0,	// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
		&SL2[0],		// выход деформаций
		&VL2[0],		// выход деформаций
		0,	// номер узла 1
		1	// номер узла 2
		);

	const int width = 9;
	std::cout << "SL cpp " << setw(width) << "avx " << setw(width) << "l1" << setw(width) << "l2" << setw(width) << "l3\n";
	for (int i = 0; i < 8; i++)
	{
		std::cout << setw(width) << SL[i] << " "
			<< setw(width) << SL2[i] << " "
			<< setw(width) << SL3[i] << " "
			<< setw(width) << SL4[i] << " "
			<< setw(width) << SL5[i] << " "
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
					<< std::endl;
	}

	cppSolver->Solve(10);

	delete forSolver;
	delete cppSolver;

	return 0;
}