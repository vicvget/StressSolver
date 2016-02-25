#include "StressStrainSolverExports.h"
#include "StressStrainFortranIterativeSolver.h"
#include "StressStrainCppIterativeSolver.h"
#include <cmath>
#include "../../AdditionalModules/fmath/Vector3.h"

#define M_PI 3.1415926535897932384626433832795
using std::setw;

// X1,Y2,Z3 rotation (airplane angles)
void MakeRotationMatrix(double UE[3], double* a, int stride)
{
	const int id = 0;
	double xc = cos(UE[id]);
	double yc = cos(UE[id + 1]);
	double zc = cos(UE[id + 2]);
	double xs = sin(UE[id]);
	double ys = sin(UE[id + 1]);
	double zs = sin(UE[id + 2]);


	double* firstRow = a;
	double* secondRow = a + stride;
	double* thirdRow = a + stride * 2;

	firstRow[0] = yc*zc;
	firstRow[1] = -yc*zs;
	firstRow[2] = ys;

	secondRow[0] = xs*ys*zc + xc*zs;
	secondRow[1] = -xs*ys*zs + xc*zc;
	secondRow[2] = -xs*yc;

	thirdRow[0] = -xc*ys*zc + xs*zs;
	thirdRow[1] = xc*ys*zs + xs*zc;
	thirdRow[2] = xc*yc;
}

int main()
{
	double params[4];
	int links[1];
	double nLinks = 0;
	double nodes[300];
	int nElements = 100;
	double gridStep = 1;
	double timeStep = 1;
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
		MakeRotationMatrix(UE, forSolver->GetRotationMatrix(i), for_stride);
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

	cppSolver->linksh4
		(
		0,	// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
		&SL[0],		// выход деформаций
		&VL[0],		// выход изм. скоростей
		rx,		// выход
		ry,		// выход 
		rz,		// выход 
		0,	// номер узла 1
		1,	// номер узла 2
		nElements		// количество
		);


	cppSolver->linksh4AVX
		(
		0,	// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
		&SL2[0],		// выход деформаций
		&VL2[0],		// выход изм. скоростей
		rx,		// выход
		ry,		// выход 
		rz,		// выход 
		0,	// номер узла 1
		1,	// номер узла 2
		nElements		// количество узлов
		);

	const int width = 9;
	std::cout << "SL l4 " << setw(width) << "lavx " << setw(width) << "l1" << setw(width) << "l2" << setw(width) << "l3\n";
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

	//cppSolver->Solve(10000);

	delete forSolver;
	delete cppSolver;

	return 0;
}