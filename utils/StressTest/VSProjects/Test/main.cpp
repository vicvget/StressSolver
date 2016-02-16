#include "StressStrainSolverExports.h"
#include "StressStrainFortranIterativeSolver.h"
#include "StressStrainCppIterativeSolver.h"
#include <cmath>

#define M_PI 3.1415926535897932384626433832795
using std::setw;

// X1,Y2,Z3 rotation (airplane angles)
void MakeRotationMatrix(double UE[3], double* a)
{
	const int id = 0;
	double xc = cos(UE[id]);
	double yc = cos(UE[id + 1]);
	double zc = cos(UE[id + 2]);
	double xs = sin(UE[id]);
	double ys = sin(UE[id + 1]);
	double zs = sin(UE[id + 2]);


	double* firstRow = a;
	double* secondRow = a + vecStride;
	double* thirdRow = a + vecStride * 2;

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
	int nNodes = 100;

	for (int i = 0; i < nNodes*3; i++)
	{
		nodes[i] = rand();

	}

	Stress::StressStrainCppIterativeSolver* cppSolver = new Stress::StressStrainCppIterativeSolver
		(
		params,
		links,
		nLinks,
		nodes,
		nNodes,
		1,
		1,
		1
		);

	StressStrainFortranIterativeSolver* hsolver = new StressStrainFortranIterativeSolver
		(
		params,
		links,
		nLinks,
		nodes,
		nNodes,
		1,
		1,
		1
		);

	double UE[3];
	for (int i = 0; i < nNodes; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			UE[j] = ((double)rand())/RAND_MAX * 2*M_PI;
		}
		MakeRotationMatrix(UE, cppSolver->_dataRotationMtx + matStride * i);
	}

	for (int i = 0; i < nNodes; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			*(cppSolver->GetElementVelocity(i)+j) = rand();
			*(cppSolver->GetElementVelocityAngular(i)+ j) = rand();
		}

	}


	__declspec(align(32)) double SL[8] = { 0 }, SL2[8] = { 0 }, VL[8] = { 0 }, VL2[8] = { 0 }, rx, ry, rz;

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
		nNodes		// количество узлов
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
		nNodes		// количество узлов
		);

	const int width = 9;
	for (int i = 0; i < vecStride; i++)
	{
		std::cout << setw(9) << SL[i] << " " << setw(9) << SL2[i] << std::endl;
		std::cout << setw(9) << VL[i] << " " << setw(9) << VL2[i] << std::endl;
	}

	delete hsolver;
	delete cppSolver;

	return 0;
}