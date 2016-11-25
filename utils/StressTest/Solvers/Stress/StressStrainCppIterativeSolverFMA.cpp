#ifndef USE_KNC

#include "StressStrainCppIterativeSolverFMA.h"
#include "../../AdditionalModules/fmath/Matrix3x3.h"
#include "../../AdditionalModules/fmath/Matrix3x4.h"

using MathHelpers::Vec3;
using MathHelpers::Vec3Ref;
using MathHelpers::MakeVec3;
using MathHelpers::Mat3;
using MathHelpers::Mat3x4;

namespace Stress
{

StressStrainCppIterativeSolverFMA::StressStrainCppIterativeSolverFMA
	(
		double* params, 
		int* links, 
		int nLinks, 
		double *gridElements, 
		int nElements, 
		double gridStep, 
		double timeStep,
		int numThreads,
		int stride
	)
	:
		StressStrainCppIterativeSolver
			(
				params,
				links, 
				nLinks, 
				gridElements, 
				nElements, 
				gridStep,
				timeStep,
				numThreads,
				stride
			)
{
	std::cout << "FMA SOLVER" << std::endl << std::flush;
}

// virtual
StressStrainCppIterativeSolverFMA::~StressStrainCppIterativeSolverFMA()
{

}

// AVX-������
void StressStrainCppIterativeSolverFMA::SolveFull(const int nIterations)
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	__m256d timeStep;
	__m256d timeStep2;
	__m256d timeStep4;
	__m256d constantD2;
	__m256d constantD6;

	timeStep = _mm256_set1_pd(_timeStep);
	timeStep2 = _mm256_set1_pd(_timeStep2);
	timeStep4 = _mm256_set1_pd(_timeStep4);
	constantD2 = _mm256_set1_pd(0.5);
	constantD6 = _mm256_set1_pd(1 / 6.0);

	FTimer test_timer;
	_iterationNumber = 0;
	_testTimer.Start(0);

	__m256d Xtmp, DXtmp, DDXtmp, tmp;
	__m256d hDDX1, hDDX2, hDDX3, sDDX;

	while (_iterationNumber != nIterations && _rotationSolver->IsValid())
	{
		_iterationNumber++;
		_nIteration++;

		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
		memcpy(_initX, _varX, sizeof(double)*_nVariables);
		memcpy(_initDX, _varDX, sizeof(double)*_nVariables);

		_testTimer.Start(3);
		#pragma omp parallel for num_threads(_numThreads)
		for (int j = 0; j < _nVariables; j += regSize)
		{
			//_hDDX1[j] = _varDDX[j] * _timeStep;
			//_varX[j] += _varDX[j] * _timeStep2;
			//_varDX[j] += _hDDX1[j] * 0.5;

			Xtmp = _mm256_load_pd(_varX + j);
			DXtmp = _mm256_load_pd(_varDX + j);
			DDXtmp = _mm256_load_pd(_varDDX + j);

			hDDX1 = _mm256_mul_pd(DDXtmp, timeStep);
			_mm256_store_pd(_hDDX1 + j, hDDX1);
			tmp = _mm256_fmadd_pd(DXtmp, timeStep2, Xtmp);
			//tmp = _mm256_add_pd(_mm256_mul_pd(DXtmp, timeStep2), Xtmp);
			_mm256_store_pd(_varX + j, tmp);
			tmp = _mm256_fmadd_pd(hDDX1, constantD2, DXtmp);
			//tmp = _mm256_add_pd(_mm256_mul_pd(hDDX1, constant), DXtmp);
			_mm256_store_pd(_varDX + j, tmp);
		}
		_testTimer.Stop(3);

		MeasuredRun(1, _rotationSolver->Solve2());
		MeasuredRun(2, CalculateForces());

		_testTimer.Start(3);
#pragma omp parallel for num_threads(_numThreads)
		for (int j = 0; j < _nVariables; j += regSize)
		{
			//_hDDX2[j] = _varDDX[j] * _timeStep;
			//_varX[j] += _hDDX1[j] * _timeStep4;
			//_varDX[j] = _initDX[j] + _hDDX2[j] * 0.5;

			Xtmp = _mm256_load_pd(_varX + j);
			DXtmp = _mm256_load_pd(_initDX + j);
			DDXtmp = _mm256_load_pd(_varDDX + j);
			hDDX1 = _mm256_load_pd(_hDDX1 + j);

			hDDX2 = _mm256_mul_pd(DDXtmp, timeStep);
			_mm256_store_pd(_hDDX2 + j, hDDX2);
			tmp = _mm256_fmadd_pd(hDDX1, timeStep4, Xtmp);
			//tmp = _mm256_add_pd(_mm256_mul_pd(hDDX1, timeStep4), Xtmp);
			_mm256_store_pd(_varX + j, tmp);
			tmp = _mm256_fmadd_pd(hDDX2, constantD2, DXtmp);
			//tmp = _mm256_add_pd(_mm256_mul_pd(hDDX2, constant), DXtmp);
			_mm256_store_pd(_varDX + j, tmp);
		}
		_testTimer.Stop(3);

		MeasuredRun(1, _rotationSolver->Solve3());
		MeasuredRun(2, CalculateForces());

		_testTimer.Start(3);
#pragma omp parallel for num_threads(_numThreads)
		for (int j = 0; j < _nVariables; j += regSize)
		{
			//_hDDX3[j] = _varDDX[j] * _timeStep;
			//_varX[j] = _initX[j] + (_initDX[j] + _hDDX2[j] * 0.5) * _timeStep;
			//_varDX[j] = _initDX[j] + _hDDX3[j];

			Xtmp = _mm256_load_pd(_initX + j);
			DXtmp = _mm256_load_pd(_initDX + j);
			DDXtmp = _mm256_load_pd(_varDDX + j);
			hDDX2 = _mm256_load_pd(_hDDX2 + j);

			hDDX3 = _mm256_mul_pd(DDXtmp, timeStep);
			_mm256_store_pd(_hDDX3 + j, hDDX3);
			tmp = _mm256_fmadd_pd(hDDX2, constantD2, DXtmp);
			//tmp = _mm256_add_pd(_mm256_mul_pd(hDDX2, constant), DXtmp);
			tmp = _mm256_fmadd_pd(tmp, timeStep, Xtmp);
			//tmp = _mm256_add_pd(_mm256_mul_pd(tmp, timeStep), Xtmp);
			_mm256_store_pd(_varX + j, tmp);
			_mm256_store_pd(_varDX + j, _mm256_add_pd(DXtmp, hDDX3));
		}
		_testTimer.Stop(3);

		MeasuredRun(1, _rotationSolver->Solve4());
		MeasuredRun(2, CalculateForces());

		_testTimer.Start(3);
#pragma omp parallel for num_threads(_numThreads)
		for (int j = 0; j < _nVariables; j += regSize)
		{
			//float sDDX = _hDDX2[j] + _hDDX3[j];
			//_varX[j] = _initX[j] + (_initDX[j] + sDDX / 6.0) * _timeStep;
			//_varDX[j] = _initDX[j] + (_hDDX1[j] + sDDX + sDDX + _varDDX[j] * _timeStep) / 6.0;

			Xtmp = _mm256_load_pd(_initX + j);
			DXtmp = _mm256_load_pd(_initDX + j);
			DDXtmp = _mm256_load_pd(_varDDX + j);
			hDDX1 = _mm256_load_pd(_hDDX1 + j);
			hDDX2 = _mm256_load_pd(_hDDX2 + j);
			hDDX3 = _mm256_load_pd(_hDDX3 + j);

			sDDX = _mm256_add_pd(hDDX2, hDDX3);
			tmp = _mm256_fmadd_pd(sDDX, constantD6, DXtmp);
			//tmp = _mm256_add_pd(_mm256_mul_pd(sDDX, constant), DXtmp);
			tmp = _mm256_fmadd_pd(tmp, timeStep, Xtmp);
			//tmp = _mm256_add_pd(_mm256_mul_pd(tmp, timeStep), Xtmp);
			_mm256_store_pd(_varX + j, tmp);
			tmp = _mm256_add_pd(_mm256_add_pd(hDDX1, sDDX), sDDX);
			tmp = _mm256_fmadd_pd(DDXtmp, timeStep, tmp);
			//tmp = _mm256_add_pd(_mm256_mul_pd(DDXtmp, timeStep), tmp);
			tmp = _mm256_fmadd_pd(tmp, constantD6, DXtmp);
			//tmp = _mm256_add_pd(_mm256_mul_pd(tmp, constant), DXtmp);
			_mm256_store_pd(_varDX + j, tmp);
		}
		_testTimer.Stop(3);

		MeasuredRun(1, _rotationSolver->Solve1());
		MeasuredRun(2, CalculateForces());

	}
	_testTimer.Stop(0);
#ifndef NOTIMER
	const int width = 16;
//	_testTimer.SetWidth(width);
	std::cout << "-----------------------------------\n";
	double t1 = _testTimer.Print(1, "Rotations: ");
	double t2 = _testTimer.Print(2, "Forces: ");
	double t3 = _testTimer.Print(3, "Integration: ");
	//_testTimer.Print(5, "Linksh:");
	//std::cout << std::setw(width) << "Summ: " << t1 + t2 + t3 << std::endl;
	_testTimer.Print(0, "Total: ");
#endif
}


void StressStrainCppIterativeSolverFMA::Solve(const int nIterations)
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	FTimer test_timer;
	_iterationNumber = 0;
	_testTimer.Start(0);
	while (_iterationNumber != nIterations && _rotationSolver->IsValid())
	{
		_iterationNumber++;
		_nIteration++;
		
		//Solve1() used only in initialization
		Solve2();
		Solve3();
		Solve4();
		Solve5();
	}
}

// virtual
void StressStrainCppIterativeSolverFMA::InitialSolve()
{
	MeasuredRun(1, _rotationSolver->InitialSolve());
	MeasuredRun(2, CalculateForces());
}


/**
* ������ ������ ������ ������ �����-�����
*/
// virtual
void StressStrainCppIterativeSolverFMA::Solve1()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	MeasuredRun(1, _rotationSolver->Solve1());
	MeasuredRun(2, CalculateForces());
}

/**
* ������ ������ ������ ������ �����-�����
*/
// virtual
void StressStrainCppIterativeSolverFMA::Solve2()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	__m256d timeStep;
	__m256d timeStep2;
	__m256d constantD2;

	timeStep = _mm256_set1_pd(_timeStep);
	timeStep2 = _mm256_set1_pd(_timeStep2);
	constantD2 = _mm256_set1_pd(0.5);


	memcpy(_initX, _varX, sizeof(double)*_nVariables);
	memcpy(_initDX, _varDX, sizeof(double)*_nVariables);

	_testTimer.Start(3);
#pragma omp parallel for num_threads(_numThreads)
	for (int j = 0; j < _nVariables; j += regSize)
	{
		//_hDDX1[j] = _varDDX[j] * _timeStep;
		//_varX[j] += _varDX[j] * _timeStep2;
		//_varDX[j] += _hDDX1[j] * 0.5;

		__m256d Xtmp = _mm256_load_pd(_varX + j);
		__m256d DXtmp = _mm256_load_pd(_varDX + j);
		__m256d DDXtmp = _mm256_load_pd(_varDDX + j);

		__m256d hDDX1 = _mm256_mul_pd(DDXtmp, timeStep);
		_mm256_store_pd(_hDDX1 + j, hDDX1);
		__m256d tmp = _mm256_fmadd_pd(DXtmp, timeStep2, Xtmp);
		//tmp = _mm256_add_pd(_mm256_mul_pd(DXtmp, timeStep2), Xtmp);
		_mm256_store_pd(_varX + j, tmp);
		tmp = _mm256_fmadd_pd(hDDX1, constantD2, DXtmp);
		//tmp = _mm256_add_pd(_mm256_mul_pd(hDDX1, constant), DXtmp);
		_mm256_store_pd(_varDX + j, tmp);
	}
	_testTimer.Stop(3);

	MeasuredRun(1, _rotationSolver->Solve2());
	MeasuredRun(2, CalculateForces());

}

/**
* ������ ������� ������ ������ �����-�����
*/
// virtual
void StressStrainCppIterativeSolverFMA::Solve3()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	__m256d timeStep;
	__m256d timeStep4;
	__m256d constantD2;

	timeStep = _mm256_set1_pd(_timeStep);
	timeStep4 = _mm256_set1_pd(_timeStep4);
	constantD2 = _mm256_set1_pd(0.5);

	_testTimer.Start(3);
#pragma omp parallel for num_threads(_numThreads)
	for (int j = 0; j < _nVariables; j += regSize)
	{
		//_hDDX2[j] = _varDDX[j] * _timeStep;
		//_varX[j] += _hDDX1[j] * _timeStep4;
		//_varDX[j] = _initDX[j] + _hDDX2[j] * 0.5;

		__m256d Xtmp = _mm256_load_pd(_varX + j);
		__m256d DXtmp = _mm256_load_pd(_initDX + j);
		__m256d DDXtmp = _mm256_load_pd(_varDDX + j);
		__m256d hDDX1 = _mm256_load_pd(_hDDX1 + j);
		__m256d hDDX2 = _mm256_mul_pd(DDXtmp, timeStep);

		_mm256_store_pd(_hDDX2 + j, hDDX2);
		__m256d tmp = _mm256_fmadd_pd(hDDX1, timeStep4, Xtmp);
		//tmp = _mm256_add_pd(_mm256_mul_pd(hDDX1, timeStep4), Xtmp);
		_mm256_store_pd(_varX + j, tmp);
		tmp = _mm256_fmadd_pd(hDDX2, constantD2, DXtmp);
		//tmp = _mm256_add_pd(_mm256_mul_pd(hDDX2, constant), DXtmp);
		_mm256_store_pd(_varDX + j, tmp);
	}
	_testTimer.Stop(3);

	MeasuredRun(1, _rotationSolver->Solve3());
	MeasuredRun(2, CalculateForces());
}

/**
* ������ ��������� ������ ������ �����-�����
*/
// virtual
void StressStrainCppIterativeSolverFMA::Solve4()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	__m256d timeStep;
	__m256d constantD2;

	timeStep = _mm256_set1_pd(_timeStep);
	constantD2 = _mm256_set1_pd(0.5);

	_testTimer.Start(3);
#pragma omp parallel for num_threads(_numThreads)
	for (int j = 0; j < _nVariables; j += regSize)
	{
		//_hDDX3[j] = _varDDX[j] * _timeStep;
		//_varX[j] = _initX[j] + (_initDX[j] + _hDDX2[j] * 0.5) * _timeStep;
		//_varDX[j] = _initDX[j] + _hDDX3[j];

		__m256d Xtmp = _mm256_load_pd(_initX + j);
		__m256d DXtmp = _mm256_load_pd(_initDX + j);
		__m256d DDXtmp = _mm256_load_pd(_varDDX + j);
		__m256d hDDX2 = _mm256_load_pd(_hDDX2 + j);

		__m256d hDDX3 = _mm256_mul_pd(DDXtmp, timeStep);
		_mm256_store_pd(_hDDX3 + j, hDDX3);
		__m256d tmp = _mm256_fmadd_pd(hDDX2, constantD2, DXtmp);
		//tmp = _mm256_add_pd(_mm256_mul_pd(hDDX2, constant), DXtmp);
		tmp = _mm256_fmadd_pd(tmp, timeStep, Xtmp);
		//tmp = _mm256_add_pd(_mm256_mul_pd(tmp, timeStep), Xtmp);
		_mm256_store_pd(_varX + j, tmp);
		_mm256_store_pd(_varDX + j, _mm256_add_pd(DXtmp, hDDX3));
	}
	_testTimer.Stop(3);

	MeasuredRun(1, _rotationSolver->Solve4());
	MeasuredRun(2, CalculateForces());
}

/**
* ������ ����� ������ ������ �����-�����
*/
// virtual
void StressStrainCppIterativeSolverFMA::Solve5()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	__m256d timeStep;
	__m256d constantD6;

	timeStep = _mm256_set1_pd(_timeStep);
	constantD6 = _mm256_set1_pd(1 / 6.0);

	_testTimer.Start(3);
#pragma omp parallel for num_threads(_numThreads)
	for (int j = 0; j < _nVariables; j += regSize)
	{
		//float sDDX = _hDDX2[j] + _hDDX3[j];
		//_varX[j] = _initX[j] + (_initDX[j] + sDDX / 6.0) * _timeStep;
		//_varDX[j] = _initDX[j] + (_hDDX1[j] + sDDX + sDDX + _varDDX[j] * _timeStep) / 6.0;

		__m256d Xtmp = _mm256_load_pd(_initX + j);
		__m256d DXtmp = _mm256_load_pd(_initDX + j);
		__m256d DDXtmp = _mm256_load_pd(_varDDX + j);
		__m256d hDDX1 = _mm256_load_pd(_hDDX1 + j);
		__m256d hDDX2 = _mm256_load_pd(_hDDX2 + j);
		__m256d hDDX3 = _mm256_load_pd(_hDDX3 + j);

		__m256d sDDX = _mm256_add_pd(hDDX2, hDDX3);
		__m256d tmp = _mm256_fmadd_pd(sDDX, constantD6, DXtmp);
		//tmp = _mm256_add_pd(_mm256_mul_pd(sDDX, constant), DXtmp);
		tmp = _mm256_fmadd_pd(tmp, timeStep, Xtmp);
		//tmp = _mm256_add_pd(_mm256_mul_pd(tmp, timeStep), Xtmp);
		_mm256_store_pd(_varX + j, tmp);
		tmp = _mm256_add_pd(_mm256_add_pd(hDDX1, sDDX), sDDX);
		tmp = _mm256_fmadd_pd(DDXtmp, timeStep, tmp);
		//tmp = _mm256_add_pd(_mm256_mul_pd(DDXtmp, timeStep), tmp);
		tmp = _mm256_fmadd_pd(tmp, constantD6, DXtmp);
		//tmp = _mm256_add_pd(_mm256_mul_pd(tmp, constant), DXtmp);
		_mm256_store_pd(_varDX + j, tmp);
	}
	_testTimer.Stop(3);

	MeasuredRun(1, _rotationSolver->Solve1());
	MeasuredRun(2, CalculateForces());
}

void StressStrainCppIterativeSolverFMA::CalculateStrains
(
size_t side,			// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
double *shiftStrains,		// ����� ����������
double *velocityStrains,	// ����� ���. ���������
size_t nodeId1,				// ����� ���� 1
size_t nodeId2					// ����� ���� 2
) const
{

	// Start AVX code
	double* pmatA01 = GetRotationMatrix(nodeId1);
	double* pmatA02 = GetRotationMatrix(nodeId2);

	// vecStride must be 4
	__m256d matA01row1 = _mm256_load_pd(pmatA01);
	__m256d matA01row2 = _mm256_load_pd(pmatA01 + vecStride);
	__m256d matA01row3 = _mm256_load_pd(pmatA01 + vecStride2);

	__m256d matA02el1;
	__m256d matA02el2;
	__m256d matA02el3;

	__m256d matA21row1;
	__m256d matA21row2;
	__m256d matA21row3;

	// matA02 column 1/3
	matA02el1 = _mm256_set1_pd(pmatA02[0]);
	matA02el2 = _mm256_set1_pd(pmatA02[vecStride]);
	matA02el3 = _mm256_set1_pd(pmatA02[vecStride2]);

	// matA01.TMul(matA02)
	matA21row1 = _mm256_mul_pd(matA01row1, matA02el1);
	matA21row1 = _mm256_fmadd_pd(matA01row2, matA02el2, matA21row1);
	matA21row1 = _mm256_fmadd_pd(matA01row3, matA02el3, matA21row1);

	// matA02 column 2/3
	matA02el1 = _mm256_set1_pd(pmatA02[1]);
	matA02el2 = _mm256_set1_pd(pmatA02[vecStride + 1]);
	matA02el3 = _mm256_set1_pd(pmatA02[vecStride2 + 1]);

	// matA01.TMul(matA02)
	matA21row2 = _mm256_mul_pd(matA01row1, matA02el1);
	matA21row2 = _mm256_fmadd_pd(matA01row2, matA02el2, matA21row2);
	matA21row2 = _mm256_fmadd_pd(matA01row3, matA02el3, matA21row2);

	// matA02 column 3/3
	matA02el1 = _mm256_set1_pd(pmatA02[2]);
	matA02el2 = _mm256_set1_pd(pmatA02[vecStride + 2]);
	matA02el3 = _mm256_set1_pd(pmatA02[vecStride2 + 2]);

	// matA01.TMul(matA02)
	matA21row3 = _mm256_mul_pd(matA01row1, matA02el1);
	matA21row3 = _mm256_fmadd_pd(matA01row2, matA02el2, matA21row3);
	matA21row3 = _mm256_fmadd_pd(matA01row3, matA02el3, matA21row3);
	// ������� A_{21} ������������

	double*  vc1 = GetRadiusVector(side);
	__m256d ivecC1 = _mm256_load_pd(vc1);
	__m256d vecDP = _mm256_sub_pd(
		_mm256_load_pd(GetElementShift(nodeId1)),
		_mm256_load_pd(GetElementShift(nodeId2))); // P1-P2

	matA02el1 = _mm256_set1_pd(vc1[0]);
	matA02el2 = _mm256_set1_pd(vc1[1]);
	matA02el3 = _mm256_set1_pd(vc1[2]);

	__m256d
		res = _mm256_fmadd_pd(matA21row1, matA02el1, ivecC1);
	res = _mm256_fmadd_pd(matA21row2, matA02el2, res);
	res = _mm256_fmadd_pd(matA21row2, matA02el2, res);

	__declspec(align(32)) double tmp[4];
	_mm256_store_pd(tmp, vecDP);
	matA02el1 = _mm256_set1_pd(tmp[0]);
	matA02el2 = _mm256_set1_pd(tmp[1]);
	matA02el3 = _mm256_set1_pd(tmp[2]);

	res = _mm256_fmadd_pd(matA01row1, matA02el1, res);
	res = _mm256_fmadd_pd(matA01row2, matA02el2, res);
	res = _mm256_fmadd_pd(matA01row3, matA02el3, res);

	_mm256_store_pd(shiftStrains, res); // �������� SL, �������� ����������

	// ������ VL
	// ��������� ������ ������� �������� ��������� ����� ����� �2-�1 � ��1
	// Vec3 VecT1 = matA01.Tmul(vecV1 - vecV2) + vecW1.Cross(vecC1) - matA12*(vecW2.Cross(vecC2));
	// vecC2 = -vecC1
	// vecC2.Cross(vecW1) = -vecC1.Cross(vecW1)

	__declspec(align(32)) double cp1[4] = { 0 };	// ��������� ������������ 
	__declspec(align(32)) double cp2[4] = { 0 };	// ��������� ������������

	CrossProduct(GetElementVelocityAngular(nodeId1), GetRadiusVector(side), cp1);	// [w1 x c1]
	CrossProduct(GetElementVelocityAngular(nodeId2), GetRadiusVector(side), cp2);	// -[w2 x c2] = [w2 x c1]

	res = _mm256_load_pd(&cp1[0]);
	__m256d vecDV = _mm256_sub_pd(
		_mm256_load_pd(GetElementVelocity(nodeId1)),
		_mm256_load_pd(GetElementVelocity(nodeId2))); // V1-V2

	matA02el1 = _mm256_set1_pd(cp2[0]);
	matA02el2 = _mm256_set1_pd(cp2[1]);
	matA02el3 = _mm256_set1_pd(cp2[2]);

	res = _mm256_fmadd_pd(matA21row1, matA02el1, res);
	res = _mm256_fmadd_pd(matA21row2, matA02el2, res);
	res = _mm256_fmadd_pd(matA21row2, matA02el2, res);

	_mm256_store_pd(tmp, vecDV);
	matA02el1 = _mm256_set1_pd(tmp[0]);
	matA02el2 = _mm256_set1_pd(tmp[1]);
	matA02el3 = _mm256_set1_pd(tmp[2]);

	res = _mm256_fmadd_pd(matA01row1, matA02el1, res);
	res = _mm256_fmadd_pd(matA01row2, matA02el2, res);
	res = _mm256_fmadd_pd(matA01row3, matA02el3, res);


	_mm256_store_pd(velocityStrains, res); // �������� VL, �������� ����������

	__m256d x1 = _mm256_load_pd(GetElementShiftAngular(nodeId1));
	__m256d x2 = _mm256_load_pd(GetElementShiftAngular(nodeId2));
	_mm256_store_pd(shiftStrains + vecStride, _mm256_sub_pd(x1, x2));

	x1 = _mm256_load_pd(GetElementVelocityAngular(nodeId1));
	x2 = _mm256_load_pd(GetElementVelocityAngular(nodeId2));
	_mm256_store_pd(velocityStrains + vecStride, _mm256_sub_pd(x1, x2));
}

}

#endif