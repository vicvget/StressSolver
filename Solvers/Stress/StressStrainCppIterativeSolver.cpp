#include "StressStrainCppIterativeSolver.h"
#include "../../AdditionalModules/fmath/Matrix3x3.h"


//#define NOLINKSH
//#define NO_INTOMSUB

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#define SQR(x) ((x) * (x))

using MathHelpers::Vec3;
using MathHelpers::Vec3Ref;
using MathHelpers::MakeVec3;
using MathHelpers::Mat3;

namespace Stress
{

StressStrainCppIterativeSolver::StressStrainCppIterativeSolver
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
		StressStrainCppSolver
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
			),
		_iterationNumber(0)
{
	CalculateRotations();
}

// virtual
StressStrainCppIterativeSolver::~StressStrainCppIterativeSolver()
{

}

// virtual
void StressStrainCppIterativeSolver::Solve
	(
		const int nIterations
	)
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	FTimer test_timer;

	_iterationNumber = 0;
	_testTimer.Start(0);
	while (_iterationNumber != nIterations)
	{
		_iterationNumber++;
		_nIteration++;
		// _timeTmp - это член класса, а тут объявлялась временная переменная с таким же названием
		_timeTmp = _iterationNumber * _timeStep;
		_time = _timeTmp - _timeStep2;

		_stageRK = 1;
		_testTimer.Start(1);
		CalculateRotations();
		_testTimer.Stop(1);
		_testTimer.Start(2);
		CalculateForces();
		_testTimer.Start(2);

		_testTimer.Start(3);
		memcpy(_initX, _varX, sizeof(double)*_nVariables);
		memcpy(_initDX, _varDX, sizeof(double)*_nVariables);

		// RK4 step 1
		#pragma omp parallel for num_threads(_numThreads)
		for (int j = 0; j < _nVariables; j++)
		{
			_hDDX1[j] = _varDDX[j] * _timeStep;
			_varX[j] += _varDX[j] * _timeStep2;
			_varDX[j] += _hDDX1[j] * 0.5;       
		}
		_stageRK = 2;
		_testTimer.Stop(3);


		_testTimer.Start(1);
		CalculateRotations();
		_testTimer.Stop(1);
		_testTimer.Start(2);
		CalculateForces();
		_testTimer.Stop(2);

		_testTimer.Start(3);
		// RK4 step 2
		#pragma omp parallel for num_threads(_numThreads)
		for (int j = 0; j < _nVariables; j++)
		{
			_hDDX2[j] = _varDDX[j] * _timeStep;
			_varX[j] += _hDDX1[j] * _timeStep4;
			_varDX[j] = _initDX[j] + _hDDX2[j] * 0.5;
		}
		_stageRK = 3;
		_testTimer.Stop(3);

		_testTimer.Start(1);
		CalculateRotations();
		_testTimer.Stop(1);
		_testTimer.Start(2);
		CalculateForces();
		_testTimer.Stop(2);
		_time = _timeTmp + _timeStep;

		_testTimer.Start(3);
		// RK4 step 3
		#pragma omp parallel for num_threads(_numThreads)
		for (int j = 0; j < _nVariables; j++)
		{      
			_hDDX3[j] = _varDDX[j] * _timeStep;
			_varX[j] = _initX[j] + (_initDX[j] + _hDDX2[j] * 0.5) * _timeStep;
			_varDX[j] = _initDX[j] + _hDDX3[j];
		}
		_stageRK = 4;
		_testTimer.Stop(3);

		_testTimer.Start(1);
		CalculateRotations();
		_testTimer.Stop(1);
		_testTimer.Start(2);
		CalculateForces();
		_testTimer.Stop(2);

		_testTimer.Start(3);
		// RK4 step 4
		#pragma omp parallel for num_threads(_numThreads)
		for (int j = 0; j < _nVariables; j++)
		{  
			double sDDX = _hDDX2[j]+_hDDX3[j];
			_varX[j] = _initX[j] + (_initDX[j] + sDDX / 6.0) * _timeStep;
			_varDX[j] = _initDX[j] + (_hDDX1[j] + sDDX + sDDX + _varDDX[j] * _timeStep) / 6.0;
		}
		_testTimer.Stop(3);

	}
	_testTimer.Stop(0);
	const int width = 16;
	_testTimer.SetWidth(width);
	std::cout << "-----------------------------------\n";
	double t1 = _testTimer.Print(1, "Rotations: ");
	double t2 = _testTimer.Print(2, "Forces: ");
	double t3 = _testTimer.Print(3, "Integration: ");
	//_testTimer.Print(5, "Linksh:");
	std::cout << std::setw(width) << "Summ: " << t1 + t2 + t3 << std::endl;
	_testTimer.Print(0, "Total: ");

}

/**
* Расчет первой стадии метода Рунге-Кутты
*/
// virtual
void StressStrainCppIterativeSolver::Solve1()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	_iterationNumber++;
	double _timeTmp = _iterationNumber * _timeStep;
	_time = _timeTmp - _timeStep * 0.5;

	_stageRK = 1;
	CalculateRotations();
	CalculateForces();
}

/**
* Расчет второй стадии метода Рунге-Кутты
*/
// virtual
void StressStrainCppIterativeSolver::Solve2()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	memcpy(_initX, _varX, sizeof(double)*_nVariables);
	memcpy(_initDX, _varDX, sizeof(double)*_nVariables);

	// RK4 step 1
	#pragma omp parallel for num_threads(_numThreads)
	for (int j = 0; j < _nVariables; j++)
	{
		_hDDX1[j] = _varDDX[j] * _timeStep;
		_varX[j] += _varDX[j] * _timeStep * 0.5;
		_varDX[j] += _hDDX1[j] * 0.5;       
	}
	_stageRK = 2;
	CalculateRotations();
	CalculateForces();
}

/**
* Расчет третьей стадии метода Рунге-Кутты
*/
// virtual
void StressStrainCppIterativeSolver::Solve3()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	// RK4 step 2
	#pragma omp parallel for num_threads(_numThreads)
	for (int j = 0; j < _nVariables; j++)
	{      
		_hDDX2[j] = _varDDX[j] * _timeStep;
		_varX[j] += _hDDX1[j] * _timeStep * 0.25;
		_varDX[j] = _initDX[j] + _hDDX2[j] * 0.5;
	}
	_stageRK = 3;
	CalculateRotations();
	CalculateForces();
}

/**
* Расчет четвертой стадии метода Рунге-Кутты
*/
// virtual
void StressStrainCppIterativeSolver::Solve4()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	_time = _timeTmp + _timeStep;

	// RK4 step 3
	#pragma omp parallel for num_threads(_numThreads)
	for (int j = 0; j < _nVariables; j++)
	{      
		_hDDX3[j] = _varDDX[j] * _timeStep;
		_varX[j] = _initX[j] + (_initDX[j] + _hDDX2[j] * 0.5) * _timeStep;
		_varDX[j] = _initDX[j] + _hDDX3[j];
	}
	_stageRK = 4;
	CalculateRotations();
	CalculateForces();
}

/**
* Расчет пятой стадии метода Рунге-Кутты
*/
// virtual
void StressStrainCppIterativeSolver::Solve5()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	for (int j = 0; j < _nVariables; j++)
	{  
		double sDDX = _hDDX2[j]+_hDDX3[j];
		_varX[j] = _initX[j] + (_initDX[j] + sDDX / 6.0) * _timeStep;
		_varDX[j] = _initDX[j] + (_hDDX1[j] + sDDX + sDDX + _varDDX[j] * _timeStep) / 6.0;
	}
}

/** Получить смещения
* @param data - массив для записи смещений как скалярного параметра
*/
// virtual
void StressStrainCppIterativeSolver::GetDisplacement
	(
		float* data
	)
{
	int internalIndex = 0;
	double displacement;
	double squareSum;

	for
		(
			int i = 0;
			i < _nElements;
			i++,
			internalIndex += 3
		)
	{		
		squareSum = 0;
		for (int j = 0; j < 3; j++, internalIndex++)
		{
			displacement = _elements[i * 3 + j] - _dataInternal[internalIndex];
			squareSum += SQR(displacement);
		}
		data[i] = (float)sqrt(squareSum);
		//data[i] = tmp;
	}
}

/** Получить напряжения по первой теории прочности
* @param data - массив для записи напряжений как скалярного параметра
*/
// virtual
void StressStrainCppIterativeSolver::GetStressesByFirstTheoryOfStrength
	(
		float* data
	)
{
	double* movements = _dataInternal;
	double coordinates[3];
	double relativeShifts[3];
	double relativeShiftsSigned[3];
	int neighbourNodeNumber;
	int shiftIndex;
	double shift;
	double fullRelativeShift;
	double maxRelativeShift = 0;

	for (int i = 0; i < _nElements; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			coordinates[j] = _dataInternal[6 * i + j];
			relativeShifts[j] = 0;
			relativeShiftsSigned[j] = 0;
		}
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				neighbourNodeNumber = _linkedElements[6 * i + j + 3 * k];
				if (neighbourNodeNumber != 0)
				{
					neighbourNodeNumber--;
					shiftIndex = j;
					shift = movements[6 * neighbourNodeNumber + shiftIndex] - coordinates[shiftIndex];
					if (fabs(shift) > fabs(relativeShifts[shiftIndex]))
					{
						relativeShiftsSigned[shiftIndex] = shift * (2 * k - 1) - _gridStep;
						if (fabs(relativeShiftsSigned[shiftIndex]) < 1e-15)
						{
							relativeShiftsSigned[shiftIndex] = 0.0;
						}
						relativeShifts[shiftIndex] = fabs(shift);
					}
				}
			}
		}
		fullRelativeShift  = relativeShiftsSigned[0];
		for (int j = 1; j < 3; j++)
		{
			if (fabs(relativeShiftsSigned[j]) > fabs(fullRelativeShift))
			{
				fullRelativeShift = relativeShiftsSigned[j];
			}
		}
		data[i] = (float)fabs(_elasticModulus * fullRelativeShift / _gridStep);
	}
}

/** Получить напряжения по von Mises
* @param data - массив для записи напряжений как скалярного параметра
*/
// virtual
void StressStrainCppIterativeSolver::GetStressesByVonMises
	(
		float* data
	)
{
	double* movements = _dataInternal;
	double relativeShiftsSigned[6];
	double maxRelativeShift = 0;

	double sigma[6];
	
	for (int nodeIndex = 0; nodeIndex < _nElements; nodeIndex++)
	{
		for (int i = 0; i < 6; i++)
		{
			relativeShiftsSigned[i] = _stress[6 * nodeIndex + i];
		}
		for (int i = 0; i < 3; i++)
		{
			relativeShiftsSigned[i] /= _gridStep;
		}
		
		double summsq = 0;

		for (int i = 0; i < 6; i++)
		{
			sigma[i] = 0;
			for (int j = 0; j < 6; j++)
			{
				sigma[i] += _lameMatrix[i][j] * relativeShiftsSigned[j];
			}
			if (i > 2)
			{
				sigma[i] *= 2;
				summsq += SQR(sigma[i]);
			}
		}
		//summsq = 0;
		data[nodeIndex] = (float)(_elasticModulus * sqrt
			(
				0.5 *
					(
						SQR(sigma[0] - sigma[1]) +
						SQR(sigma[1] - sigma[2]) +
						SQR(sigma[0] - sigma[2]) +
						6 * summsq
					)
			));
	}
}
#define NEW_FL
#ifdef NEW_FL
void StressStrainCppIterativeSolver::CalculateForces()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	static int it = 0;

	if (_isFirstIteration)
	{
		_isFirstIteration = false;
		_stageRK = 0;
		CalculateRotations();
		_stageRK = 1;
	}

	//double E = _elasticModulusScaled * _gridStep;
	double elasticModulus = _elasticModulusScaled;

	__declspec(align(32)) double strains[8], velocityStrains[8];

	for (size_t elementId1 = 0; elementId1 < _nElements; elementId1++)
	{
		
		double* accelerationVector = GetElementAcceleration(elementId1);

		memset(GetElementAcceleration(elementId1), 0, sizeof(double)*vecStride2);
		memset(GetElementStress(elementId1), 0, sizeof(double)*vecStride2);

		_testTimer.Start(5);

		// обход 6 связанных элементов x-,y-,z-,x+,y+,z+
		for (size_t side = 0; side < 6; side++)
		{
			size_t elementId2 = _linkedElements[side];
			if (elementId2)
			{
				elementId2--;
#ifdef ALIGNED_MEM
				CalculateStrainsAVX(side, strains, velocityStrains, elementId1, elementId2);
#endif
				CalculateStrains(side, strains, velocityStrains, elementId1, elementId2);

				Vec3Ref linear_strains = MakeVec3(&strains[0]);
				Vec3Ref angular_strains = MakeVec3(&strains[0] + vecStride);
				Vec3Ref linear_vstrains = MakeVec3(&velocityStrains[0]);
				Vec3Ref angular_vstrains = MakeVec3(&velocityStrains[0] + vecStride);

				for (size_t component = 0; component < 3; component++)
				{
					GetElementStress(elementId1)[component] += strains[component];
					GetElementStress(elementId1)[component + vecStride] += strains[component + vecStride];
				}

				// сила и момент из полученных деформаций
				Vec3 force = -linear_vstrains * _dampingFactorLinear - linear_strains * elasticModulus * _gridStep;
				Vec3 torque = -angular_vstrains * _dampingFactorAngular - angular_strains * elasticModulus * _gridStep3;

				// расчет суммарных сил
				Mat3 matA01(GetRotationMatrix(elementId1));
				Vec3Ref vR = MakeVec3(GetRadiusVector(side));
				Vec3 vAcc = matA01*force;
				Vec3 vM = vR.Cross(force) + torque;

				for (int i = 0; i < 3; i++)
				{
					accelerationVector[i] += vAcc[i];
					accelerationVector[i + vecStride] += vM[i];
				}

			}
		}
		_testTimer.Stop(5);

		//Vec3Ref force = MakeVec3(&forces[0]);
		//Vec3Ref torque = MakeVec3(&forces[0] + vecStride);
		
	}
	ApplyBoundary(); // модифицирует силы и моменты
	ApplyMass();	 // вычисляет ускорения делением сил на массы и моментов на моменты инерции
}
#else
void StressStrainCppIterativeSolver::pravsubfl() 
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	static int it = 0;

	if (_isFirstIteration)
	{
		_isFirstIteration = false;
		_stageRK = 0;
		std::cout << "#############################################" << std::endl;
		std::cout << "pravsubfl _stageRK == 0 isFirstIteration" << std::endl;
		CalculateRotations();
		_stageRK = 1;
	}

	//double E = _elasticModulusScaled * _gridStep;
	double elasticModulus = _elasticModulusScaled;


	double sigma[6];

	int i, j, k, l, N2;
	double DX, DY, DZ, SILA[6];

	double X1[3], X2[3], SL[6], VL[6], A[36], C[36];

	int nLinks = 0;
	//memset(A,0,sizeof(double)*36);
	//memset(C,0,sizeof(double)*36);
	//memset(_mtx.ia,0,sizeof(int)*(_mtx._size+1));

//#pragma omp parallel for private \
//	( \
//		NADU, \
//		j, \
//		i, \
//		k, \
//		l, \
//		DX, \
//		DY, \
//		DZ, \
//		X1, \
//		X2, \
//		SL, \
//		VL, \
//		A, \
//		C, \
//		N2, \
//		SILA, \
//		sigma, \
//		stressStrainMatrix \
//	) \
//	num_threads(_numThreads)

	for (j = 0; j < _nElements; j++)
	{
		int accOffset = _nElements * vecStride2 * 2 + j * _vecStride2;

		for (i = 0; i < vecStride2; i++) // можно вынести
			_dataInternal[accOffset + i] = 0.0;

		for (i = 0; i < 6; i++) // можно вынести
			_stress[6 * j + i] = 0;

		for (i = 0; i < 6; i++)
		{
			N2 = _linkedElements[j * 6 + i] - 1;

			if (N2 >= 0)
			{
				nLinks++;
				DX = (_elements[j * 3] - _elements[N2 * 3]) * 0.5;
				DY = (_elements[j * 3 + 1] - _elements[N2 * 3 + 1]) * 0.5;
				DZ = (_elements[j * 3 + 2] - _elements[N2 * 3 + 2]) * 0.5;

				X1[0] =- DX;
				X1[1] =- DY;
				X1[2] =- DZ;

				X2[0] = DX;
				X2[1] = DY;
				X2[2] = DZ;

				/*
				if((j == 4) && (N2 == 11))
				{
				int k = 0;
				}
				*/
#ifdef NOLINKSH
				for(int i = 0; i < 6; i++)
				{
					SL[i]=0;
					VL[i]=0;
				}
#else
				_testTimer.Start(5);
				linksh3(X1, X2, SL, VL, A, C, j, N2, _nElements);
				//linksh(X1, X2, SL, VL, A, C, j, N2, _nElements);
				//links_timeStep2(X1, X2, SL, VL, A, C, j, N2, _nElements);
				
				//std::cout << "iter xx = " << _nIteration << " j=" << j << " N2=" << N2 << std::endl;
				//for(int ii = 0; ii < 6; ii++)
				//{
				//	std::cout << std::setw(12) << std::setprecision(4) << SL[ii] << ' ';
				//}
				//std::cout << std::endl;
				//for(int ii = 0; ii < 6; ii++)
				//{
				//	std::cout << std::setw(12) << std::setprecision(4) << VL[ii] << ' ';
				//}
				//std::cout << std::endl << std::endl;

//#ifdef _DEBUG
//				double SL1[6], VL1[6], A1[36], C1[36];
//				linksh(X1, X2, SL1, VL1, A1, C1, j, N2, _nElements);
//				double eps = 1e-7;
//				int debug, id;
//				for(int iI = 0; iI < 6; iI++)
//				{
//					if(fabs(SL1[iI]-SL[iI]) > eps)
//						debug = 1;
//					if(fabs(VL1[iI]-VL[iI]) > eps)
//						debug = 2;
//					for(int jJ = 0; jJ < 6; jJ++)
//					{
//						id = iI*6+jJ;
//						if(fabs(A[iI*6+jJ]-A1[iI*6+jJ]) > eps)
//							debug = 3;
//						if(fabs(C[iI*6+jJ]-C1[iI*6+jJ]) > eps)
//							debug = 4;
//					}
//				}
//#endif
				_testTimer.Stop(5);
				//Genmatl(j,N2,A,C);
				//std::cout << "Node: " << j << " EOF Genmatl\n";
#endif

#ifdef NOLINKSH
				for (k = 0; k < 6; k++)
				{
					if(k == 0)
					{
						SL[k]=(fabs(GR1[j*6+k]-GR1[(N2-1)*6+k])-_gridStep);
						if(j < (N2-1))
							SL[k] = -SL[k];
						VL[k]=(GR1[j*6+_nElements*6+k]-GR1[(N2-1)*6+_nElements*6+k]);
					}
					else
					{
						SL[k]=VL[k]=0;
					}
				}
#endif

				for (k = 0; k < 6; k++)
				{
					_stress[6 * j + k] += SL[k];
				}

				for (k = 0; k < 6; k++)
				{
					/*sigma[k] = 0;
					for (l = 0; l < 6; l++)
					{
						sigma[k] += stressStrainMatrix[k][l] * SL[l];
					}*/
					sigma[k] = SL[k];
				}

				for (k = 0; k < 6; k++)
				{
					//SILA[k] =- (sigma[k] * E + VL[k] * DM);
					SILA[k] = 0;
					if (k < 3)
					{
						SILA[k] = - VL[k] * _dampingFactorLinear;
					}
					else
					{
						SILA[k] = - VL[k] * _dampingFactorAngular;
					}
					if (k < 3)
					{
						SILA[k] -= sigma[k] * elasticModulus * _gridStep; // сила
					}
					else
					{
						SILA[k] -= sigma[k] * elasticModulus * _gridStep3; // момент
					}
				}
				// a = A*F
				for (k = 0; k < 6; k++)
				{
					for (l = 0; l < 3; l++)
					{
						double a1 = A[6 * k + l];
						double a2 = (k < 3) ? -A[6 * k + l + 3] : A[6 * k + l + 3];
						_dataInternal[accOffset + l] += SILA[k] * a1;
						_dataInternal[accOffset + l + vecStride] += SILA[k] * a2;
					}
				}
				int dbg = 1;
			}
		}
	}
	ApplyBoundary();
	ApplyMass();
}
#endif

void StressStrainCppIterativeSolver::ApplyBoundary()
{
	//работает с типами из frm_provider
	//BCT_stresstrainBoundaryForce = 3,
	//BCT_stresstrainBoundarySealing = 4,

	double* dataPointer = _dataInternal + _nElements * vecStride2 * 2;
	vector<BoundaryParams>::iterator it = _boundaryParamsSet.begin();

	while (it != _boundaryParamsSet.end())
	{
		switch (it->GetKind())
		{
		case 3:
			it->ApplyForceBoundary(dataPointer);
			break;

		case 4:
			it->ApplySealedBoundary(dataPointer);
			break;

		default:
			break;
		}
		it++;
	}
}

void StressStrainCppIterativeSolver::ApplyMass()
{
#ifdef _DEBUG
	//_controlfp(0, EM_ZERODIVIDE);
	//_control87(~_EM_ZERODIVIDE, _MCW_EM);
#endif
	double* accelerations = _dataInternal + _nElements * vecStride2 * 2;

	for (int i = 0; i < _nElements; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			accelerations[vecStride2 * i + j] /= _cellMass;
			accelerations[vecStride2 * i + j + vecStride] /= _cellMass * _gridStep * _gridStep / 6;
		}
	}
}
}