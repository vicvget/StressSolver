#include "StressStrainFortranIterativeSolver.h"

#include "../../Fcore/Exceptions/fcExceptions.h"


#include <fstream>


using std::ofstream;


//#define NOLINKSH
//#define NO_INTOMSUB
#define TIMER

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#define SQR(x) ((x) * (x))


StressStrainFortranIterativeSolver::StressStrainFortranIterativeSolver
	(
		double* params,
		int* links,
		int nLinks,
		double *nodes,
		int nNodes,
		double gridStep,
		double timeStep,
		int numThreads,
		int stride
	)
	:
		StressStrainFortranSolver
			(
				params,
				links,
				nLinks,
				nodes,
				nNodes,
				gridStep,
				timeStep,
				numThreads,
				stride
			),
		_iterationNumber(0)
{
}

// virtual
void StressStrainFortranIterativeSolver::Solve
	(
		const int nIterations
	)
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	//SetZeroVelocities();
	FTimer test_timer;
	//_testTimer.Allocate(4);
	double H2 = H * 0.5;
	double H4 = H * 0.25;

	_iterationNumber = 0;
#ifdef TIMER
	_testTimer.Start(0);
#endif

	//if (_nIteration == 0)
	//{
	//	CreateSystemOfLinearEquationsForStiffness();
	//	if (_readIco)
	//	{
	//		// TODO: readico
	//	}
	//}

	while (_iterationNumber != nIterations)
	{
		//		if (_iterationNumber == 0)
		//		CreateStiffnessMatrix("WTF.txt");
		_iterationNumber++;
		_nIteration++;
		double TZ = _iterationNumber * H;
		T = TZ - H2;

		METS = 1;
#ifdef TIMER
		_testTimer.Start(1);
#endif		
		intomsub();
#ifdef TIMER
		_testTimer.Stop(1);
		_testTimer.Start(2);
#endif		
		pravsubfl();
#ifdef TIMER
		_testTimer.Start(2);
		_testTimer.Start(3);
#endif		
		memcpy(_initX, _varX, sizeof(double)*_nVariables);
		memcpy(_initDX, _varDX, sizeof(double)*_nVariables);

		// RK4 step 1
		#pragma omp parallel for num_threads(NumThreads)
		for (int j = 0; j < _nVariables; j++)
		{
			_hDDX1[j] = _varDDX[j] * H;
			_varX[j] += _varDX[j] * H2;
			_varDX[j] += _hDDX1[j] * 0.5;
		}
		METS = 2;

#ifdef TIMER
		_testTimer.Stop(3);
		_testTimer.Start(1);
#endif
		intomsub();

#ifdef TIMER
		_testTimer.Stop(1);
		_testTimer.Start(2);
#endif

		pravsubfl();

#ifdef TIMER
		_testTimer.Stop(2);
		_testTimer.Start(3);
#endif		// RK4 step 2
		#pragma omp parallel for num_threads(NumThreads)
		for (int j = 0; j < _nVariables; j++)
		{
			_hDDX2[j] = _varDDX[j] * H;
			_varX[j] += _hDDX1[j] * H4;
			_varDX[j] = _initDX[j] + _hDDX2[j] * 0.5;
		}
		METS = 3;

#ifdef TIMER
		_testTimer.Stop(3);
		_testTimer.Start(1);
#endif
		intomsub();

#ifdef TIMER
		_testTimer.Stop(1);
		_testTimer.Start(2);
#endif
		pravsubfl();
#ifdef TIMER
		_testTimer.Stop(2);
		_testTimer.Start(3);
#endif
		T = TZ + H;
		// RK4 step 3
		#pragma omp parallel for num_threads(NumThreads)
		for (int j = 0; j < _nVariables; j++)
		{
			_hDDX3[j] = _varDDX[j] * H;
			_varX[j] = _initX[j] + (_initDX[j] + _hDDX2[j] * 0.5) * H;
			_varDX[j] = _initDX[j] + _hDDX3[j];
		}
		METS = 4;
#ifdef TIMER
		_testTimer.Stop(3);
		_testTimer.Start(1);
#endif
		intomsub();
#ifdef TIMER
		_testTimer.Stop(1);
		_testTimer.Start(2);
#endif
		pravsubfl();
#ifdef TIMER
		_testTimer.Stop(2);
		_testTimer.Start(3);
#endif
		// RK4 step 4
		#pragma omp parallel for num_threads(NumThreads)
		for (int j = 0; j < _nVariables; j++)
		{
			double sDDX = _hDDX2[j] + _hDDX3[j];
			_varX[j] = _initX[j] + (_initDX[j] + (_hDDX1[j] + _hDDX2[j] + _hDDX3[j]) / 6.0) * H;
			_varDX[j] = _initDX[j] + (_hDDX1[j] + sDDX + sDDX + _varDDX[j] * H) / 6.0;
		}
#ifdef TIMER
		_testTimer.Stop(3);
#endif
		if (_writeIco)
		if (_nIteration == _nWriteIteration)
		{
			// TODO: writeico
		}
	}
#ifdef TIMER
	_testTimer.Stop(0);
	_testTimer.Print(0, "Total:");
	double t1 = _testTimer.Print(1, "Intomsub:");
	double t2 = _testTimer.Print(2, "Pravsubfl:");
	double t3 = _testTimer.Print(3, "Integration:");
	_testTimer.Print(5, "Linksh:");
	std::cout << "Sum: " << t1+t2+t3 << std::endl;
#endif

}

/**
* Расчет первой стадии метода Рунге-Кутты
*/
// virtual
void StressStrainFortranIterativeSolver::Solve1()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	_iterationNumber++;
	double TZ = _iterationNumber * H;
	T = TZ - H * 0.5;

	METS = 1;
	intomsub();
	pravsubfl();
}

/**
* Расчет второй стадии метода Рунге-Кутты
*/
// virtual
void StressStrainFortranIterativeSolver::Solve2()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	memcpy(_initX, _varX, sizeof(double)*_nVariables);
	memcpy(_initDX, _varDX, sizeof(double)*_nVariables);

	// RK4 step 1
	#pragma omp parallel for num_threads(NumThreads)
	for (int j = 0; j < _nVariables; j++)
	{
		_hDDX1[j] = _varDDX[j] * H;
		_varX[j] += _varDX[j] * H * 0.5;
		_varDX[j] += _hDDX1[j] * 0.5;
	}
	METS = 2;
	intomsub();
	pravsubfl();
}

/**
* Расчет третьей стадии метода Рунге-Кутты
*/
// virtual
void StressStrainFortranIterativeSolver::Solve3()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	// RK4 step 2
	#pragma omp parallel for num_threads(NumThreads)
	for (int j = 0; j < _nVariables; j++)
	{
		_hDDX2[j] = _varDDX[j] * H;
		_varX[j] += _hDDX1[j] * H * 0.25;
		_varDX[j] = _initDX[j] + _hDDX2[j] * 0.5;
	}
	METS = 3;
	intomsub();
	pravsubfl();
}

/**
* Расчет четвертой стадии метода Рунге-Кутты
*/
// virtual
void StressStrainFortranIterativeSolver::Solve4()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	T = TZ + H;

	// RK4 step 3
	#pragma omp parallel for num_threads(NumThreads)
	for (int j = 0; j < _nVariables; j++)
	{
		_hDDX3[j] = _varDDX[j] * H;
		_varX[j] = _initX[j] + (_initDX[j] + _hDDX2[j] * 0.5) * H;
		_varDX[j] = _initDX[j] + _hDDX3[j];
	}
	METS = 4;
	intomsub();
	pravsubfl();
}

/**
* Расчет пятой стадии метода Рунге-Кутты
*/
// virtual
void StressStrainFortranIterativeSolver::Solve5()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	for (int j = 0; j < _nVariables; j++)
	{
		double sDDX = _hDDX2[j] + _hDDX3[j];
		_varX[j] = _initX[j] + (_initDX[j] + sDDX / 6.0) * H;
		_varDX[j] = _initDX[j] + (_hDDX1[j] + sDDX + sDDX + _varDDX[j] * H) / 6.0;
	}
}

//void StressStrainFortranIterativeSolver::urejl4s(double* a, double *ug, double *om, int mets, const int id)
//{
//	int K=id*3;
//	int KA=id*9;
//	
//	if (mets==0) 
//	{
//		_rotationSolver->Update(K);
//		_rotationSolver->FLAGRY[id]=0;
//		_copym->Copy(a,_rotationSolver->A1+KA,id,T);
//		
//		return;
//	}
//	
//	double UY=_rotationSolver->_varX[K];
//	UY=fabs(DMOD_c(UY,2*M_PI));           // DMOD - ???????
//
//	if ((mets==1) && (UY>1.28))
//	{
//
//		_rotationSolver->Update(K);
//		_rotationSolver->FLAGRY[id]=1.0;
//		_rotationSolver->IFLAGRY=1;
//		
//		_copym->Copy(a,_rotationSolver->A1+KA,id,T);
//	} 
//	_rotationSolver->Update(K,mets);
//	_rotationSolver->UpdateR(K,om,H);
//	_rotationSolver->UpdateR2(K,mets);
//
//	_rotationSolver->UpdateMtx(K, a);
//	MatrixMul(_rotationSolver->A1+KA,a);
//}
//
//void StressStrainFortranIterativeSolver::MatrixMul(double *a1,double *a2)
//{
//	double az[9];
//
//	for (int i=0;i<9;i++)
//		az[i]=a2[i];
//
//	for (int j=0;j<3;j++)
//		for (int i=0;i<3;i++)  
//			a2[j*3+i]=a1[j*3]*az[i]+a1[j*3+1]*az[3+i]+a1[j*3+2]*az[6+i];
//}


/** Получить смещения
* @param data - массив для записи смещений как скалярного параметра
*/
// virtual
void StressStrainFortranIterativeSolver::GetDisplacement
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
			displacement = _nodes[i * 3 + j] - _dataInternal[internalIndex];
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
void StressStrainFortranIterativeSolver::GetStressesByFirstTheoryOfStrength
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
				neighbourNodeNumber = MLINK[6 * i + j + 3 * k];
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
		fullRelativeShift = relativeShiftsSigned[0];
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
void StressStrainFortranIterativeSolver::GetStressesByVonMises
	(
		float* data
	)
{
	double* movements = _dataInternal;
	//double coordinates[6];
	//double relativeShifts[6];
	double relativeShiftsSigned[6];
	//int neighbourNodeNumber;
	//double shift;
	double maxRelativeShift = 0;

	double mtx[6][6];
	double sigma[6];

	FindStressStrainMatrix(mtx);
	for (int nodeIndex = 0; nodeIndex < _nElements; nodeIndex++)
	{
		/*for (int j = 0; j < 6; j++)
		{
			coordinates[j] = _dataInternal[6 * i + j];
			relativeShifts[j] = 0;
			relativeShiftsSigned[j] = 0;
		}
		for (int j = 0; j < 3; j++)
		{
			int countCoordinateNeighbour = 0;
			for (int k = 0; k < 2; k++)
			{
				neighbourNodeNumber = MLINK[6 * i + j + 3 * k];
				if (neighbourNodeNumber != 0)
				{
					countCoordinateNeighbour++;
					neighbourNodeNumber--;
					shift = movements[6 * neighbourNodeNumber + j] - coordinates[j];
					relativeShiftsSigned[j] += shift * (2 * k - 1);
					relativeShiftsSigned[j] -= _gridStep;
					if (fabs(relativeShiftsSigned[j]) < 1e-15)
					{
						relativeShiftsSigned[j] = 0;
					}
				}
			}
			if (countCoordinateNeighbour == 2)
			{
				relativeShiftsSigned[j] *= 0.5;
			}
		}
		for (int j = 0; j < 3; j++)
		{
			int countCoordinateNeighbour = 0;
			for (int l = 0; l < 3; l++)
			{
				if (l == 1)
				{
					for (int k = 0; k < 2; k++)
					{
						neighbourNodeNumber = MLINK[6 * i + j + 3 * k];
						if (neighbourNodeNumber != 0)
						{
							countCoordinateNeighbour++;
							neighbourNodeNumber--;
							shift = movements[6 * neighbourNodeNumber + 3 + l] - coordinates[3 + l];
							relativeShiftsSigned[3 + l] += shift * (2 * k - 1);
							if (fabs(relativeShiftsSigned[3 + l]) < 1e-15)
							{
								relativeShiftsSigned[3 + l] = 0;
							}
						}
					}
				}
			}
			if (countCoordinateNeighbour > 0)
			{
				relativeShiftsSigned[j] /= countCoordinateNeighbour;
			}
		}
		for (int i = 0; i < 6; i++)
		{
			relativeShiftsSigned[i] = (i < 3) ? relativeShiftsSigned[i] / _gridStep : tan(relativeShiftsSigned[i] / 2);
		}*/
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
				sigma[i] += mtx[i][j] * relativeShiftsSigned[j];
			}
			if (i > 2)
			{
				sigma[i] *= 2;
				summsq += SQR(sigma[i]);
			}
		}
		//summsq = 0;
		data[nodeIndex] = static_cast<float>
			(
				_elasticModulus * sqrt
					(
						0.5 *
							(
								SQR(sigma[0] - sigma[1]) +
								SQR(sigma[1] - sigma[2]) +
								SQR(sigma[0] - sigma[2]) +
								6.0 * summsq
							)
					)
			);
		if (data[nodeIndex] > 1e6)
		{
			int debug = 0;
		}
	}
}

// пересчет второй производной температуры
void StressStrainFortranIterativeSolver::pravsubfl()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	double* GR1 = _dataInternal;
	static int it = 0;

	if (_isFirstIteration)
	{
		_isFirstIteration = false;
		METS = 0;
		std::cout << "#############################################" << std::endl;
		std::cout << "pravsubfl METS == 0 isFirstIteration" << std::endl;
		intomsub();
		METS = 1;
	}

	//double E = _elasticModulusScaled * _gridStep;
	double elasticModulus = _elasticModulusScaled;
	double G = _elasticModulusScaled / 2;
	double dampingFactorForward = 0.01 * 2 * _dampingFactor * sqrt(_elasticModulusScaled * _cellMass * _gridStep);
	double dampingFactorRotatory = 10 * _dampingFactor * sqrt(2 * _elasticModulusScaled * _cellMass * _gridStep / 3) *
		_gridStep * _gridStep;
	double sigma[6];
	double stressStrainMatrix[6][6];
	//double PUAS=0;
	//double DH=0.1;

	FindStressStrainMatrix(stressStrainMatrix);

	int NADU;
	int i, j, k, l, N2;
	double DX, DY, DZ, SILA[6];

	double X1[3], X2[3], SL[6], VL[6], A[36], C[36];
	//double SL1[6], VL1[6], A1[36], C1[36];

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
	//	num_threads(NumThreads)
	for (j = 0; j < _nElements; j++)
	{
		// NAD=j*6;
		NADU = _nElements * 2 * 6 + j * 6;

		for (i = 0; i < 6; i++) // можно вынести
		{
			GR1[NADU + i] = 0.0;
			_stress[6 * j + i] = 0;
		}

		for (i = 0; i < 6; i++)
		{
			N2 = MLINK[j * 6 + i] - 1;

			if (N2 >= 0)
			{
				nLinks++;
				DX = (_nodes[j * 3] - _nodes[N2 * 3]) * 0.5;
				DY = (_nodes[j * 3 + 1] - _nodes[N2 * 3 + 1]) * 0.5;
				DZ = (_nodes[j * 3 + 2] - _nodes[N2 * 3 + 2]) * 0.5;

				X1[0] = -DX;
				X1[1] = -DY;
				X1[2] = -DZ;

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
				for (int i = 0; i < 6; i++)
				{
					SL[i] = 0;
					VL[i] = 0;
				}
#else
#ifdef TIMER
				_testTimer.Start(5);
#endif		
				linksh3(X1, X2, SL, VL, A, C, j, N2, _nElements);
				//linksh(X1, X2, SL, VL, A, C, j, N2, _nElements);
				//linksh2(X1, X2, SL, VL, A, C, j, N2, _nElements);


				//std::cout << "iter = " << _nIteration << " j=" << j << " N2=" << N2 << std::endl;
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
#ifdef TIMER
				_testTimer.Stop(5);
#endif
				//Genmatl(j,N2,A,C);
				//std::cout << "Node: " << j << " EOF Genmatl\n";
#endif

				for (k = 0; k < 6; k++)
				{
#ifdef NOLINKSH
					if (k == 0)
					{
						SL[k] = (fabs(GR1[j * 6 + k] - GR1[(N2 - 1) * 6 + k]) - _gridStep);
						if (j < (N2 - 1))
							SL[k] = -SL[k];
						VL[k] = (GR1[j * 6 + _nElements * 6 + k] - GR1[(N2 - 1) * 6 + _nElements * 6 + k]);
					}
					else
					{
						SL[k] = VL[k] = 0;
					}
#endif
				}

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
						SILA[k] = -VL[k] * dampingFactorForward;
					}
					else
					{
						SILA[k] = -VL[k] * dampingFactorRotatory;
					}
					if (k < 3)
					{
						SILA[k] -= sigma[k] * elasticModulus * _gridStep; // сила
					}
					else
					{
						//SILA[k] -= sigma[k] * E * _gridStep;
						SILA[k] -= sigma[k] * elasticModulus * _gridStep * _gridStep * _gridStep; // сила по напряжению
						//SILA[k] -= SL[k] * E * _gridStep; // ???
					}
				}
				for (k = 0; k < 6; k++)
				{
					for (l = 0; l < 6; l++)
					{
						double a = A[6 * k + l];

						if ((k < 3) && (l > 2))
						{
							a = -a;
						}
						//GR1[NADU + l] += SILA[k] * A[6 * k + l];
						GR1[NADU + l] += SILA[k] * a;
					}
				}
				int dbg = 1;
			}
		}
	}
	//std::cout << "            nLinks = " << nLinks << std::endl;

	//for(int i = 0; i < _nElements*6; i++)
	//{
	//	//for(int j = 0; j < 6; j++)
	//	//{
	//		std::cout << "\nEQUATION: " << i+1 << std::endl;
	//		std::cout << "\nADDRES: ";
	//		for(int k = _mtx.ia[i]; k < _mtx.ia[i+1]; k++)
	//		{
	//			std::cout << _mtx.ja[k-1] << ' ';				
	//		}
	//		std::cout << "\nELEM: ";
	//		std::cout << std::endl;
	//		for(int k = _mtx.ia[i]; k < _mtx.ia[i+1]; k++)
	//		{
	//			std::cout << _mtx.a[k-1] << ' ';
	//		}
	//		std::cout << std::endl;
	////	}
	//}

	ApplyBoundary();
	ApplyMass();
}

void StressStrainFortranIterativeSolver::FindStressStrainMatrix
	(
		double (&matrix)[6][6]
	)
{
	const double lambda = _poissonRatio / (1.0 - SQR(_poissonRatio));
	const double mu = 0.5 / (1.0 + _poissonRatio);

	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			matrix[i][j] = 0;
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			matrix[i][j] = lambda;
			if (i == j)
			{
				matrix[i][j] += mu * 2.0;
			}
		}
	}
	for (int i = 3; i < 6; i++)
	{
		matrix[i][i] = mu;
	}
}

void StressStrainFortranIterativeSolver::ApplyBoundary()
{
	//работает с типами из frm_provider
	//BCT_stresstrainBoundaryForce = 3,
	//BCT_stresstrainBoundarySealing = 4,

	double* dataPointer = _dataInternal + _nElements * 6 * 2;
	vector<BoundaryParams>::iterator it = _boundaryParamsSet.begin();

	while (it != _boundaryParamsSet.end())
	{
		//std::cout << "Boundary nodes count: " << it->GetNodesCount() << std::endl;
		switch (it->GetKind())
		{
		case 3:
			it->ApplyForceBoundary(dataPointer);
			/*
			it->CorrectForceBoundary
				(
					_dataInternal + _nElements*6*2,
					_dataInternal + _nElements*6,
					_dataInternal,
					_nodes,
					_elasticModulusScaled,
					_dampingFactor
				);
			*/
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

void StressStrainFortranIterativeSolver::ApplyMass()
{
#ifdef _DEBUG
	//_controlfp(0, EM_ZERODIVIDE);
	//_control87(~_EM_ZERODIVIDE, _MCW_EM);
#endif
	double* accelerations = _dataInternal + _nElements * 6 * 2;

	for (int i = 0; i < _nElements; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			accelerations[6 * i + j] /= _cellMass;
		}
		for (int j = 3; j < 6; j++)
		{
			accelerations[6 * i + j] /= _cellMass * _gridStep * _gridStep / 6;
		}
	}
}

//void StressStrainFortranIterativeSolver::AddBoundary(int* boundaryNodesIndices, 
//	int numberOfBoundaryNodes, int bcKind, double* bcParams)
//{
//	BoundaryParams bp(bcKind, bcParams, boundaryNodesIndices, numberOfBoundaryNodes);
//	_boundaryParamsSet.push_back(bp);
//}
//
//void StressStrainFortranIterativeSolver::ChangeBoundaryParam(const int bcNumber,
//	const int bcParamNumber, const double bcParamValue)
//{
//	_boundaryParamsSet[bcNumber].Params[bcParamNumber] = bcParamValue / _stiffScale;
//}
//
//double StressStrainFortranIterativeSolver::GetBoundaryParam(const int bcNumber,
//	const int bcParamNumber)
//{
//	return _boundaryParamsSet[bcNumber].Params[bcParamNumber];
//}
//

void DumpRhs(const double* rhs, int count, const string& fileNameText, const string& fileNameBinary)
{
	ofstream ofs1(fileNameBinary, ofstream::binary);
	ofstream ofs2(fileNameText);
	if (ofs1.is_open())
	{
		ofs1.write(reinterpret_cast<const char*>(rhs), sizeof(double)*count);
		ofs1.close();
	}
	else
	{
		std::cout << "ERROR OPEN FILE\n";
	}
#ifdef DUMPTEXT
	if (ofs2.is_open())
	{
		for (int i = 0; i < count; i++)
		{
			ofs2 << rhs[i] << std::endl;
		}
		ofs2.close();
	}
	else
	{
		std::cout << "ERROR OPEN FILE\n";
	}
#endif
}

/**
* Сформировать СЛАУ на основе матрицы жесткости
* @param stiffnessMatrixFileName - общее наименование файлов для записи формируемой матрицы жесткости
* в ряде форматов
* @param stiffnessMatrixFileFormats - список форматов файлов, в которые будет записана
* формируемая матрица жесткости
* @param stiffnessRHSFileName - общее наименование файлов для записи формируемого вектора
* правых частей СЛАУ на основе матрицы жесткости
* @param writeStiffnessRHSVectorToTextFileFlag - флаг, следует ли записывать формируемый вектор
* правых частей СЛАУ на основе матрицы жесткости в текстовый файл
* @param stiffnessEpsilon - число, использующееся в численном нахождении производной
* при формировании матрицы жесткости
*/
void StressStrainFortranIterativeSolver::CreateSystemOfLinearEquationsForStiffness
	(
		const string& stiffnessMatrixFileName,
		StiffnessMatrixFileFormats stiffnessMatrixFileFormats,
		const string& stiffnessRHSFileName,
		bool writeStiffnessRHSVectorToTextFileFlag,
		double stiffnessEpsilon
	)
try
{
	const IsSealedFlagsList isSealedFlagsList = FormIsSealedFlagsList();

	CreateStiffnessMatrix
		(
			isSealedFlagsList,
			stiffnessMatrixFileName,
			stiffnessMatrixFileFormats,
			stiffnessEpsilon
		);
	CreateStiffnessRHSVector
		(
			isSealedFlagsList,
			stiffnessRHSFileName,
			writeStiffnessRHSVectorToTextFileFlag
		);
}
catch (const exceptions::CoreException& exception)
{
	std::cout << exception.ToString() << std::endl;
}

/**
* Сформировать список наборов флагов закрепленных степеней свободы элементов сеточного представления
* @return формируемый список наборов флагов закрепленных степеней свободы элементов сеточного представления
*/
IsSealedFlagsList StressStrainFortranIterativeSolver::FormIsSealedFlagsList() const
{
	IsSealedFlagsList isSealedFlagsList(_nElements);

	for (const auto& boundaryParams : _boundaryParamsSet)
	{
		if (boundaryParams.GetKind() == 4)
		{
			boundaryParams.FormIsSealedFlagsList(isSealedFlagsList);
		}
	}

	return isSealedFlagsList;
}

/**
* Сформировать матрицу жесткости
* @param isSealedFlagsList - список наборов флагов закрепленных степеней свободы
* элементов сеточного представления
* @param stiffnessMatrixFileName - общее наименование файлов для записи формируемой матрицы жесткости
* в ряде форматов
* @param stiffnessMatrixFileFormats - список форматов файлов, в которые будет записана
* формируемая матрица жесткости
* @param stiffnessEpsilon - число, использующееся в численном нахождении производной
* при формировании матрицы жесткости
*/
void StressStrainFortranIterativeSolver::CreateStiffnessMatrix
	(
		const IsSealedFlagsList& isSealedFlagsList,
		const string& stiffnessMatrixFileName,
		StiffnessMatrixFileFormats stiffnessMatrixFileFormats,
		double stiffnessEpsilon
	)
{
	std::cout << "Create Stiffness Matrix" << std::endl;

	StiffnessMatrixType stiffnessMatrix = CreateStiffnessMatrix(isSealedFlagsList, stiffnessEpsilon);

	std::cout << "Dump Stiffness Matrix" << std::endl;
	stiffnessMatrix.WriteToFile(stiffnessMatrixFileName, stiffnessMatrixFileFormats);

	//StiffnessMatrixPackedInCSCFormat stiffnessMatrixCSC;

	//stiffnessMatrixCSC.ReadFromBinaryFile("StiffnessMatrixPackedCSCF.bin");
	//stiffnessMatrixCSC.WriteToTextFile("StiffnessMatrixPackedCSCF2.txt");
}

/**
* Сформировать матрицу жесткости
* @param isSealedFlagsList - список наборов флагов закрепленных степеней свободы
* элементов сеточного представления
* @param stiffnessEpsilon - число, использующееся в численном нахождении производной
* при формировании матрицы жесткости
* @return сформированная матрица жесткости
*/
StressStrainFortranIterativeSolver::StiffnessMatrixType StressStrainFortranIterativeSolver::CreateStiffnessMatrix
	(
		const IsSealedFlagsList& isSealedFlagsList,
		double stiffnessEpsilon
	)
{
	StiffnessMatrixType stiffnessMatrix(_nElements);

	FindStiffnessMatrix(stiffnessMatrix, isSealedFlagsList, stiffnessEpsilon);

	return stiffnessMatrix;
}

/**
* Заполнить матрицу жесткости
* @param stiffnessMatrix - матрица жесткости, которая заполняется значениями
* @param isSealedFlagsList - список наборов флагов закрепленных степеней свободы
* элементов сеточного представления
* @param stiffnessEpsilon - число, использующееся в численном нахождении производной
* при формировании матрицы жесткости
*/
void StressStrainFortranIterativeSolver::FindStiffnessMatrix
	(
		StiffnessMatrixType& stiffnessMatrix,
		const IsSealedFlagsList& isSealedFlagsList,
		double stiffnessEpsilon
	)
{
	for (int nodeId = 0; nodeId < _nElements; ++nodeId)
	{
		FindStiffnessMatrixForAdjElements(stiffnessMatrix, nodeId, isSealedFlagsList, stiffnessEpsilon);
	}
}

/**
* Заполнить матрицу жесткости значениями всех соседних элементов к элементу с данным индексом
* @param stiffnessMatrix - матрица жесткости, которая заполняется значениями
* @param nodeId - индекс данного элемента
* @param isSealedFlagsList - список наборов флагов закрепленных степеней свободы
* элементов сеточного представления
* @param stiffnessEpsilon - число, использующееся в численном нахождении производной
* при формировании матрицы жесткости
*/
void StressStrainFortranIterativeSolver::FindStiffnessMatrixForAdjElements
	(
		StiffnessMatrixType& stiffnessMatrix,
		int nodeId,
		const IsSealedFlagsList& isSealedFlagsList,
		double stiffnessEpsilon
	)
{
	FindStiffnessMatrixForSealedElement(stiffnessMatrix, nodeId, isSealedFlagsList);
	FindStiffnessMatrixForElement(stiffnessMatrix, nodeId, nodeId, isSealedFlagsList, stiffnessEpsilon);
	for (int adjIndex = 0; adjIndex < FreedomsCount; ++adjIndex)
	{
		int adjNodeId = MLINK[nodeId * FreedomsCount + adjIndex] - 1;

		if (adjNodeId >= 0)
		{
			FindStiffnessMatrixForElement
				(
					stiffnessMatrix,
					nodeId,
					adjNodeId,
					isSealedFlagsList,
					stiffnessEpsilon
				);
		}
	}
}

/**
* Заполнить матрицу жесткости значениями соседнего элемента к элементу с данным индексом
* @param stiffnessMatrix - матрица жесткости, которая заполняется значениями
* @param nodeId - индекс данного элемента
* @param adjNodeId - индекс соседнего элемента
* @param isSealedFlagsList - список наборов флагов закрепленных степеней свободы
* элементов сеточного представления
* @param stiffnessEpsilon - число, использующееся в численном нахождении производной
* при формировании матрицы жесткости
*/
void StressStrainFortranIterativeSolver::FindStiffnessMatrixForElement
	(
		StiffnessMatrixType& stiffnessMatrix,
		int nodeId,
		int adjNodeId,
		const IsSealedFlagsList& isSealedFlagsList,
		double stiffnessEpsilon
	)
{
	const IsSealedFlags& isSealedFlags = isSealedFlagsList[nodeId];

	for (int dof = 0; dof < FreedomsCount; ++dof)
	{
		if (!isSealedFlags[dof])
		{
			FindStiffnessMatrixForDof
				(
					stiffnessMatrix,
					nodeId,
					dof,
					adjNodeId,
					isSealedFlagsList,
					stiffnessEpsilon
				);
		}
	}
}

/**
* Заполнить матрицу жесткости значениями элемента, некоторые степени свободы которого закреплены
* @param stiffnessMatrix - матрица жесткости, которая заполняется значениями
* @param nodeId - индекс элемента, некоторые степени свободы которого закреплены
* @param isSealedFlagsList - список наборов флагов закрепленных степеней свободы
* элементов сеточного представления
*/
void StressStrainFortranIterativeSolver::FindStiffnessMatrixForSealedElement
	(
		StiffnessMatrixType& stiffnessMatrix,
		int nodeId,
		const IsSealedFlagsList& isSealedFlagsList
	)
{
	const IsSealedFlags& isSealedFlags = isSealedFlagsList[nodeId];

	for (int dof = 0; dof < FreedomsCount; ++dof)
	{
		if (isSealedFlags[dof])
		{
			stiffnessMatrix.FillElement(nodeId, dof, nodeId, dof, 1.0);
		}
	}
}

/**
* Заполнить матрицу жесткости значениями соседнего элемента к элементу с данным индексом
* и данной степенью свободы
* @param stiffnessMatrix - матрица жесткости, которая заполняется значениями
* @param nodeId - индекс данного элемента
* @param dof - индекс данной степени свободы
* @param adjNodeId - индекс соседнего элемента
* @param isSealedFlagsList - список наборов флагов закрепленных степеней свободы
* элементов сеточного представления
* @param stiffnessEpsilon - число, использующееся в численном нахождении производной
* при формировании матрицы жесткости
*/
void StressStrainFortranIterativeSolver::FindStiffnessMatrixForDof
	(
		StiffnessMatrixType& stiffnessMatrix,
		int nodeId,
		int dof,
		int adjNodeId,
		const IsSealedFlagsList& isSealedFlagsList,
		double stiffnessEpsilon
	)
{
	const double* forces = &_dataInternal[_nElements * 2 * FreedomsCount];
	const double* adjForces = &forces[adjNodeId * FreedomsCount];
	const IsSealedFlags& isSealedFlags = isSealedFlagsList[adjNodeId];

	UpdateForcesForNode(nodeId, dof, adjNodeId, stiffnessEpsilon);
	for (int adjDof = 0; adjDof < 6; ++adjDof)
	{
		if (!isSealedFlags[adjDof])
		{
			const double adjForce = adjForces[adjDof];

			// TODO: вставить использование RealComparing
			if (fabs(adjForce) > 1e-10)
			{
				const double stiffness = adjForce / stiffnessEpsilon;

				stiffnessMatrix.FillElement(adjNodeId, adjDof, nodeId, dof, stiffness);
			}
		}
	}
}

/**
* Сформировать массив, содержащий элементы матрицы поворота, на основе массива углов Эйлера
* @param rotationMatrix - формируемый массив, содержащий элементы матрицы поворота
* @param eulerAngles - массив углов Эйлера
*/
static
void FormRotationMatrix
	(
		double* rotationMatrix,
		const double (&eulerAngles)[3]
	)
{
	const double xc = cos(eulerAngles[0]);
	const double yc = cos(eulerAngles[1]);
	const double zc = cos(eulerAngles[2]);
	const double xs = sin(eulerAngles[0]);
	const double ys = sin(eulerAngles[1]);
	const double zs = sin(eulerAngles[2]);

	rotationMatrix[0] = yc * zc;
	rotationMatrix[1] = -yc * zs;
	rotationMatrix[2] = ys;
	rotationMatrix[3] = xs * ys * zc + xc * zs;
	rotationMatrix[4] = -xs * ys * zs + xc * zc;
	rotationMatrix[5] = -xs * yc;
	rotationMatrix[6] = -xc * ys * zc + xs * zs;
	rotationMatrix[7] = xc * ys * zs + xs * zc;
	rotationMatrix[8] = xc * yc;
}

/**
* Сформировать массив, содержащий элементы матрицы поворота на некоторый угол
* @param rotationMatrix - формируемый массив, содержащий элементы матрицы поворота
* @param angle - угол, на который производится поворот
* @param dof - индекс координатной оси, вокруг которой производится поворот
*/
static
void FormRotationMatrix
	(
		double* rotationMatrix,
		double angle,
		int dof
	)
{
	double eulerAngles[3]{};

	eulerAngles[dof] = angle;
	FormRotationMatrix(rotationMatrix, eulerAngles);
}

/**
* Обновить значения массива ускорений для элемента с данным индексом, изменив значение приложенной
* по направлению changedDof силы к элементу с индексом changedNodeId
* @param changedNodeId - индекс элемента, к которому прикладывается сила
* @param changedDof - индекс степени свободы, по которой прикладывается сила
* @param nodeId - индекс данного элемента
* @param stiffnessEpsilon - число, использующееся в численном нахождении производной
* при формировании матрицы жесткости
*/
void StressStrainFortranIterativeSolver::UpdateForcesForNode
	(
		int changedNodeId,
		int changedDof,
		int nodeId,
		double stiffnessEpsilon
	)
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	double* rotationMatrix = _dataRotationMtx + changedNodeId * 9;
	double rotationMatrixBackup[9];

	memcpy(rotationMatrixBackup, rotationMatrix, sizeof(double) * 9);
	if (changedDof > 2)
	{
		FormRotationMatrix(rotationMatrix, stiffnessEpsilon, changedDof - 3);
	}

	double* GR1 = _dataInternal;
	double& changedNodeDofValue = _dataInternal[changedNodeId * FreedomsCount + changedDof];
	const double nodeDofValueBackup = changedNodeDofValue;

	changedNodeDofValue += stiffnessEpsilon;

	const double elasticModulus = _elasticModulusScaled;
	const double G = _elasticModulusScaled / 2;
	const double dampingFactorForward = 0.01 * 2 * _dampingFactor *
		sqrt(_elasticModulusScaled * _cellMass * _gridStep);
	const double dampingFactorRotatory = 10 * _dampingFactor *
		sqrt(2 * _elasticModulusScaled * _cellMass * _gridStep / 3) * _gridStep * _gridStep;
	double sigma[6];

	int NADU;
	int i, j, k, l, N2;
	double DX, DY, DZ, SILA[6];

	double X1[3], X2[3], SL[6], VL[6], A[36], C[36];

	int nLinks = 0;
	j = nodeId;

	NADU = _nElements * 2 * 6 + j * 6;

	for (i = 0; i < 6; i++) // можно вынести
	{
		GR1[NADU + i] = 0.0;
	}

	for (i = 0; i < 6; i++)
	{
		N2 = MLINK[j * 6 + i] - 1;

		if (N2 >= 0)
		{
			nLinks++;
			DX = (_nodes[j * 3] - _nodes[N2 * 3]) * 0.5;
			DY = (_nodes[j * 3 + 1] - _nodes[N2 * 3 + 1]) * 0.5;
			DZ = (_nodes[j * 3 + 2] - _nodes[N2 * 3 + 2]) * 0.5;

			X1[0] = -DX;
			X1[1] = -DY;
			X1[2] = -DZ;

			X2[0] = DX;
			X2[1] = DY;
			X2[2] = DZ;

			linksh3(X1, X2, SL, VL, A, C, j, N2, _nElements);
			//linksh(X1, X2, SL, VL, A, C, j, N2, _nElements);

			for (k = 0; k < 6; k++)
			{
				sigma[k] = SL[k];
			}

			for (k = 0; k < 6; k++)
			{
				SILA[k] = 0;
				if (k < 3)
				{
					SILA[k] = -VL[k] * dampingFactorForward;
				}
				else
				{
					SILA[k] = -VL[k] * dampingFactorRotatory;
				}
				if (k < 3)
				{
					SILA[k] -= sigma[k] * elasticModulus * _gridStep; // сила
				}
				else
				{
					//SILA[k] -= sigma[k] * E * _gridStep;
					//SILA[k] -= sigma[k] * elasticModulus * _gridStep * _gridStep * _gridStep; // момент
					SILA[k] -= sigma[k] * elasticModulus * _gridStep * _gridStep;
					//SILA[k] -= SL[k] * E * _gridStep; // ???
				}
			}
			for (k = 0; k < 6; k++)
			{
				for (l = 0; l < 6; l++)
				{
					double a = A[6 * k + l];

					if ((k < 3) && (l > 2))
					{
						a = -a;
					}
					//GR1[NADU + l] += SILA[k] * A[6 * k + l];
					GR1[NADU + l] += SILA[k] * a;
				}
			}
			int dbg = 1;
		}
	}
	//std::cout << "NODE=" << nodeId << std::endl;
	//std::cout << "CHANGE NODE=" << changedNodeId << std::endl;
	//std::cout << "CHANGE NODE DOF=" << changedDof << std::endl;
	//for (l = 0; l < 6; l++)
	//{
	//	std::cout << GR1[NADU + l] << ' ';
	//}
	//std::cout << std::endl;
	memcpy(rotationMatrix, rotationMatrixBackup, sizeof(double) * 9);
	changedNodeDofValue = nodeDofValueBackup;
}

/**
* Сформировать вектор правых частей СЛАУ на основе матрицы жесткости
* @param isSealedFlagsList - список наборов флагов закрепленных степеней свободы
* элементов сеточного представления
* @param stiffnessRHSFileName - общее наименование файлов для записи формируемого вектора
* правых частей СЛАУ на основе матрицы жесткости
* @param writeStiffnessRHSVectorToTextFileFlag - флаг, следует ли записывать формируемый вектор
* правых частей СЛАУ на основе матрицы жесткости в текстовый файл
*/
void StressStrainFortranIterativeSolver::CreateStiffnessRHSVector
	(
		const IsSealedFlagsList& isSealedFlagsList,
		const string& stiffnessRHSFileName,
		bool writeStiffnessRHSVectorToTextFileFlag
	)	const
{
	StiffnessRHSVector stiffnessRHSVector = CreateStiffnessRHSVector(isSealedFlagsList);

	WriteStiffnessRHSVectorToBinaryFile(stiffnessRHSFileName + ".bin", stiffnessRHSVector);
	if (writeStiffnessRHSVectorToTextFileFlag)
	{
		WriteStiffnessRHSVectorToTextFile(stiffnessRHSFileName + ".txt", stiffnessRHSVector);
	}
}

/**
* Сформировать вектор правых частей СЛАУ на основе матрицы жесткости
* @param isSealedFlagsList - список наборов флагов закрепленных степеней свободы
* элементов сеточного представления
* @return сформированный вектор правых частей СЛАУ на основе матрицы жесткости, заполняемый значениями
*/
StiffnessRHSVector StressStrainFortranIterativeSolver::CreateStiffnessRHSVector
	(
		const IsSealedFlagsList& isSealedFlagsList
	)	const
{
	StiffnessRHSVector stiffnessRHSVector(_nElements * FreedomsCount);

	FindStiffnessRHSVector(stiffnessRHSVector);
	ApplySealedElementsToStiffnessRHSVector(stiffnessRHSVector, isSealedFlagsList);
	InvertStiffnessRHSVector(stiffnessRHSVector);

	return stiffnessRHSVector;
}

/**
* Заполнить значениями вектор правых частей СЛАУ на основе матрицы жесткости
* @param stiffnessRHSVector - вектор правых частей СЛАУ на основе матрицы жесткости, заполняемый значениями
*/
void StressStrainFortranIterativeSolver::FindStiffnessRHSVector
	(
		StiffnessRHSVector& stiffnessRHSVector
	)	const
{
	for (const auto& boundaryParams : _boundaryParamsSet)
	{
		if (boundaryParams.GetKind() == 3)
		{
			boundaryParams.ApplyForceBoundary(stiffnessRHSVector.data());
		}
	}
}