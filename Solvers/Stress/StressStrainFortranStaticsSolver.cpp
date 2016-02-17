#include "StressStrainFortranStaticsSolver.h"


//#define NO_LINKSH
//#define NO_INTOMSUB

#ifndef M_PI

#define M_PI 3.1415926535897932384626433832795

#endif

#define MAX(x, y) ((x) > (y) ? x : y)
#define IS_ZERO(x) (fabs(x) < 1e-10)


bool SolveSLES
	(
		int* ia,
		int* ja,
		double* b,
		double* a,
		double* x,
		int n
	);

StressStrainFortranStaticsSolver::StressStrainFortranStaticsSolver
	(
		double* params, 
		int* links, 
		int nLinks, 
		double *nodes, 
		int nNodes, 
		double gridStep, 
		double timeStep,
		int numThreads
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
				3
			)
{
	InitSolveLES();
}

// virtual
StressStrainFortranStaticsSolver::~StressStrainFortranStaticsSolver()
{

}

// virtual
void StressStrainFortranStaticsSolver::Solve
	(
		const int nIterations
	)
{
	double* GR1 = _dataInternal;
	double H2 = H * 0.5;
	double H4 = H * 0.25;
	int N = _nElements * 6;

	int iteration = 0;
	int nIterationModifier = 1;

	T = (iteration + 1) * H;
	iteration++;
	double TZ = T;
	T = T - H;
	T = T + H2;
	METS = 1;
	intomsub();
	pravsubfl();
	_mtx.Sort();
	_mtx.RemoveNearestToZero();
	this->SolveLes();
	//_mtx.Print();
	//_mtx.PrintPortrait();
	for (int i = 0; i < _nElements*6; i++)
	{
		//_dataInternal[i] += (_mtx.x[i] / (_elasticModulus / (_gridStep * _gridStep)));
		_dataInternal[i] += (_mtx.x[i] / (_elasticModulus * _gridStep));
		//std::cout << _dataInternal[i] << std::endl;
	}
	int dbg = 0;
}

/**
* Расчет первой стадии метода Рунге-Кутты
*/
// virtual
void StressStrainFortranStaticsSolver::Solve1()
{
	// ничего не надо делать
}

/**
* Расчет второй стадии метода Рунге-Кутты
*/
// virtual
void StressStrainFortranStaticsSolver::Solve2()
{
	// ничего не надо делать
}

/**
* Расчет третьей стадии метода Рунге-Кутты
*/
// virtual
void StressStrainFortranStaticsSolver::Solve3()
{
	// ничего не надо делать
}

/**
* Расчет четвертой стадии метода Рунге-Кутты
*/
// virtual
void StressStrainFortranStaticsSolver::Solve4()
{
	// ничего не надо делать
}

/**
* Расчет пятой стадии метода Рунге-Кутты
*/
// virtual
void StressStrainFortranStaticsSolver::Solve5()
{
	// ничего не надо делать
}


/** Получить смещения
* @param data - массив для записи смещений как скалярного параметра
*/
// virtual
void StressStrainFortranStaticsSolver::GetDisplacement
	(
		float* data
	)
{
	// TODO: implement
	throw "Not Implemented";
}

/** Получить напряжения по первой теории прочности
* @param data - массив для записи напряжений как скалярного параметра
*/
// virtual
void StressStrainFortranStaticsSolver::GetStressesByFirstTheoryOfStrength
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
			coordinates[j] = _mtx.x[6 * i + j];
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
					shift = _mtx.x[6 * neighbourNodeNumber + shiftIndex] - coordinates[shiftIndex];
					if (fabs(shift) > fabs(relativeShifts[shiftIndex]))
					{
						relativeShiftsSigned[shiftIndex] = shift * (2 * k - 1);
						//relativeShifts[shiftIndex] = fabs(shift);
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
		data[i] = (float)(fullRelativeShift / (_gridStep * _gridStep));
	}
	int dbg = 0;
}

/** Получить напряжения по von Mises
* @param data - массив для записи напряжений как скалярного параметра
*/
// virtual
void StressStrainFortranStaticsSolver::GetStressesByVonMises
	(
		float* data
	)
{
	// TODO: implement
	throw "Not Implemented";
}

void StressStrainFortranStaticsSolver::MakeKmtx1d
	(
		int n, 
		int* &ia,
		int* &ja,
		double* &a
	)
{
	int nElements = 2*n-1;
	a[0]=a[nElements-1]=1;
	a[1]=a[nElements-2]=-1;
	ja[0]=1;ja[1]=2;
	ja[nElements-1]=n;
	ia[0]=1;ia[1]=3;
	for(int i=2, id=3; i < n; i++, id+=2)
	{
		ia[i]=2*i+1;		
		ja[id-1]=i;
		ja[id]=i+1;
		a[id-1]	=2;
		a[id]	=-1;
	}
	ia[n-1] = 2*(n-1)+1;
	ia[n] = ia[n-1]+1;
}

//void StressStrainFortranStaticsSolver::MakeKmtx1d(int n, 
//		int* &ia, int* &ja, double* &a)
//{
//	int nElements = 3*n-2;
//	a = new double[nElements];
//	ia = new int[n+1];
//	ja = new int[nElements];
//	
//	a[0]=a[nElements-1]=1;
//	a[1]=a[nElements-2]=-1;
//	ja[0]=0;ja[1]=1;
//	ja[nElements-2]=n-2;ja[nElements-1]=n-1;
//	ia[0]=0;ia[1]=2;
//	for(int i=1, id=2; i < n-1; i++, id+=3)
//	{
//		ia[i]=3*i-1;		
//		ja[id]	=i-1;
//		ja[id+1]=i;
//		ja[id+2]=i+1;	
//		a[id]	=-1;
//		a[id+1]	=2;
//		a[id+2]	=-1;
//	}
//	ia[n-1] = 3*(n-1)-1;
//	ia[n] = ia[n-1]+2;
//}

void StressStrainFortranStaticsSolver::InitSolveLES()
{
	if(_isFirstSolution)
	{
		_mtx.Allocate3d(_nElements);

		_isFirstSolution = false;
	}
}
	
bool StressStrainFortranStaticsSolver::SolveLes()
{
	std::cout << "Init solution\n";
	ApplyBoundary();
	bool res = false;
#ifdef USE_MKL
	if(!_pardisoSolver.Init(_mtx._size, _mtx.ia, _mtx.ja, _mtx.a, -2, 0, 0))
	{
		throw "PARDISO INIT ERROR!";
	}
	std::cout << "Solution\n";
	res = _pardisoSolver.Solve(_mtx.rhs, _mtx.x);
	std::cout << "EOF solution\n";
#else
	throw "MKL IS NOT USED (#define USE_MKL)";
#endif
	return res;
}

void StressStrainFortranStaticsSolver::ApplyBoundary()
{
	vector<BoundaryParams>::iterator it = _boundaryParamsSet.begin();
	
	while (it != _boundaryParamsSet.end())
	{
		switch (it->GetKind())
		{
		case 3:
			ApplyForceBoundary(*it);
			break;

		case 4:
			ApplySealedBoundary(*it);
			break;

		default:
			break;
		}
		it++;
	}
}

void StressStrainFortranStaticsSolver::ApplyForceBoundary
	(
		const BoundaryParams& boundaryParams
	)
{
	const int* nodes;
	int nodesCount;

	boundaryParams.GetNodes(nodes, nodesCount);
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < nodesCount; j++)
		{
			_mtx.ApplyRhsTo3dNode(nodes[j] - 1, i, boundaryParams.GetParam(i) / nodesCount); // равномерное распределение
		}
	}	
}

void StressStrainFortranStaticsSolver::ApplySealedBoundary
	(
		const BoundaryParams& boundaryParams
	)
{
	const int* nodes;
	int nodesCount;

	boundaryParams.GetNodes(nodes, nodesCount);
	for (int i = 0; i < 6; i++)
	{
		if (boundaryParams.GetParam(i) < 0)
		{
			for (int j = 0; j < nodesCount; j++)
			{
				_mtx.Fix3dNode(nodes[j] - 1, i);
			}
		}
	}
}

void PrintMatrix
	(
		int rowsCount,
		int columnsCount,
		const double* matrix,
		const string& matrixName
	)
{
	std::cout << std::setw(10) << matrixName;
	for (int j = 0; j < columnsCount; j++)
	{
		std::cout << std::setw(10) << j;
	}
	std::cout << std::endl;
	for (int i = 0; i < rowsCount; i++)
	{
		std::cout << std::setw(10) << i;
		for (int j = 0; j < columnsCount; j++)
		{
			std::cout << std::setw(10) << matrix[rowsCount * i + j];
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void StressStrainFortranStaticsSolver::pravsubfl() 
{   
	double* GR1 = _dataInternal;
    static int it = 0;
	if(_isFirstIteration)
	{
		_isFirstIteration = false;
		METS=0;
		intomsub();
		METS=1;
	}

	double E = _elasticModulusScaled * _gridStep;
	double DM = 2 * _dampingFactor * sqrt(_elasticModulusScaled * _cellMass);
	//double PUAS=0;
	//double DH=0.1;
	
	int NADU;
	int i,j,k,N2;
	double DX,DY,DZ,SILA;

	double X1[3],X2[3],SL[6],VL[6],A[36],C[36];
	//memset(A,0,sizeof(double)*36);
	//memset(C,0,sizeof(double)*36);
	memset(_mtx.ia,0,sizeof(int)*(_mtx._size+1));

	#pragma omp parallel for private(NADU,j,i,DX,DY,DZ,X1,X2,SL,VL,A,C,N2,SILA) num_threads(NumThreads)
	for (j=0;j<_nElements;j++)
	{    
		_stress[j] = 0;
		// NAD=j*6;
		NADU=_nElements*2*6+j*6;         

		for (i=0;i<6;i++) // можно вынести
		{
			GR1[NADU+i]=0.0;
		}

		for (i=0;i<6;i++)
		{
			N2=MLINK[j*6+i]-1;

			if (N2>=0)
			{
				DX=(_nodes[j*3]-_nodes[N2*3])*0.5;
				DY=(_nodes[j*3+1]-_nodes[N2*3+1])*0.5;
				DZ=(_nodes[j*3+2]-_nodes[N2*3+2])*0.5;

				X1[0]=-DX;
				X1[1]=-DY;
				X1[2]=-DZ;

				X2[0]=DX;
				X2[1]=DY;
				X2[2]=DZ;
				
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
				linksh(X1,X2,SL,VL,A,C,j,N2,_nElements);
				
				Genmatl(j,N2,A,C);
				//std::cout << "Node: " << j << " EOF Genmatl\n";
#endif
				//double SL1[6], VL1[6];

				for (k=0;k<6;k++)
				{
					//if(k == 0)
					//{
					//	SL1[k]=(fabs(GR1[j*6+k]-GR1[(N2-1)*6+k])-_gridStep);
					//	if(j < (N2-1))
					//		SL1[k] = -SL1[k];
					//	VL1[k]=(GR1[j*6+_nElements*6+k]-GR1[(N2-1)*6+_nElements*6+k]);
					//}
					//else
					//{
					//	SL1[k]=VL1[k]=0;
					//}
					//SILA=-(SL1[k]*E+VL1[k]*DM);
#ifdef NOLINKSH
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
#endif
					SILA=-(SL[k]*E+VL[k]*DM);					
					
					for (int l=0;l<6;l++)
					{
#ifdef NOLINKSH
						if(l==0)					
						{
							GR1[NADU+l]+=SILA;
						}
#else
							GR1[NADU+l]+=SILA*A[6*l+k];
#endif
					}
					double stress = SL[k]/_gridStep*_elasticModulusScaled;
					if(stress > _stress[j])
						_stress[j] = stress;
				}
				int dbg=1;
			}
		}
		
	}
	PrintMatrix(6, 6, A, "A");
	PrintMatrix(6, 6, C, "C");
	
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
}

void AddElement
	(
		int* ia,
		int* ja,
		double* a,
		int nEq,
		int nDof,
		double am,
		int nEqs
	)
{
	if (nEq == 6)
	{
		int test = 0;
	}

	for(int i = ia[nEq]; i < ia[nEq + 1]; i++)
	{
		if (ja[i - 1] == nDof + 1)
		{
			a[i - 1] = a[i - 1] + am;
			return;
		}
	}

	if(ia[nEq+1] == 0)
		ia[nEq+1]=ia[nEq];

	int m1=0,m2=0;

	if(nEq < nEqs-1)
	{
		if(ia[nEq+2] != 0)
		{
			m1 = ia[nEq+1]-1;
			int i = 0;

			for(int i = 2; i < (nEqs-nEq)+1; i++)
			{
				if(ia[nEq+i] == 0)
					break;
				m2 = ia[nEq+i]-1;
				ia[nEq+i]+=1;
			}
			for(int i = m2-1; i >= m1; i--)
			{
				a[i+1]=a[i];
				ja[i+1]=ja[i];
			}
		}
	}
	a[ia[nEq+1]-1]=am;
	ja[ia[nEq+1]-1]=nDof+1;
	ia[nEq+1]+=1;
}

void StressStrainFortranStaticsSolver::Genmatl(int n1, int n2, double* a, double* c)
{
	double am;
	_mtx.ia[0] = 1;
	int k1 = n1 * 6;
	int k2 = n2 * 6;

	for (int k = 0; k < 6; k++)
	{
		for (int n = 0; n < 6; n++)
		if (!IS_ZERO(a[n * 6 + k]))
		{
			for (int j = 0; j < 6; j++)
			{
				// TODO: k1 повторяется
				if (!IS_ZERO(a[n * 6 + j]) && (k1 + j >= k1 + k))
				{
					am = a[n * 6 + j] * a[n * 6 + k];
					AddElement(_mtx.ia, _mtx.ja, _mtx.a, k1 + k, k1 + j, am, _mtx._size);
				}
			}
			for (int j = 0; j < 6; j++)
			{
				if (!IS_ZERO(c[n * 6 + j]) && (k2 + j >= k1 + k))
				{
					am = c[n * 6 + j] * a[n * 6 + k];
					AddElement(_mtx.ia, _mtx.ja, _mtx.a, k1 + k, k2 + j, am, _mtx._size);
				}
			}
		}
	}
}