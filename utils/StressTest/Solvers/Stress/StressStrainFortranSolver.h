#ifndef StressStrainFortranSolverH

#define StressStrainFortranSolverH


#include "StressStrainSolver.h"
#include "BoundaryParams.h"

#include <cmath>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include "FTimer.h"


using std::string;
using std::vector;


#define MAX(x, y) ((x) > (y) ? x : y)


double DMOD(const double d1, const double d2);

struct Fue
{
	double* R;
	double* R1;
	double* R2;
	double* R3;
	double* UE;
	double* U1;

	double* R1Z;
	double* RZ;
	double* A1;
	double* FLAGRY;
	int IFLAGRY;

	Fue	(const int nNodes)
	{
		R = new double [nNodes*9];
		R1 = new double [nNodes*9];
		R2 = new double [nNodes*9];
		R3 = new double [nNodes*9];
		UE = new double [nNodes*9];
		U1 = new double [nNodes*9];
		R1Z = new double [nNodes*9];
		RZ =  new double [nNodes*9];
		A1 =  new double [nNodes*9];
		FLAGRY =  new double [nNodes*3];
	}

	~Fue ()
	{
		delete [] R;
		delete [] R1;
		delete [] R2;
		delete [] R3;
		delete [] UE;
		delete [] U1;
		delete [] R1Z;
		delete [] RZ;
		delete [] A1;
		delete [] FLAGRY;
	}

	void UpdateR(const int id, const double* om, const double H)
	{
		double COS2=cos(UE[id+1]);
		double TG2=tan(UE[id+1]);
		double SIN3=sin(UE[id+2]);
		double COS3=cos(UE[id+2]);
		double OM1=om[0];
		double OM2=om[1];
		double OM3=om[2];

		/*      
		C     Контроль точности
		IF(COS2.EQ.0.0)THEN
		PRINT *,'COS2 РАВЕН 0 , УРАВНЕНИЕ 2',JJ
		STOP
		ENDIF*/


		double CU1=(-SIN3*OM1-COS3*OM2)/COS2;
		double CU2=0.0;
		double CU3=(OM2*COS3+OM1*SIN3)*TG2;
		double CU0=400.0/(H*H);

		/*
		if (DABS(CU1)>CU0)
		#print *,' Мал шаг для интегрирования 1-го уравнения Эйлера ',jj
		#,'4 COS2 OM1 OM2',COS2,OM1,OM2

		IF(DABS(CU3)>CU0)
		#print *,' Мал шаг для интегрирования 3-го уравнения Эйлера ',jj
		#,'TG2 OM1 OM2',TG2,OM1,OM2
		*/

		for (int j=id;j<id+3;j++)
			RZ[j] = R[j];	
		R[id]=(OM1*COS3-OM2*SIN3)/COS2*H;
		R[id+1]=(OM1*SIN3+OM2*COS3)*H;
		R[id+2]=((OM2*SIN3-OM1*COS3)*TG2+OM3)*H;
	}

	void UpdateR2(const int id, const int mets)
	{

		for(int j = id; j < id+3; j++)
		{
			switch(mets)
			{
			case 1:
				R1Z[j]=R1[j];
				R1[j]=R[j];
				break;
			case 2:
				R2[j]=R[j];
				break;
			case 3:
				R3[j]=R[j];
				break;
			}
		}
	}

	void UpdateMtx(const int id, double* a)
	{
		double xc=cos(UE[id]);
		double yc=cos(UE[id+1]);
		double zc=cos(UE[id+2]);
		double xs=sin(UE[id]);
		double ys=sin(UE[id+1]);
		double zs=sin(UE[id+2]);

		a[0]=yc*zc;
		a[1]=-yc*zs;
		a[2]=ys;
		a[3]=xs*ys*zc+xc*zs;
		a[4]=-xs*ys*zs+xc*zc;
		a[5]=-xs*yc;
		a[6]=-xc*ys*zc+xs*zs;
		a[7]=xc*ys*zs+xs*zc;
		a[8]=xc*yc;	
	}

	void Update(const int id, const int mets = 0)
	{
		for (int j=id;j<id+3;j++)
		{
			switch(mets)
			{
				case 0:
					R[j]=0.0; 
					R1[j]=0.0;
					R2[j]=0.0;
					R3[j]=0.0;
					UE[j]=0.0;
					U1[j]=0.0;
				break;
				case 1:
					UE[j]=U1[j]+(R1[j]+(R2[j]+R3[j])*2+R[j])/6.0;
					U1[j]=UE[j];
				break;
				case 2:
				case 3:
					UE[j]=U1[j]+R1[j]*0.5;
				break;
				case 4:
					UE[j]=U1[j]+R3[j];
				break;
			}
		}
	}

};

struct Copym
{
	bool* isFirstCopy;
	//double T;
    //COMMON/UGR/ KGR1,KGR2,T
	double* TM;
	double* AZ;
	double T0;
      
	Copym(const int nNodes)
	{
		isFirstCopy = new bool [nNodes];
		for(int i = 0 ; i < nNodes; i++)
			isFirstCopy[i] = true;
		TM = new double [nNodes];
		AZ = new double [nNodes*9];
		T0 = 0;
	}
	
	~Copym()
	{
		delete [] isFirstCopy;
		delete [] TM;
		delete [] AZ;
	}

	void Copy(double* A, double* A1, const int N, const double T)
	{
		if(isFirstCopy[N] || T==T0)
		{
			memcpy(A1,A,9*sizeof(double));
			isFirstCopy[N] = false;
			TM[N]=T;
		}
		else 
		{
			const int KA=N*9;
			if(TM[N] != T)
			{
				memcpy(A1,A,9*sizeof(double));
				TM[N]=T;
				memcpy(AZ+KA,A,9*sizeof(double));
			}
			else
			{
				memcpy(A1,AZ+KA,9*sizeof(double));
			}
		}
	}
};


/** Унаследованный код на ФОРТРАНе
*/
class StressStrainFortranSolver
	:
		public StressStrainSolver
{
public:
	
	// создает объект с заданными параметрами
	StressStrainFortranSolver
		(
			double* params,
			int* links,
			int nLinks, 
			double *nodes,
			int nNodes,
			double gridStep, 
			double timeStep,
			int numThreads
		);

	virtual
	~StressStrainFortranSolver();


#pragma region overriden

	virtual
	void AddLinks
		(
			const int* links
		);

	virtual
	void AddBoundary
		(
			int* boundaryNodesIndices, 
			int numberOfBoundaryNodes,
			int bcKind,
			double* bcParams
		);
	
	void AddPartialBoundary
		(
			int* boundaryNodesIndicesInPart, 
			int numberOfNodesInPart,
			int numberOfNodes,
			int bcKind,
			double* bcParams
		);

	virtual
	void ChangeBoundaryParam
		(
			const int bcNumber,
			const int bcParamNumber,
			const double bcParamValue
		);

	virtual
	double GetBoundaryParam
		(
			const int bcNumber,
			const int bcParamNumber
		);

	virtual
	void UpdateBuffer
		(
			double scaleFactor
		);

#pragma endregion


protected:

	bool _isFirstSolution;
	vector<BoundaryParams> _boundaryParamsSet;

	//double* GRA; // массив матриц поворота
// кэш-массивы, используемые тольк ов solve
	double *_initX;
	double *_initDX;
	double *_hDDX1;
	double *_hDDX2;
	double *_hDDX3;
// кэш-массивы

	double* _varX;	// смещения до неизвестных
	double* _varDX; // смещения до производных
	double* _varDDX;// смещения до 2 производных

	double *R;
	double *RZ;
	double *R1Z;
	int* MLINK;
	double* _nodes;
	
	Fue* _fue;
	Copym* _copym;

	int _nNodes;
	int _nVariables;
	int _nIteration;
	int METS;
	double T; // время
	double TZ; // вспомогательная функция для сохранения промежуточного значения времени
	double H; // шаг по времени
	bool _isFirstIteration;
	//double* GR1; // массив результатов
	double _gridStep; // шаг сетки
	double _elasticModulus; // модуль упругости
	double _elasticModulusScaled; // отмасштабированный модуль упругости
	double _dampingFactor; // коэффициент демпфирования
	double _density; // плотность материала
	double _cellMass; // масса ячейки
	double _stiffScale; // масштабирование
	double _poissonRatio; // коэффициент Пуассона
	int NumThreads;
	FTimer _testTimer;

	void intomsub();

	void urejl4s
		(
			double* a,
			double *ug,
			double *om,
			int mets,
			const int id
		);

	void multm2
		(
			double *a1,
			double *a2
		);

	void linksh
		(
			double *X1,
			double *X2,
			double *SL,
			double *VL,
			double *A,
			double *C,
			int j,
			int N2,
			int KUZ
		);

	void linksh2
	(
		double *cVec1,	//координаты точки связи элемента 1
		double *cVec2,	//координаты точки связи элемента 2
		double *SL,		// выход деформаций
		double *VL,		// выход изм. скоростей
		double *A,		// выход
		double* C,		// выход 
		int nodeId1,	// номер узла 1
		int nodeId2,	// номер узла 2
		int nNodes		// количество узлов
	);

	void linksh3
	(
		double *cVec1,	//координаты точки связи элемента 1
		double *cVec2,	//координаты точки связи элемента 2
		double *SL,		// выход деформаций
		double *VL,		// выход изм. скоростей
		double *A,		// выход
		double* C,		// выход 
		int nodeId1,	// номер узла 1
		int nodeId2,	// номер узла 2
		int nNodes		// количество узлов
	);

	void CompareAngles
		(
			double* A1,
			double* A2,
			double* SL
		);

	void FindStressStrainMatrix();

	/**
	* Получить векторный параметр
	* @param coordinatesData - массив для записи векторного параметра
	* @param scaleFactor - масштабный коэффициент
	*/
	void GetVectorParameter
		(
			float* coordinatesData,
			double scaleFactor = 1.0
		)	const;

};
#endif