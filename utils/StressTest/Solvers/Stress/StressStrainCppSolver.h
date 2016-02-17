#pragma once

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

namespace Stress
{

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

	int vecStride, vecStride2, matStride;
	Fue	(const int nNodes, int stride):
		vecStride(stride),
		vecStride2(stride * 2),
		matStride(stride * 3)
	{
		const size_t matSize = nNodes*matStride*sizeof(double);
		R =  (double*)_aligned_malloc(matSize,alignment);
		R1 = (double*)_aligned_malloc(matSize,alignment);
		R2 = (double*)_aligned_malloc(matSize,alignment);
		R3 = (double*)_aligned_malloc(matSize,alignment);
		UE = (double*)_aligned_malloc(matSize,alignment);
		U1 = (double*)_aligned_malloc(matSize,alignment);
		R1Z = (double*)_aligned_malloc(matSize,alignment);
		RZ =  (double*)_aligned_malloc(matSize,alignment);
		A1 =  (double*)_aligned_malloc(matSize,alignment);
		FLAGRY = (double*)_aligned_malloc(nNodes*vecStride*sizeof(double),alignment);
		memset(R,0,matSize);
		memset(R1,0,matSize);
		memset(R2,0,matSize);
		memset(R3,0,matSize);
		memset(UE,0,matSize);
		memset(U1,0,matSize);
		memset(R1Z,0,matSize);
		memset(RZ,0,matSize);
		memset(A1,0,matSize);
		memset(FLAGRY,0,nNodes*vecStride*sizeof(double));
	}

	~Fue ()
	{
		_aligned_free(R);
		_aligned_free(R1);
		_aligned_free(R2);
		_aligned_free(R3);
		_aligned_free(UE);
		_aligned_free(U1);
		_aligned_free(R1Z);
		_aligned_free(RZ);
		_aligned_free(A1);
		_aligned_free(FLAGRY);
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
		C     �������� ��������
		IF(COS2.EQ.0.0)THEN
		PRINT *,'COS2 ����� 0 , ��������� 2',JJ
		STOP
		ENDIF*/


		double CU1=(-SIN3*OM1-COS3*OM2)/COS2;
		double CU2=0.0;
		double CU3=(OM2*COS3+OM1*SIN3)*TG2;
		double CU0=400.0/(H*H);

		/*
		if (DABS(CU1)>CU0)
		#print *,' ��� ��� ��� �������������� 1-�� ��������� ������ ',jj
		#,'4 COS2 OM1 OM2',COS2,OM1,OM2

		IF(DABS(CU3)>CU0)
		#print *,' ��� ��� ��� �������������� 3-�� ��������� ������ ',jj
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


		double* firstRow = a;
		double* secondRow = a + vecStride;
		double* thirdRow = a + vecStride*2;
		
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

	void Update(const int id, const int mets = 0)
	{
		if (id == 280)
			int debug = 0;
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
      
	int vecStride, vecStride2, matStride;
	Copym(const int nNodes, int stride) :
		vecStride(stride),
		vecStride2(stride * 2),
		matStride(stride * 3)
	{
		isFirstCopy = new bool [nNodes];
		for(int i = 0 ; i < nNodes; i++)
			isFirstCopy[i] = true;
		TM = new double [nNodes];
		AZ = new double [nNodes*matStride];
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
			memcpy(A1,A,matStride*sizeof(double));
			isFirstCopy[N] = false;
			TM[N]=T;
		}
		else 
		{
			const int KA=N*matStride;
			if(TM[N] != T)
			{
				memcpy(A1,A,matStride*sizeof(double));
				TM[N]=T;
				memcpy(AZ+KA,A,matStride*sizeof(double));
			}
			else
			{
				memcpy(A1,AZ+KA,matStride*sizeof(double));
			}
		}
	}
};


/** �������������� ��� �� ��������
*/
class StressStrainCppSolver
	:
		public StressStrainSolver
{
public:
	
	// ������� ������ � ��������� �����������
	StressStrainCppSolver
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
		);

	virtual
	~StressStrainCppSolver();


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
			double scale
		);

#pragma endregion


//protected:

	bool _isFirstSolution;
	vector<BoundaryParams> _boundaryParamsSet;

	
// ���-�������, ������������ ����� �� solve
	double *_initX;
	double *_initDX;
	double *_hDDX1;
	double *_hDDX2;
	double *_hDDX3;
// ���-�������

	double* _varX;		// �������� �� �����������
	double* _varDX;		// �������� �� �����������
	double* _varDDX;	// �������� �� 2 �����������

	double *R;
	double *RZ;
	double *R1Z;
	int* _linkedElements;
	double* _elements;
	//double* _rotationMatrices; // ������ ������ ��������

	Fue* _fue;			// ��������� ��� ��������� ������
	Copym* _copym;		// ��������� ��� ����������� ???

	int _nVariables;	// ���������� �����������
	int _nIteration;	// ����� ��������
	int _stageRK;		// ��� ����� ����� 4
	double _time;		// �����
	double _timeTmp;	// ��������������� ������� ��� ���������� �������������� �������� �������
	double _timeStep;	// ��� �� �������
	double _timeStep2;	// ��� �� ������� /2
	double _timeStep4;	// ��� �� ������� /4
	bool _isFirstIteration;
	
	//double* GR1; // ������ �����������
	double _gridStep; // ��� �����
	double _gridStep2; // ��� ����� � ��������
	double _gridStep3; // ��� ����� � ����
	double _elasticModulus; // ������ ���������
	double _shearModulusScaled; // ������ ������
	double _elasticModulusScaled; // ������������������ ������ ���������
	double _dampingFactor; // ����������� �������������
	double _dampingFactorAngular; // ����������� ����������� �������� �������������
	double _dampingFactorLinear; // ����������� ����������� ��������� �������������
	double _density; // ��������� ���������
	double _cellMass; // ����� ������
	double _stiffScale; // ���������������
	double _poissonRatio; // ����������� ��������
	double _lameMatrix[6][6]; // ������� �������� �� ����������� � �����������

	int _numThreads;
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

	void MatrixMul
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
		double *cVec1,	//���������� ����� ����� �������� 1
		double *cVec2,	//���������� ����� ����� �������� 2
		double *SL,		// ����� ����������
		double *VL,		// ����� ���. ���������
		double *A,		// �����
		double* C,		// ����� 
		int nodeId1,	// ����� ���� 1
		int nodeId2,	// ����� ���� 2
		int nNodes		// ���������� �����
	);

	void linksh3
	(
		double *cVec1,	//���������� ����� ����� �������� 1
		double *cVec2,	//���������� ����� ����� �������� 2
		double *SL,		// ����� ����������
		double *VL,		// ����� ���. ���������
		double *A,		// �����
		double* C,		// ����� 
		int nodeId1,	// ����� ���� 1
		int nodeId2,	// ����� ���� 2
		int nNodes		// ���������� �����
	);

	void linksh4
	(
		size_t side,
		double *SL,		// ����� ����������
		double *VL,		// ����� ���. ���������
		double& rx,		// �����
		double& ry,		// ����� 
		double& rz,		// ����� 
		int nodeId1,	// ����� ���� 1
		int nodeId2,	// ����� ���� 2
		int nNodes		// ���������� �����
	);

	void linksh4AVX
	(
		size_t side,
		double *SL,		// ����� ����������
		double *VL,		// ����� ���. ���������
		double& rx,		// �����
		double& ry,		// ����� 
		double& rz,		// ����� 
		int nodeId1,	// ����� ���� 1
		int nodeId2,	// ����� ���� 2
		int nNodes		// ���������� �����
	);

	void FindStressStrainMatrix();
//protected:

};
};