#pragma once
#include <iostream>

using std::cout;
using std::string;

#define SQR(x) ((x) * (x))
#define MAX(x, y) ((x) > (y) ? x : y)

#define MeasuredRun(TIMER, COMMAND) \
	_testTimer.Start(TIMER); \
	COMMAND; \
	_testTimer.Stop(TIMER);

#define ALIGNMENT 64 // KNC

enum DataType
{
	DT_Shifts=0,
	DT_Speeds=1,
	DT_Accelerations=2
};

//typedef T double;
class StressStrainSolver
{
	//size_t _stride;

protected:

	double _velocitySum[3];
	double _velocitySumSingle[3];
	double _velocitySumX[3];
	double _velocitySumY[3];
	double _velocitySumZ[3];

public:

	size_t vecStride;
	size_t vecStride2;
	size_t matStride;

	StressStrainSolver
		(
			int const nNodes,
			int stride
		);

	void SetZeroVelocitiesCache();
	void SetZeroVelocities();
	void SetZeroVelocitiesX();
	void SetZeroVelocitiesY();
	void SetZeroVelocitiesZ();

	void CheckVelocitySumm();
	bool Check(double* _velocitySum);
	double GetSquareSummOfVelocities();
	virtual
	~StressStrainSolver();
	
	void InitIco ( const string& fileName,
		bool readIco,
		bool writeIco,
		int nWriteIteration);

	bool ReadIco(const string& fileName);
	void WriteIco(const string& fileName) const;

	float* GetMemoryPointer() const;

	int GetMemorySize() const;

	/** ���������� ��������� �� ������ ���������� ������
	* @param dataType - ��� ������ (�����������, ��������, ���������)
	* @return - ��������� �� ������ ������
	*/
	double* GetDataInternal(DataType dataType) const;
	
	/** ���������� ��������� �� ������ ������ ��������
	* @return - ��������� �� ������ ������ ��������
	*/
	double* GetDataRotaionMtx() const;

	/** ���������� ������
	* @param links - ������ 6*_nElements � ���������� ������ (1 - ���� �����, 0 - ��� �����)
	*/
	virtual
	void AddLinks
		(
			const int* links
		)
		= 0;

	/** �������� ��������� ������� ����������
	* @param boundaryNodesIndices - ������� ��������� �����
	* @param numberOfNodesInPartialBoundary - ���������� ��������� ����� � �������� �������
	* @param numberOfNodesInFullBoundary - ���������� ��������� ����� � ������ �������
	* @param bcKind - ��� ��������� ������ (=1 - ����, ���������� ����������� �� ������, 
	* �������� �������� �� 6 ��������� - 3 �������������� ���������� ������� ���� 
	* �� ������������ ���� � ��� ������� ������������ ����)
	* @param bcParams - ���������� ��������� �������
	* =1: 6 ���������� - ���������� ���� Fx,Fy,Fz,Mx,My,Mz
	*/
	virtual
	void AddPartialBoundary
		(
			int* boundaryNodesIndices, 
			int numberOfNodesInPartialBoundary ,
			int numberOfNodesInFullBoundary ,
			int bcKind,
			double* bcParams
		)
		= 0;

	/** �������� ��������� �������
	* @param boundaryNodesIndices - ������� ��������� �����
	* @param numberOfBoundaryNodes - ���������� ��������� �����
	* @param bcKind - ��� ��������� ������ (=1 - ����, ���������� ����������� �� ������, 
	* �������� �������� �� 6 ��������� - 3 �������������� ���������� ������� ���� 
	* �� ������������ ���� � ��� ������� ������������ ����)
	* @param bcParams - ���������� ��������� �������
	* =1: 6 ���������� - ���������� ���� Fx,Fy,Fz,Mx,My,Mz
	*/
	virtual
	void AddBoundary
		(
			int* boundaryNodesIndices, 
			int numberOfBoundaryNodes,
			int bcKind,
			double* bcParams
		)
		= 0;

	/** �������� ������� ���������� �������
	* @param bcNumber - ����� ���������� ������� (� ������� ����������),
	* @param bcParamNumber - ����� ��������� � ���������� ��������� �������
	* @param bcParamValue - ����� �������� ���������
	*/
	virtual
	void ChangeBoundaryParam
		(
			const int bcNumber,
			const int bcParamNumber,
			const double bcParamValue
		)
		= 0;
	
	/** �������� �������� �� ������ ���������� ������� � �� �������
	* @param bcNumber - ����� ���������� �������
	* @param bcParamNumber - ����� ���������
	*/
	virtual
	double GetBoundaryParam
		(
			const int bcNumber,
			const int bcParamNumber
		)
		= 0;

	virtual
	float UpdateBufferWithOutput() = 0;

	virtual
	void UpdateBuffer() = 0;

	/** ������ ���������� �������
	*/
	virtual
	void InitialSolve
		(
		)
		= 0;


	/** ������ ���������� ��������
	* @param nIteratons - ���������� ��������
	*/
	virtual
	void Solve
		(
			const int nIteratons
		)
		= 0;

	/**
	* ������ ������ ������ ������ �����-�����
	*/
	virtual
	void Solve1()
		= 0;

	/**
	* ������ ������ ������ ������ �����-�����
	*/
	virtual
	void Solve2()
		= 0;

	/**
	* ������ ������� ������ ������ �����-�����
	*/
	virtual
	void Solve3()
		= 0;

	/**
	* ������ ��������� ������ ������ �����-�����
	*/
	virtual
	void Solve4()
		= 0;

	/**
	* ������ ����� ������ ������ �����-�����
	*/
	virtual
	void Solve5()
		= 0;

	double GetData
		(
			int nNode,
			int dir,
			int type
		);
	
	/** �������� ��������� ��������
	* @param dataVector - ������ ��� ������ ���������� ���������
	*/
	virtual
	void GetVectorParameter(float* dataVector);

	/** �������� ��������� ��������
	* @param data - ������ ��� ������ ���������� ���������
	*/
	virtual
	void GetScalarParameter
		(
			float* data
		);

	/** �������� ��������
	* @param data - ������ ��� ������ �������� ��� ���������� ���������
	*/
	virtual
	void GetDisplacement
		(
			float* data
		)
		= 0;

	/** �������� ���������� �� ������ ������ ���������
	* @param data - ������ ��� ������ ���������� ��� ���������� ���������
	*/
	virtual
	void GetStressesByFirstTheoryOfStrength
		(
			float* data
		)
		= 0;

	/** �������� ���������� �� von Mises
	* @param data - ������ ��� ������ ���������� ��� ���������� ���������
	*/
	virtual
	void GetStressesByVonMises
		(
			float* data
		)
		= 0;

	/** �������� ���������� �� X
	* @param data - ������ ��� ������ ���������� ��� ���������� ���������
	*/
	virtual
		void GetStressesX
		(
		float* data
		)
		= 0;

	/** �������� ���������� �� Y
	* @param data - ������ ��� ������ ���������� ��� ���������� ���������
	*/
	virtual
		void GetStressesY
		(
		float* data
		)
		= 0;

	/** �������� ���������� �� Z
	* @param data - ������ ��� ������ ���������� ��� ���������� ���������
	*/
	virtual
		void GetStressesZ
		(
		float* data
		)
		= 0;

	/** �������� ���������� �� XY
	* @param data - ������ ��� ������ ���������� ��� ���������� ���������
	*/
	virtual
		void GetStressesXY
		(
		float* data
		)
		= 0;

	/** �������� ���������� �� XY
	* @param data - ������ ��� ������ ���������� ��� ���������� ���������
	*/
	virtual
		void GetStressesXZ
		(
		float* data
		)
		= 0;

	/** �������� ���������� �� XY
	* @param data - ������ ��� ������ ���������� ��� ���������� ���������
	*/
	virtual
		void GetStressesYZ
		(
		float* data
		)
		= 0;
public:

	float*	_data;						// ����������� ��� ������ � ���� (dataScalar,dataVector)
	float*	_dataVector;				// ����������� ��� ������ � ���� (x,y,z) (��������� � ������� _data)
	int		_dataSize;					// ������ ������ ��� ������
	double* _dataInternal;				// ����������� ��� ������� (X,DX,DDX)
	double* _dataRotationMtx;			// ������� �������� ��� �������
	double* _elementStressFactorCache;	// ������������ ���
	double* _buffer;					// ������ ��� KNC
	double* _coordinates;				// �������� ���������� ���������
	double* _stress;					// ����������
	int _nElements;						// ����� ���������


	// ������ ��� ��������� ������������� ��������
	virtual double  GetElementDisplacement(size_t elementId) const;
	virtual double* GetElementGridCoordinates(size_t elementId) const;
	virtual double* GetElementStress				(size_t elementId) const;
	virtual double* GetElementStressAngular			(size_t elementId) const;
	virtual double* GetElementShift					(size_t elementId) const;
	virtual double* GetElementVelocity				(size_t elementId) const;
	virtual double* GetElementAcceleration			(size_t elementId) const;
	virtual double* GetElementShiftAngular			(size_t elementId) const;
	virtual double* GetElementVelocityAngular		(size_t elementId) const;
	virtual double* GetElementAccelerationAngular	(size_t elementId) const;
	virtual double* GetRotationMatrix				(size_t elementId) const;

	string _uid;
	virtual void SetUid(const string& uid);

	virtual double* GetElementStressFactors(size_t elementId) const;
	virtual void CalculateStrains
		(
		size_t side,			// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
		double *shiftStrains,		// ����� ����������
		double *velocityStrains,	// ����� ���. ���������
		size_t nodeId1,				// ����� ���� 1
		size_t nodeId2					// ����� ���� 2
		) const = 0;
};

