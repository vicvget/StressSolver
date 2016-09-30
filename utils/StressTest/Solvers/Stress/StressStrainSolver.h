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
	string _fileName;
	bool _readIco;
	bool _writeIco;
	int _nWriteIteration;
public:

	size_t vecStride;
	size_t vecStride2;
	size_t matStride;

	StressStrainSolver
		(
			int const nNodes,
			int stride
		);


	void SetZeroVelocities();

	virtual
	~StressStrainSolver();
	
	void InitIco ( const string& fileName,
		bool readIco,
		bool writeIco,
		int nWriteIteration);

	bool ReadIco(const char* fileName);
	void WriteIco(const char* fileName) const;

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
	float UpdateBuffer
		(
			double scale = 1
		)
		= 0;

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

public:

	float*	_data;				// ����������� ��� ������ � ���� (x,y,z,p)
	double* _dataInternal;		// ����������� ��� ������� (X,DX,DDX)
	double* _dataRotationMtx;	// ������� �������� ��� �������
	double* _buffer;			// ������ ��� KNC

	double* _elements;			// �������� ���������� ���������

	int _dataSize;
	double* _stress;			// ����������
	int _nElements;				// ����� ���������


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
};

