#pragma once
#include <iostream>

using std::cout;
using std::string;

enum DataType
{
	DT_Shifts=0,
	DT_Speeds=1,
	DT_Accelerations=2
};

#define ALIGNED_MEM
namespace stress
{
#define MAX(x, y) ((x) > (y) ? x : y)
	double DMOD(const double d1, const double d2);

//#ifdef ALIGNED_MEM
//	const size_t vecStride = 4;				// �������� ��������
//#else
//	const size_t vecStride = 3;				// �������� ��������
//#endif
//
//	const size_t matStride = vecStride * 3; // �������� ������ �������� (3x3, 3x4)
//	const size_t vecStride2 = vecStride * 2; // �������� �������� ����������� ��� ����((3+0)X, (3+0)R)
//	const size_t alignment = 32;
}

#define MAX(x, y) ((x) > (y) ? x : y)
double DMOD(const double d1, const double d2);

//#ifdef ALIGNED_MEM
//const size_t vecStride = 4;				// �������� ��������
//#else
//const size_t vecStride = 3;				// �������� ��������
//#endif

//const size_t matStride = vecStride * 3; // �������� ������ �������� (3x3, 3x4)
//const size_t vecStride2 = vecStride * 2; // �������� �������� ����������� ��� ����((3+0)X, (3+0)R)
const size_t alignment = 32;


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
	void UpdateBuffer
		(
			double scale = 1
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


public:

	float*	_data;				// ����������� ��� ������ � ���� (x,y,z,p)
	double* _dataInternal;		// ����������� ��� ������� (X,DX,DDX)
	double* _dataRotationMtx;	// ������� �������� ��� �������

	double* _elements;			// �������� ���������� ���������

	int _dataSize;
	double* _stress;			// ����������
	int _nElements;				// ����� ���������


	virtual double* GetElementGridCoordinates(int elementId) const;;
	virtual double* GetElementStress(int elementId) const;
	virtual double* GetElementStressAngular(int elementId) const;
	virtual double* GetElementShift(int elementId) const;
	virtual double* GetElementVelocity(int elementId) const;
	virtual double* GetElementAcceleration(int elementId) const;
	virtual double* GetElementShiftAngular(int elementId) const;
	virtual double* GetElementVelocityAngular(int elementId) const;
	virtual double* GetElementAccelerationAngular(int elementId) const;
	virtual double* GetRotationMatrix(int elementId) const;
};

