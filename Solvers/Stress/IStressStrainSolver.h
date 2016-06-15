#pragma once
#include "IMemory.h"

/** ������������ ����� ��� �������� ���
*/
class IStressStrainSolver: public IMemory
{
public:
	virtual ~IStressStrainSolver() = 0;

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

	/** ���������� �������� ������ � �������� ���������������� ��������
	* @param scale - ������� ��������
	*/
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


protected:

	float* _data;
	double* _dataInternal;
	int _dataSize;
	double* _stress;
	int _nNodes;

};

#endif