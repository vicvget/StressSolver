#ifndef ARRAYVERTEXTOGRIDFORMATTER_H
#define ARRAYVERTEXTOGRIDFORMATTER_H
//#pragma GCC visibility push(hidden)
#include "MeshDataProvider.h"
#include "BoundaryCondition.h"
#include "Vertex.h"

/*
����� ������ ��� �������������� ������� �������� ���� 
Vertex � ����� ���������� � ������� ���� OccRectilinearGrid
*/
class ArrayVertexToGridFormatter
	:
		public MeshDataProvider
{
protected:

	vector<Vertex*> _array;
	vector<BoundaryCondition*> _bc;
	//��������� ����� �����
	double startPointX{};
	double startPointY{};
	double startPointZ{};
	//
	int nx{};
	int ny{};
	int nz{};

public:

	ArrayVertexToGridFormatter
		(
			const vector<Vertex*>& myArray
		);

	ArrayVertexToGridFormatter
		(
			const vector<Vertex*>& Array,
			const vector<BoundaryCondition*>& BC
		);

	/**
	* ������� ��������������� ������ ���������� �������� ���� Vertex � �����
	* @return ��������� 1 ���� �������������� �������
	*/
	int Transform();

private:

	/**
	* ������� ���������� ������ �����
	* @param sx 
	* @param sy  ������� ���� ��������������� ������ �� ����
	* @param sz 
	* @return ��������� 1 ���� ������ ������
	*/
	int GetBBox(double *sx,double *sy, double *sz);

	/**
	* ������� ���������� ������ �������� �����
	* @return ��������� 1 ���� ������ �������� ������
	*/
	int GetSize();

};

//#pragma GCC visibility pop
#endif