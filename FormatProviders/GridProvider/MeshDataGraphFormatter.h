#ifndef MESHDATAGRAPHFORMATTER_H
#define MESHDATAGRAPHFORMATTER_H
//#pragma GCC visibility push(hidden)
#include "Vertex.h"
#include "BoundaryCondition.h"
#include "MeshDataProvider.h"

/*
*����� ��� �������� ������ ��� ���������� �����
*/
class MeshDataGraphFormatter:public MeshDataProvider
{
protected:
	vector<Vertex *> graph;
	vector<BoundaryCondition *> bc;
public:
	MeshDataGraphFormatter(OccRectilinearGrid* newgrid);
	MeshDataGraphFormatter(OccRectilinearGrid *newgrid,vector<BoundaryCondition *> Surfaces);
	~MeshDataGraphFormatter(void);
	/**
	* ������� ��������������� ����� � ������ ���������� �������� ���� Vertex
	* @return ��������� 1 ���� �������������� �������
	*/
	int Vectorize();
	/**
	* ������� ���������� �����
	* @return ������ ���������� �� ������� ���� Vertex
	*/
	vector<Vertex *> GetArrayOfVertex();
	/**
	* ������� ���������� ��������� ������� ��� �����
	* @return ������ ���������� �� ������� ���� BoundaryCondition
	*/
	vector<BoundaryCondition *> GetArrayOfBC();
};
//#pragma GCC visibility pop
#endif