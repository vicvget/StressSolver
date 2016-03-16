#ifndef MESHDATAGRAPHFORMATTER_H
#define MESHDATAGRAPHFORMATTER_H
//#pragma GCC visibility push(hidden)
#include "Vertex.h"
#include "BoundaryCondition.h"
#include "MeshDataProvider.h"

/*
*Класс для создания данных для разделения сетки
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
	* Функция преобразовывает сетку в массив указателей объектов типа Vertex
	* @return возвращет 1 если преобразование успешно
	*/
	int Vectorize();
	/**
	* Функция возвращает сетку
	* @return вектор указателей на объекты типа Vertex
	*/
	vector<Vertex *> GetArrayOfVertex();
	/**
	* Функция возвращает граничные условия для сетки
	* @return вектор указателей на объекты типа BoundaryCondition
	*/
	vector<BoundaryCondition *> GetArrayOfBC();
};
//#pragma GCC visibility pop
#endif