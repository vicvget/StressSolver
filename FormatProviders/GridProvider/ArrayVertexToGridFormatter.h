#ifndef ARRAYVERTEXTOGRIDFORMATTER_H
#define ARRAYVERTEXTOGRIDFORMATTER_H
//#pragma GCC visibility push(hidden)
#include "MeshDataProvider.h"
#include "BoundaryCondition.h"
#include "Vertex.h"

/*
Класс слежит для преобразования вектора объектов типа 
Vertex в сетку хранящуюся в объекте типа OccRectilinearGrid
*/
class ArrayVertexToGridFormatter
	:
		public MeshDataProvider
{
protected:

	vector<Vertex*> _array;
	vector<BoundaryCondition*> _bc;
	//начальная точка сетки
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
	* Функция преобразовывает массив указателей объектов типа Vertex в сетку
	* @return возвращет 1 если преобразование успешно
	*/
	int Transform();

private:

	/**
	* Функция возвращает размер сетки
	* @param sx 
	* @param sy  размеры куба ограничивающего объект по осям
	* @param sz 
	* @return возвращет 1 если размер найден
	*/
	int GetBBox(double *sx,double *sy, double *sz);

	/**
	* Функция возвращает размер элемента сетки
	* @return возвращет 1 если размер элемента найден
	*/
	int GetSize();

};

//#pragma GCC visibility pop
#endif