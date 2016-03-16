#ifndef OccRectilinearGrid_H

#define OccRectilinearGrid_H

//#pragma GCC visibility push(hidden)
#include "MeshLayer.h"

#include <vector>


using std::vector;


/**
* Класс хранит трехмерное сеточное представление объекта в виде 0 и 1
*/
class OccRectilinearGrid
{
public:

	// Конструкторы и деструктор

	OccRectilinearGrid
		(
			int nl,
			int nr,
			int nel
		);

	OccRectilinearGrid
		(
			OccRectilinearGrid* basicGrid
		);

	~OccRectilinearGrid();

	// Селекторы

	//вернуть размерность сетки
	void GetGridSize
		(
			int *nx,
			int *ny,
			int *nz
		)	const;

	//вернуть начальную точку сетки
	void GetStartPoint
		(
			double *x,
			double *y,
			double *z
		)	const;

	//вернуть размер элемента сетки
	void GetElementSize
		(
			double *x,
			double *y,
			double *z
		)	const;

	//вернуть количество слоев в Grid
	size_t GetNGrid() const;

	/**
	* Вернуть значение элемента сетки
	* @param x - положение элемента по оси Ox
	* @param y - положение элемента по оси Oy
	* @param z - положение элемента по оси Oz
	* @return значение элемента по данным координатам
	*/
	int GetElement
		(
			int x,
			int y,
			int z
		)	const;

	/**
	* Получить размеры сетки (количество элменетов сетки, отложенных вдоль каждой из осей)
	* @param size - массив из 3-х элементов, равных размерам сетки вдоль соответствующих осей
	*/
	void GetSizes
		(
			int sizes[3]
		)	const;

	/**
	* Получить данные о минимальной и максимальной координатах сетки вдоль оси Ox
	* @param minCoordinate - минимальная координата сетки вдоль оси Ox
	* @param maxCoordinate - максимальная координата сетки вдоль оси Ox
	*/
	void GetMinMaxCoordinatesX
		(
			double& minCoordinate,
			double& maxCoordinate
		)	const;

	/**
	* Получить данные о минимальной и максимальной координатах сетки вдоль оси Oy
	* @param minCoordinate - минимальная координата сетки вдоль оси Oy
	* @param maxCoordinate - максимальная координата сетки вдоль оси Oy
	*/
	void GetMinMaxCoordinatesY
		(
			double& minCoordinate,
			double& maxCoordinate
		)	const;

	/**
	* Получить данные о минимальной и максимальной координатах сетки вдоль оси Oz
	* @param minCoordinate - минимальная координата сетки вдоль оси Oz
	* @param maxCoordinate - максимальная координата сетки вдоль оси Oz
	*/
	void GetMinMaxCoordinatesZ
		(
			double& minCoordinate,
			double& maxCoordinate
		)	const;

	/**
	* Получить данные о геометрических размерах сетки
	* @param minCoordinates - массив из 3-х элементов, равных минимальным координатам сетки
	* @param maxCoordinates - массив из 3-х элементов, равных максимальным координатам сетки
	* @param elementSizes - массив из 3-х элементов, равных размерам элемента сетки вдоль соответствующих осей
	*/
	void GetGeometricSizes
		(
			double minCoordinates[3],
			double maxCoordinates[3],
			double elementSizes[3]
		)	const;

	/**
	* Получить положение узла в сетке по его номеру
	* @param nodeNumber - номер узла
	* @param position - найденная позиция узла в сетке
	* @return признак, существует ли узел с данным номером
	*/
	bool GetNodePosition
		(
			int nodeNumber,
			int position[3]
		)	const;

	/**
	* Получить номер узла в сетке по его положению
	* @param position - предполагаемая позиция узла в сетке
	* @param nodeNumber - номер найденного узла
	* @return признак, существует ли узел с данным положением
	*/
	bool GetNodeByPosition
		(
			const int position[3],
			int& nodeNumber
		)	const;

	/**
	* Получить координаты узла по его положению
	* @param position - предполагаемая позиция узла в сетке
	* @param coordinates - координаты найденного узла
	* @return признак, существует ли узел с данным положением
	*/
	void GetNodeCoordinatesByPosition
		(
			const int position[3],
			double coordinates[3]
		)	const;

	/**
	* Получить положение узла сетки по его координатам
	* @param supposedCoordinates - предполагаемые координаты узла
	* @param position - положение узла в сетке
	*/
	void GetNodePositionByCoordinates
		(
			const double supposedCoordinates[3],
			int position[3]
		)	const;

	/**
	* Определить, может ли принадлежать точка с данной координатой
	* ограничивающему параллелепипеду сетки
	* @param direction - индекс направления оси, по которой берется координата
	* @param coordinate - координата вдоль оси
	*/
	bool IsCorrectCoordinate
		(
			int direction,
			double coordinate
		)	const;

	// Перегрузка оператора []
	const MeshLayer& operator[]
		(
			int index
		)	const;

	// Перегрузка оператора []
	MeshLayer& operator[]
		(
			int index
		);


	// Модификаторы

	//установить начальную точку сетки
	void SetStartPoint
		(
			double x,
			double y,
			double z
		);

	//установить размер элемента сетки
	void SetElementSize
		(
			double x,
			double y,
			double z
		);

	//Добавить слой
	void AddLayer
		(
			MeshLayer *layer
		);

	// пронумеровать сетку
	int EnumerateMesh();

	// отменить нумерацию сетки
	int DeenumerateMesh();

	//работа с ячейкой сетки

	//fields

protected:

	//fields
	int nx{};

	int ny{};

	int nz{};

	//начальная точка сетки
	double startPointX{};

	double startPointY{};

	double startPointZ{};

	//размеры элемента
	double dx{};

	double dy{};

	double dz{};

	vector<MeshLayer *> Grid;

	//methods

	/**
	* Вернуть значение элемента сетки
	* (без проверок на попадание в ограничивающий параллелепипед)
	* @param x - положение элемента по оси Ox
	* @param y - положение элемента по оси Oy
	* @param z - положение элемента по оси Oz
	* @return значение элемента по данным координатам
	*/
	int GetElementInternal
		(
			int x,
			int y,
			int z
		)	const;

};

//#pragma GCC visibility pop
#endif