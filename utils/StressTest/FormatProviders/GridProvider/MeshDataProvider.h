#ifndef MeshDataProvider_H

#define MeshDataProvider_H


//#pragma GCC visibility push(hidden)
#include "OccRectilinearGrid.h"
#include "SurfaceParam.h"
#include "BoundaryCondition.h"
#include "BoundaryNormal.h"
#include "../GridProvider/RLCHeader.h"
//#include "../ProviderFrm/FrundFacade/MphParamsSet/FreeSolverGridParams.h"
#include "../../AdditionalModules/AuxiliaryModules/SmartPointers.h"

#include <fstream>


using std::ifstream;


/** Класс для получения данных о сетке
*/
class MeshDataProvider
{
public:

	// Конструкторы и деструктор

	MeshDataProvider() = default;
	
	MeshDataProvider
		(
			OccRectilinearGrid* grid
		);

	MeshDataProvider
		(
			OccRectilinearGrid* grid,
			const vector<BoundaryCondition*>& surfaces
		);

	MeshDataProvider
		(
			OccRectilinearGrid* grid,
			const vector<BoundaryCondition*>& surfaces,
			const vector<BoundaryNormal*>& normals
		);

	MeshDataProvider
		(
			const string& gridFile
		);

	virtual
	~MeshDataProvider() = default;


	// Селекторы

	//
	OccRectilinearGrid* GetGrid();

	//
	const vector<BoundaryCondition*>& GetBC() const;

	// Получить нормали
	const vector<BoundaryNormal*>& GetNormals() const;

	// Получить грид парамс для свободного решателя
//	FreeSolverGridParams* GetFreeSolverGridParams() const;

	//
	bool IsGridLoaded() const;

	//
	int GetCount() const;

	//
	int GetNodesNumber() const;

	/**
	* Получить количество точек в ограничивающем сетку параллелепипеде
	* @return количество точек в ограничивающем сетку параллелепипеде
	*/
	int GetPointsCountInBox() const;

	/**
	* Получить количество граничных поверхностей в сетке
	* @return количество граничных поверхностей в сетке
	*/
	size_t GetBondaryConditionsCount() const;

	/**
	* Получить размеры сетки в узлах по каждому из координатных направлений
	* @param xSize - размер сетки вдоль оси Ox
	* @param ySize - размер сетки вдоль оси Oy
	* @param zSize - размер сетки вдоль оси Oz
	*/
	void GetGridSize
		(
			int& xSize,
			int& ySize,
			int& zSize
		)	const;

	/**
	* Получить координаты начальной точки сетки
	* @param xStartPoint - координата начальной точки по оси Ox
	* @param yStartPoint - координата начальной точки по оси Oy
	* @param zStartPoint - координата начальной точки по оси Oz
	*/
	void GetStartPoint
		(
			double& xStartPoint,
			double& yStartPoint,
			double& zStartPoint
		)	const;

	/**
	* Получить размер элемента (узла) сетки по каждому из координатных направлений
	* @param xElementSize - размер элемента (узла) сетки вдоль оси Ox
	* @param yElementSize - размер элемента (узла) сетки вдоль оси Oy
	* @param zElementSize - размер элемента (узла) сетки вдоль оси Oz
	*/
	void GetElementSize
		(
			double& xElementSize,
			double& yElementSize,
			double& zElementSize
		)	const;

	/**
	* Получить номер элемента сетки с данными координатами
	* @param xCoordinate - координата элемента по оси Ox
	* @param yCoordinate - координата элемента по оси Oy
	* @param zCoordinate - координата элемента по оси Oz
	*/
	int GetPoint
		(
			int xCoordinate,
			int yCoordinate,
			int zCoordinate
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


	// Модификаторы

	//нумерация элементов для решателя
	int EnumerationMesh();

	/**
	* Освободить память, выделенную под данные сетки
	*/
	void ReleaseMeshData();

	/** Загружает нормали
	* @ifs - входной поток
	*/
	void LoadNormals(ifstream &ifs);

	/** Загружает грид парамс для свободного решателя
	* @ifs - входной поток
	*/
	//void LoadFreeSolverGridParams(ifstream &ifs);


	// Статические члены
	
	/**
	* Функция создания объекта класса
	*/
	DEFINE_CREATE_OBJECT_MEMBER_FUNCTION(MeshDataProvider, ReleaseMeshData)


protected:

	//количество элементов в сетке
	int _nodesCount{};

	//размеры элемента
	double _dx{};
	double _dy{};
	double _dz{};

	OccRectilinearGrid* _grid{};

	vector<BoundaryCondition*> _surfaces;

	// вектор нормалей
	vector<BoundaryNormal*> _normals;

	//FreeSolverGridParams* _freeSolverGridParams{};

};


//#pragma GCC visibility pop
#endif