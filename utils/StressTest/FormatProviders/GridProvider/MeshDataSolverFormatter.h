#ifndef MeshDataSolverFormatter_H

#define MeshDataSolverFormatter_H


//#pragma GCC visibility push(hidden)
#include "Vertex.h"
#include "MeshDataProvider.h"
#include "NodeLinks.h"


/**
* Класс для создания данных для разделения сетки
*/
class MeshDataSolverFormatter
	:
		public MeshDataProvider
{
public:

	// Конструкторы и деструктор

	/**
	* Конструктор
	* @param gridFileName - наименование файла со структурой сетки
	*/
	MeshDataSolverFormatter
		(
			const string& gridFileName
		);

	/**
	* Конструктор
	* @param grid - структура регулярной ортогональной сетки
	*/
	MeshDataSolverFormatter
		(
			OccRectilinearGrid* grid
		);


	// Селекторы

	/**
	* Получить шаг сетки
	* @return шаг сетки
	*/
	double GetStep() const;
	
	/**
	* Формирует список ссылок на соседние узлы для всех узлов сетки
	* @param links - сформированный список ссылок на соседние узлы для всех узлов сетки
	*/
	void ExportGrid
		(
			Links& links
		)	const;
		
	/**
	* Формирует массив связей и кординат во внешние массивы
	* массивы выделяются в этой же функции
	* поэтому надо следить за памятью и удалять их после того как решатель отработает
	* или скопирует их себе
	*
	* @param mlink - массив связей, хранит по 2 узла: [nodeId1, nodeId2]; [nodeId1, nodeId2]; ...
	* @param uzk - масив координат, хранит по 3 координаты: [x1,y1,z1],[x2,y2,z2],...
	*/
	void ExportGrid
		(
			int*& mlink,
			double*& uzk
		);

	/**
	* Формирует массив связей во внешние массивы
	* массивы выделяются в этой же функции
	*
	* @param mlink - массив связей, хранит по 2 узла: [nodeId1, nodeId2]; [nodeId1, nodeId2]; ...
	*/
	void ExportGrid
		(
			int*& mlink
		);

	/**
	* Формирует массив индексов граничных узлов нагревателей
	* @param mbon - массив граничных узлов нагревателей
	* @param kbon - количество элементов массива
	* @param surfaceIds - фильтр по поверхностям
	*/
	void ExportBoundary
		(
			int* &mbon,
			size_t& kbon,
			const vector<string>& surfaceIds
		)	const;

	void ExportBoundary
		(
			int* &mbon,
			size_t& kbon
		)	const;

	void ExportBoundary
		(
			const vector<string>& surfaceIds,
			vector<int>& boundaryNodesIndices
		)	const;


	// Статические члены
	
	/**
	* Функция создания объекта класса
	*/
	DEFINE_CREATE_OBJECT_MEMBER_FUNCTION(MeshDataSolverFormatter, ReleaseMeshData)

};


#undef DEFINE_CREATE_OBJECT_MEMBER_FUNCTION


//#pragma GCC visibility pop

#endif