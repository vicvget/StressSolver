#ifndef MESH_DATA_CLEARER_H

#define MESH_DATA_CLEARER_H


#include "NodeLinks.h"
#include "UsingPolicy.h"
#include "OccRectilinearGrid.h"
#include "BoundaryCondition.h"
#include "BoundaryNormal.h"

#include <string>
#include <vector>


using std::string;
using std::vector;


// признак узла, который не принадлежит ни одной подобласти
#define NOT_A_MEMBER -1


// список индексов узлов
typedef vector<int> NodesIndexes;

// список признаков принадлежности узлов той или иной подобласти сетки
typedef vector<int> NodesBelongingList;


/**
* Класс для очистки сетки от изолированных подобластей и/или точек
*/
class MeshDataCleaner  
{
public:

	// Конструкторы и деструктор
	
	/**
	* Конструктор по умолчанию
	*/
	MeshDataCleaner();

	/**
	* Деструктор
	*/
	~MeshDataCleaner();


	// Модификаторы

	/**
	* Установить наименование файла со структурой сетки
	* @param gridFileName - наименование файла со структурой сетки
	*/
	void SetGridFileName
		(
			const string& gridFileName
		);

	/**
	* Установить тип политики удаления изолированных подобластей сетки
	* @param deletionPolicyType - устанавливаемый тип политики удаления подобластей сетки
	*/
	void SetDeletionPolicyType
		(
			DeletionPolicyType deletionPolicyType
		);

	/**
	* Очистить сетку от изолированных областей
	* @return признак успешной (true) или неуспешной (false) очистки
	*/
	bool CleanGrid() const;


protected:

	// наименование файла со структурой сетки
	string _gridFileName;

	// тип политики удаления изолированных подобластей сетки
	DeletionPolicyType _deletionPolicyType;


	// Функции очистки сетки

	/**
	* Очистить сетку от изолированных областей
	* @param deletionPolicyType - тип политики удаления изолированных подобластей сетки
	* @param grid - сетка, для которой производится очистка
	* @param boundaryConditions - список граничных условий сетки, для которой производится очистка
	*/
	static
	void CleanGrid
		(
			DeletionPolicyType deletionPolicyType,
			OccRectilinearGrid* grid,
			BoundaryConditions& boundaryConditions,
			vector<BoundaryNormal*>& boundaryNormals
		);


	// Функции по определению несвязных подобластей сетки

	/**
	* Найти несвязанные подобласти сетки
	* @param links - список ссылок на соседние узлы для всех узлов сетки
	* @param nodesBelongingList - список признаков принадлежности узлов той или иной подобласти сетки
	* @return - количество найденных подобластей сетки
	*/
	static
	int FindGridSubdomains
		(
			const Links& links,
			NodesBelongingList& nodesBelongingList
		);

	/**
	* Найти свободный (не принадлежащий ни одной подобласти) узел в сетке
	* @param nodesBelongingList - список признаков принадлежности узлов той или иной подобласти сетки
	* @param nodeIndex - индекс найденного узла
	* @return - признак, найден ли узел (true) или не найден (false)
	*/
	static
	bool FindFreeNode
		(
			NodesBelongingList& nodesBelongingList,
			int& nodeIndex
		);

	/**
	* Заполнить очередную подобласть сетки
	* @param links - список ссылок на соседние узлы для всех узлов сетки
	* @param nodesBelongingList - список признаков принадлежности узлов той или иной подобласти сетки
	* @param subdomainIndex - номер подобласти, которую необходимо заполнить
	* @param startNodeIndex - индекс узла, с которого следует начать заполнение подобласти
	*/
	static
	void FillSubdomain
		(
			const Links& links,
			NodesBelongingList& nodesBelongingList,
			int subdomainIndex,
			int startNodeIndex
		);

	/**
	* Найти размеры подобластей сетки
	* @param nodesBelongingList - список признаков принадлежности узлов той или иной подобласти сетки
	* @param subdomainsCount - количество подобластей сетки
	* @param subdomainsSizesList - сформированный массив размеров подобластей сетки
	*/
	static
	void FindSubdomainsSizes
		(
			const NodesBelongingList& nodesBelongingList,
			int subdomainsCount,
			SubdomainsSizesList& subdomainsSizesList
		);


	// Применить политику удаления

	/**
	* Применить политику удаления изолированных подобластей
	* @param applyingPolicyType - применяемый тип политики удаления изолированных подобластей
	* @param sizesList - список размеров изолированных подобластей сетки
	* @param belongingList - список флагов принадлежности изолированных подобластей чистой сетке
	*/
	static
	void ApplyDeletionPolicy
		(
			DeletionPolicyType applyingPolicyType,
			const SubdomainsSizesList& sizesList,
			SubdomainsBelongingList& belongingList
		);


	// Преобразовать сетку и граничные подобласти, оставив только те узлы,
	// которые принадлежат отмеченным подобластям

	/**
	* Оставить в сетке только те элементы, которые принадлежат отмеченным подобластям
	* @param grid - сетка, над которой производится преобразование
	* @param nodesBelongingList - массив принадлежности узлов той или иной подобласти сетки
	* @param subdomainsBelongingList - список признаков принадлежности подобластей сетке
	*/
	static
	void ConvertGrid
		(
			OccRectilinearGrid* grid,
			const NodesBelongingList& nodesBelongingList,
			const SubdomainsBelongingList& subdomainsBelongingList
		);

	/**
	* Найти новые значения индексов узлов сетки
	* @param grid - сетка, над которой производится преобразование
	* @param newNodesIndexes - список новых индексов сетки
	*/
	static
	void FindNewNodesIndexes
		(
			const OccRectilinearGrid* grid,
			NodesIndexes& newNodesIndexes
		);

	/**
	* Оставить в граничном условии сетки только те элементы, которые принадлежат отмеченным подобластям
	* @param boundaryCondition - граничное условие, над которым производится преобразование
	* @param newNodesIndexes - список новых индексов сетки
	*/
	static
	void ConvertBoundaryCondition
		(
			BoundaryCondition*& boundaryCondition,
			const NodesIndexes& newNodesIndexes
		);

	/**
	* Оставить в граничных условиях сетки только те элементы, которые принадлежат отмеченным подобластям
	* @param boundaryConditions - список граничных условий, над которым производится преобразование
	* @param newNodesIndexes - список новых индексов сетки
	*/
	static
	void ConvertBoundaryConditions
		(
			BoundaryConditions& boundaryConditions,
			const NodesIndexes& newNodesIndexes
		);

	/**
	* Оставить в граничных нормалях сетки только те элементы, которые принадлежат отмеченным подобластям, а для оставленных пересчитать индексы
	* @param boundaryConditions - список граничных условий, над которым производится преобразование
	* @param newNodesIndexes - список новых индексов сетки
	*/
	static void ConvertBoundaryNormals(vector<BoundaryNormal*> &boundaryNormals, const NodesIndexes& newNodesIndexes);
};


#endif