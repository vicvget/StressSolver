#ifndef GRID_SELECTED_AREA_H

#define GRID_SELECTED_AREA_H


#include "MeshDataProvider.h"

#include <vector>
#include <array>


/**
* Выделенная область сетки в виде параллелепипеда
*/
class GridSelectedArea final
{
public:

	// Типы данных

	// список индексов узлов сетки
	using PointList = std::vector<int>;

	// массив из трех списков индексов узлов
	using PointLists = std::array<PointList, 3>;


	// Сформировать выделенную область сетки

	/**
	* Выделить область сетки в виде паралелепипеда с заданным центром и размерами
	* @param meshDataProvider - провайдер доступа к данным сетки
	* @param centerPosition - положение центра выделяемой области
	* @param areaSizes - размеры выделяемой области (массив длин ребер параллелепипеда) в узлах сетки
	*/
	void Form
		(
			const MeshDataProvider* meshDataProvider,
			const int centerPosition[3],
			const int areaSizes[3]
		);


	// Селекторы

	/**
	* Получить список узлов, входящих в выделенную область
	* @return список узлов, входящих в выделенную область
	*/
	const PointList& GetAreaPoints() const;


	// Запись в файл / чтение из файла

	/**
	* Записать файл с выделенной областью сетки (в виде параллелепипеда)
	* @param selectedAreaFileName - наименование файла с выделенной областью сетки
	* @return признак успешной (true) или неуспешной (false) записи файла
	*/
	bool ToFile
		(
			const std::string& selectedAreaFileName
		)	const;

	/**
	* Считать файл с выделенной областью сетки (в виде параллелепипеда)
	* @param selectedAreaFileName - наименование файла с выделенной областью сетки
	* @return признак успешного (true) или неуспешного (false) чтения файла
	*/
	bool FromFile
		(
			const std::string& selectedAreaFileName
		);

	/**
	* Создать файл с выделенной областью сетки (в виде параллелепипеда)
	* @param selectedAreaFileName - наименование файла с выделенной областью сетки
	* @param meshDataProvider - провайдер доступа к данным сетки
	* @param centerPosition - положение центра выделенной области в сетке
	* @param areaSizes - размеры выделенной области (массив длин ребер параллелепипеда) в узлах сетки
	* @return признак успешного (true) или неуспешного (false) создания файла
	*/
	static
	bool CreateFile
		(
			const std::string& selectedAreaFileName,
			const MeshDataProvider* meshDataProvider,
			const int centerPosition[3],
			const int areaSizes[3]
		);


private:

	// список узлов, входящих в выделенную область
	PointList _areaPoints;

	// списки узлов, входящих в подобласти выделенной области,
	// являющихся пересечением координатных плоскостей с данной областью
	PointLists _subAreasPoints;


	// Вспомогательные функции формирования выделенной области

	/**
	* Выделить область сетки в виде паралелепипеда с заданным центром и размерами
	* @param meshDataProvider - провайдер доступа к данным сетки
	* @param centerPosition - положение центра выделяемой области
	* @param areaSizes - размеры выделяемой области (массив длин ребер параллелепипеда) в узлах сетки
	* @param selectedArea - список номеров точек сетки, составляющих выделяемую область
	* @param selectedSubAreas - массив списков номеров точек, являющихся сечениями выделяемой области,
	* перпендикулярными координатным осям и проходящими через центр выделяемой области
	*/
	static
	void SelectArea
		(
			const MeshDataProvider* meshDataProvider,
			const int centerPosition[3],
			const int areaSizes[3],
			PointList& selectedArea,
			PointLists& selectedSubAreas
		);


	// Вспомогательные функции записи в файл / чтения из файла

	/**
	* Записать в файл список номеров узлов сетки
	* @return признак успешной (true) или неуспешной (false) записи файла
	*/
	static
	bool WritePointList
		(
			ofstream& outputFile,
			const PointList& pointList
		);

	/**
	* Считать в файл список номеров узлов сетки
	* @return признак успешного (true) или неуспешного (false) чтения файла
	*/
	static
	bool ReadPointList
		(
			ifstream& inputFile,
			PointList& pointList
		);

	/**
	* Записать файл с выделенной областью сетки (в виде параллелепипеда)
	* @param selectedAreaFileName - наименование файла с выделенной областью сетки
	* @param selectedArea - список номеров точек сетки, составляющих выделенную область
	* @param selectedSubAreas - массив списков номеров точек, являющихся сечениями выделенной области,
	* перпендикулярными координатным осям и проходящими через центр выделенной области
	* @return признак успешной (true) или неуспешной (false) записи файла
	*/
	static
	bool WriteFile
		(
			const std::string& selectedAreaFileName,
			const PointList& selectedArea,
			const PointLists& selectedSubAreas
		);

	/**
	* Считать файл с выделенной областью сетки (в виде параллелепипеда)
	* @param selectedAreaFileName - наименование файла с выделенной областью сетки
	* @param selectedArea - список номеров точек сетки, составляющих выделенную область
	* @param selectedSubAreas - массив списков номеров точек, являющихся сечениями выделенной области,
	* перпендикулярными координатным осям и проходящими через центр выделенной области
	* @return признак успешного (true) или неуспешного (false) чтения файла
	*/
	static
	bool ReadFile
		(
			const std::string& selectedAreaFileName,
			PointList& selectedArea,
			PointLists& selectedSubAreas
		);

};


#endif // GRID_SELECTED_AREA_H