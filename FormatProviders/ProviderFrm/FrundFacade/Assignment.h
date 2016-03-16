#ifndef ASSIGNMENT_H

#define ASSIGNMENT_H


#include "FModel.h"

#include "../../../FormatProviders/ResProvider/OldAddresses.h"

#include <string>
#include <map>

// All members of facade classes should be public !!!

typedef pair<int, FJointChar> JointCharPair;
typedef map<int, FJointChar> JointChars;

typedef vector<FJoint> Joints;

/**
* Класс назначения телу геометрии в виде заранее подготовленной сетки
* и связанных с этим назначением формированием файла "sforces.dat"
*
* @author Izmailov Timur
*/
class Assignment
{
public:

	Assignment();

	~Assignment();
	
	/** Назначение телу новой геометрии
	* @param geoFileName - наименование файла ".geo", содержащего сетку разбиения тела
	* @param bodyNumber - номер тела для которого назначается новая геометрия в виде сетки
	* @return 0, если успешно
	*/
	int Assign
		(
			const string& geoFileName,
			int bodyNumber
		);


private:

	// модель, в которой все назначается
	//FModel _model;
	FGeometry _geometry;

	// Новая геометрия тела (хранится в "*.geo"-файле)
	FGeometry* newGeometry;

	// Тело, для которого назначается новая геометрия
	FBody _body;

	// Класс для доступа к адресам элементов модели
	OldAddresses oldAddresses;


	/** Поиск файла с моделью
	* @param modelFileName - наименование файла ".frm" модели
	* @return 0, если успешно
	*/
	static
	int FindModelFile
		(
			string& modelFileName
		);

	/** Инициация данного тела и соединительных элементов (с их характеристиками)
	* @param modelFileName - наименование файла ".frm" модели
	* @return 0, если успешно
	*/
	int Init
		(
			const string& modelFileName,
			int bodyNumber,
			JointChars& jointChars,
			Joints& joints
		);

	/** Нахождение точки, ближайшей к данной
	* @param point - заданная точка
	* @patam nearestPoint - ближайшая точка
	*/
	void FindNearestPoint
		(
			const FGeometryPoint& point,
			int& nearestPointNumber,
			FGeometryPoint& nearestPoint
		)	const;

	/** Преобразование кооординат точки по матрице поворота и вектору сдвига
	* @param point - координаты точки до преобразования
	* @param transformMatrix - вектор сдвига + матрица поворота
	* @param transformedPoint - координаты точки после преобразования
	* @return 0, если успешно
	*/
	static
	void TransformCoordinates
		(
			const FGeometryPoint& cmNode,
			const FGeometryPoint& point,
			const float (&transformMatrix)[12],
			FGeometryPoint& transformedPoint
		);

	/** Нахождение точки, ближайшей к данной
	* @param point - заданная точка
	* @patam nearestPoint - ближайшая точка
	*/
	int FormSForces
		(
			const JointChars& jointChars,
			const Joints& joints
		);

};


#endif // ASSIGNMENT_H