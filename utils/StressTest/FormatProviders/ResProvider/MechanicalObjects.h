#ifndef MECHANICAL_OBJECTS_H

#define MECHANICAL_OBJECTS_H


// модуль для определения интерфейсов классов механических объектов


#include "Types.h"

#include <map>


using std::map;
using std::pair;


// тип объекта
enum class ObjectType
{
	No, // тип отсутствует
	Body, // тело
	ConnectionElement, // соединительный элемент
	RotationMatrix // матрица трансформации
};


// тип движения
enum class MotionType
{
	Linear, // линейное движение
	Angular // угловое движение
};


// тип направления движения
enum class DirectionType
{
	X, // движение вдоль оси X
	Y, // движение вдоль оси Y
	Z // движение вдоль оси Z
};


// степени свободы
enum class FreedomType
{
	NoFreedom, // степень свободы отсутствует
	ShiftX, // перемещение вдоль оси X
	ShiftY, // перемещение вдоль оси Y
	ShiftZ, // перемещение вдоль оси Z
	RotationX, // вращение вокруг оси X
	RotationY, // вращение вокруг оси Y
	RotationZ // вращение вокруг оси Z
};


// величина
enum class Value
{
	Shift, // перемещение
	Velocity, // скорость
	Acceleration // ускорение
};


// компоненты характеристики соединительного элемента
enum class CharacteristicComponent
{
	ElasticForce, // упругая сила
	DampingForce, // демпфирующая сила
	Strain, // относительная деформация
	StrainRate // скорость деформации
};


// Адреса степеней свободы для тела
struct BodyDofOffsetsCache
{
	int PositionOffsets[6];
	int RotationOffset;
};


// степени свободы:
// ключ - номер степени свободы
// значение - смещение
typedef pair<FreedomType, Integer> Freedom;
typedef map<FreedomType, Integer> Freedoms;


// тела
// ключ - номер тела
// значение - список степеней свободы
typedef pair<Integer, Freedoms> Body;
typedef map<Integer, Freedoms> Bodies;


// соединительные элементы
// ключ - номер соединительного элемента
// значение - список степенией свободы
typedef pair<Integer, Freedoms> ConnectionElement;
typedef map<Integer, Freedoms> ConnectionElements;


// кэш степеней свободы для тел:
// ключ - номер тела
// значение - маска степеней свободы
typedef pair<Integer, BodyDofOffsetsCache> BodyDofOffsetsCachePair;
typedef map<Integer, BodyDofOffsetsCache> BodiesDofOffsetsCache;


// матрицы поворота:
// ключ - номер соответствующего тела
// значение - смещение
typedef pair<Integer, Integer> RotationMatrix;
typedef map<Integer, Integer> RotationMatrixes;


/**
* Сформировать степень свободы по типу движения и направлению
* @param motion - тип движения
* @param direction - тип направления движения
* @return сформированная степень свободы
*/
FreedomType ConstructFreedom
	(
		MotionType motion,
		DirectionType direction
	);


#endif // MECHANICAL_OBJECTS_H