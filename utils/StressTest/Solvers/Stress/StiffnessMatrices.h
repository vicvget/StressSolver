#ifndef STIFFNESS_MATRICES_H

#define STIFFNESS_MATRICES_H

#include "StiffnessMatrix.h"
#include "StiffnessMatrixPackedInCoordinateFormat.h"
#include "StiffnessMatrixPackedInCSCFormat.h"
#include "../../AdditionalModules/Macros/ArgumentsListMacros.h"

#include <tuple>


/**
* Класс кортежа матриц жесткости
*/
template
	<
		class... Matrices // классы матриц жесткости
	>
class StiffnessMatricesTuple
{
public:

	// Типы

	// Синоним для типа используемого кортежа - std::tuple
	using Tuple = std::tuple<Matrices...>;


	// Конструкторы и деструктор

	/**
	* Конструктор
	* @param nodesCount количество узлов сеточного представления тела
	*/
	explicit
	StiffnessMatricesTuple
		(
			size_t nodesCount
		);

	/**
	* Конструктор перемещения
	* @param otherStiffnessMatricesTuple - другой кортеж матриц жесткости, из которого производится перемещение
	*/
	StiffnessMatricesTuple
		(
			StiffnessMatricesTuple&& otherStiffnessMatricesTuple
		);


	// Селекторы

	/**
	* Получить кортеж матриц жесткости
	* @return кортеж матриц жесткости
	*/
	const Tuple& GetTuple() const;

	/**
	* Получить кортеж матриц жесткости
	* @return кортеж матриц жесткости
	*/
	Tuple& GetTuple();


	// Методы кортежа матриц жесткости

	// Объявление метода кортежа матриц жесткости
#	define STIFFNESS_MATRICES_FUNCTION_DECLARATION(Function, ReturnType, ...) \
	 \
	ReturnType Function \
		( \
			CREATE_PARAMETERS_LIST(__VA_ARGS__) \
		);

#	define STIFFNESS_MATRICES_FUNCTION STIFFNESS_MATRICES_FUNCTION_DECLARATION

#	include "StiffnessMatricesFunctions.h"

#	undef STIFFNESS_MATRICES_FUNCTION

#	undef STIFFNESS_MATRICES_FUNCTION_DECLARATION


private:

	// кортеж матриц жесткости
	Tuple _tuple;

};


// используемые классы матриц жесткости
#define USED_STIFFNESS_MATRICES \
	/* StiffnessMatrix, */ \
	/* StiffnessMatrixPackedInCoordinateFormat, */ \
	StiffnessMatrixPackedInCSCFormat


// Явная конкретизация StiffnessMatrices

// для типов StiffnessMatrix и StiffnessMatrixPackedInCoordinateFormat
extern
template class StiffnessMatricesTuple<USED_STIFFNESS_MATRICES>;


// синоним для StiffnessMatricesTuple, специализированного используемыми классами матриц жесткости
using StiffnessMatrices = StiffnessMatricesTuple<USED_STIFFNESS_MATRICES>;


#endif // STIFFNESS_MATRICES_H