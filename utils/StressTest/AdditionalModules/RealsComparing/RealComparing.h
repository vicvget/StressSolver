#ifndef REAL_COMPARING_H

#define REAL_COMPARING_H


#include "RealComparer.h"
#include "../../AdditionalModules/Macros/ArgumentsListMacros.h"


/**
* Класс-обертка над компаратором вещественных чисел разного типа.
* Предоставляет возможность сравнения вещественных чисел с помощью статического компаратора
*/
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences = true // признак, передавать ли параметры в функции сравнения по ссылке
	>
class RealComparing
{
public:


	// Конструкторы и деструктор

	/**
	* Конструктор по умолчанию
	*/
	RealComparing() = delete;

	/**
	* Конструктор копирования
	*/
	RealComparing
		(
			const RealComparing& otherRealComparing
		)
		= delete;

	/**
	* Деструктор
	*/
	~RealComparing() = delete;


	// Типы

	// тип параметра, передаваемого в функции сравнения
	typedef typename RealComparer<Real, UseReferences>::RealArgument RealArgument;


	// Функции класса-обертки

	// Объявления функций класса-обертки над компаратором

	// Объявление функции компаратора
#	define REAL_COMPARING_FUNCTION_DECLARATION(Function, ReturnType, ...) \
	 \
	static \
	ReturnType Function \
		( \
			CREATE_PARAMETERS_LIST(__VA_ARGS__) \
		);

#	define REAL_COMPARING_FUNCTION REAL_COMPARING_FUNCTION_DECLARATION

#	include "RealComparingFunctions.h"

#	undef REAL_COMPARING_FUNCTION

#	undef REAL_COMPARING_FUNCTION_DECLARATION


protected:

	// статический компаратор, используемый для сравнений
	static
	RealComparer<Real, UseReferences> _realComparer;

};


// класс-обертка, предназначенный для сравнения данных вещественных чисел типа float
typedef RealComparing<float, false> FloatComparing;

// класс-обертка, предназначенный для сравнения данных вещественных чисел типа double
typedef RealComparing<double, false> DoubleComparing;

// класс-обертка, предназначенный для сравнения данных вещественных чисел типа long double
typedef RealComparing<long double, true> LongDoubleComparing;


#endif