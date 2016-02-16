#include "RealComparing.h"


// template class RealComparing

// Реализация функций компаратора

// Реализация функции компаратора
#define REAL_COMPARING_FUNCTION_IMPLEMENTATION(Function, ReturnType, ...) \
 \
template \
	< \
		typename Real, \
		bool UseReferences \
	> \
ReturnType RealComparing<Real, UseReferences>::Function \
	( \
		CREATE_PARAMETERS_LIST(__VA_ARGS__) \
	) \
{ \
	return _realComparer.Function \
		( \
			CREATE_ARGUMENTS_LIST(__VA_ARGS__) \
		); \
}

#define REAL_COMPARING_FUNCTION REAL_COMPARING_FUNCTION_IMPLEMENTATION

#include "RealComparingFunctions.h"

#undef REAL_COMPARING_FUNCTION

#undef REAL_COMPARING_FUNCTION_IMPLEMENTATION

// статический компаратор, используемый для сравнений
template
	<
		typename Real, // тип вещественного числа
		bool UseReferences // признак, передавать ли параметры в функции сравнения по ссылке
	>
// static
RealComparer<Real, UseReferences> RealComparing<Real, UseReferences>::_realComparer;


// Явное инстанцирование класса RealComparing для вещественнозначных типов

// типы, использующие передачу параметров по ссылке в функциях сравнения

// для типа float
template
class RealComparing<float, true>;

// для типа double
template
class RealComparing<double, true>;

// для типа long double
template
class RealComparing<long double, true>;


// типы, использующие передачу параметров по значению в функциях сравнения

// для типа float
template
class RealComparing<float, false>;

// для типа double
template
class RealComparing<double, false>;

// для типа long double
template
class RealComparing<long double, false>;