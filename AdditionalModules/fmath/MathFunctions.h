#ifndef MATH_FUNCTIONS_H

#define MATH_FUNCTIONS_H


#include "../Macros/CompilerFlags.h"

#include <cmath>


/**
* Структура для хранения константного значения 0.5
*/
template
	<
		typename Real // тип вещественного числа
	>
class RealHalf
{
public:

	// константное значение
	static 
	const Real value;

};


// декларация явных специализаций константного значения 0.5
// не работает в Visual C++ (ошибка С2734)
#if !defined(MS_VC_COMPILER)

// декларация явной специализации константного значения 0.5 для разных типов
#define REAL_HALF_CONSTANT_DECLARATION(Real) \
 \
template <> \
const Real RealHalf<Real>::value;


// декларация явной специализации константного значения 0.5 для типа float
REAL_HALF_CONSTANT_DECLARATION(float)

// декларация явной специализации константного значения 0.5 для типа double
REAL_HALF_CONSTANT_DECLARATION(double)

// декларация явной специализации константного значения 0.5 для типа long double
REAL_HALF_CONSTANT_DECLARATION(long double)

#endif


/**
* Задание одностороннего соответствия между вещественнозначными типами и целыми
*/
template
	<
		typename Real // тип вещественного числа
	>
class RealToInteger;


// Специализация задания одностороннего соответствия между вещественнозначными типами и целыми для разных типов
#define REAL_TO_INTEGER_SPECIALIZATION(Real, Integer) \
 \
template <> \
class RealToInteger<Real> \
{ \
public: \
 \
	typedef Integer type; \
 \
};


// Специализация задания одностороннего соответствия между типами float и int
REAL_TO_INTEGER_SPECIALIZATION(float, int)

// Специализация задания одностороннего соответствия между типами double и long long int
REAL_TO_INTEGER_SPECIALIZATION(double, long long int)

// Специализация задания одностороннего соответствия между типами long double и long long int
REAL_TO_INTEGER_SPECIALIZATION(long double, long long int)


#undef REAL_TO_INTEGER_SPECIALIZATION


// работа с функцией std::round
#if !defined(MS_VC_COMPILER)

// в обычном случае просто используем библиотечную функцию std::round
using std::round;

#else

// в случае MS Visual C++ пишем собственную реализацию std::round,
// т.к. она до сих пор не реализована в стандартной библиотеке, предоставляемой Microsoft

/**
* Функция округления вещественного числа до ближайшего целого
* @param value - округляемое значение
* @return - результат округления
*/
template
	<
		typename Real // тип вещественного числа
	>
typename RealToInteger<Real>::type round
	(
		Real value
	);

#endif


#endif // MATH_FUNCTIONS_H