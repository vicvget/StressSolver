#include "MathFunctions.h"


// template class RealHalf

// заполнение константного значения 0.5 для разных типов
#define REAL_HALF_CONSTANT_DEFINITION(Real, constantValue) \
 \
template <> \
const Real RealHalf<Real>::value = constantValue;


// значение константного значения 0.5 для типа float
REAL_HALF_CONSTANT_DEFINITION(float, 0.5f)

// значение константного значения для типа double
REAL_HALF_CONSTANT_DEFINITION(double, 0.5)

// значение константного значения для типа long double
REAL_HALF_CONSTANT_DEFINITION(long double, 0.5L)


#undef REAL_HALF_CONSTANT_DEFINITION


// реализация функции round при использовании MS VC++ 
#if defined(MS_VC_COMPILER)

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
	)
{
	return static_cast<typename RealToInteger<Real>::type>
		(
			std::floor(value + RealHalf<Real>::value)
		);
}


// явное инстанцирование функции округления вещественного числа до ближайшего целого для разных типов
#define ROUND_FUNCTION_INSTANTIATION(Real) \
 \
template \
RealToInteger<Real>::type round<Real> \
	( \
		Real value \
	);


// явное инстанцирование функции округления вещественного числа до ближайшего целого для типа float
ROUND_FUNCTION_INSTANTIATION(float)

// явное инстанцирование функции округления вещественного числа до ближайшего целого для типа double
ROUND_FUNCTION_INSTANTIATION(double)

// явное инстанцирование функции округления вещественного числа до ближайшего целого для типа long double
ROUND_FUNCTION_INSTANTIATION(long double)


#undef ROUND_FUNCTION_INSTANTIATION

#endif