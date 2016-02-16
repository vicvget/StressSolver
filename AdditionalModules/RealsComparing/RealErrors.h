#ifndef REAL_ERRORS_H

#define REAL_ERRORS_H


#include "../../AdditionalModules/Macros/CompilerFlags.h"


/**
* Класс, содержащий значения минимальных погрешностей сравнения
* для разных типов вещественных чисел
*/
template
	<
		typename Real // тип вещественного числа
	>
class RealErrors
{
public:

	// минимальная абсолютная погрешность сравнений
	static
	const Real absoluteError;

	// минимальная относительная погрешность сравнений
	static
	const Real relativeError;


protected:

	RealErrors();

	RealErrors(const RealErrors&);

	~RealErrors();

};


// декларация явных специализаций погрешностей сравнений
// не работает в Visual C++ (ошибка С2734)
#if !defined(MS_VC_COMPILER)

// декларация явных специализаций погрешностей сравнений для разных типов
#define REAL_ERRORS_CONSTANTS_DECLARATION(Real) \
 \
template <> \
const Real RealErrors<Real>::absoluteError; \
 \
template <> \
const Real RealErrors<Real>::relativeError;


// декларация явных специализаций погрешностей сравнений для типа float
REAL_ERRORS_CONSTANTS_DECLARATION(float)

// декларация явных специализаций погрешностей сравнений для типа double
REAL_ERRORS_CONSTANTS_DECLARATION(double)

// декларация явных специализаций погрешностей сравнений для типа long double
REAL_ERRORS_CONSTANTS_DECLARATION(long double)


#undef REAL_ERRORS_CONSTANTS_DECLARATION

#endif


#endif // REAL_ERRORS_H