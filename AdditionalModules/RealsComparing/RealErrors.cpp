#include "RealErrors.h"


// template class RealErrors

// заполнение констант погрешностей сравнений значениями
#define REAL_ERRORS_CONSTANTS_DEFINITION(Real, absoluteErrorValue, relativeErrorValue) \
 \
template <> \
const Real RealErrors<Real>::absoluteError = absoluteErrorValue; \
 \
template <> \
const Real RealErrors<Real>::relativeError = relativeErrorValue;


// заполнение констант значениями для типа float
REAL_ERRORS_CONSTANTS_DEFINITION(float, 1e-5f, 1e-6f)

// заполнение констант значениями для типа float
REAL_ERRORS_CONSTANTS_DEFINITION(double, 1e-10, 1e-11)

// заполнение констант значениями для типа float
REAL_ERRORS_CONSTANTS_DEFINITION(long double, 1e-15L, 1e-16L)


#undef REAL_ERRORS_CONSTANTS_DEFINITION