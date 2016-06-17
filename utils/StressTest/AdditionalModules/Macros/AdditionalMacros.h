#pragma once
#include<cmath>
namespace DoubleComparing
{
	template <class Real>
	bool IsZero(Real number)
	{
		if (number == static_cast<Real>(0.0L))
		{
			return true;
		}

		return fabs(number) < 1e-10;
	}
}

#include "ArgumentsCountMacro.h"
#include "ControlMacros.h"
#include "StringMacros.h"


// вызвать макрос с заданными параметрами
#define INVOKE_MACRO(macro, ...) \
	macro(__VA_ARGS__)


// развернуть параметр макроса
#define EXPAND_SINGLE(x) \
	x


// развернуть все параметры макроса
#define EXPAND_MULTIPLE(...) \
	__VA_ARGS__


// объединить параметры макроса
#define COMBINE(a, b) \
	a b


// добавить запятую перед аргументом
#define ADD_COMMA_SINGLE(x) \
	, x

// добавить запятую перед списком аргументов
#define ADD_COMMA_MULTIPLE(...) \
	STATEMENT_IF(GET_ARGUMENTS_COUNT(__VA_ARGS__), COMMA) __VA_ARGS__

