#ifndef ADDITIONAL_MACROS_H

#define ADDITIONAL_MACROS_H


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


#endif // ADDITIONAL_MACROS_H