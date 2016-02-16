#ifndef ARGUMENTS_COUNT_MACRO_H

#define ARGUMENTS_COUNT_MACRO_H


#include "CompilerFlags.h"
#include "AdditionalMacros.h"


// последовательность из 10-ти первых натуральных чисел, записанных в обратном порядке, и нуля
#define SEQUENCE_N() \
	10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0


#if !defined(MS_VC_COMPILER)

	// реализация макроса получения количества параметров макроса
#	define GET_ARGUMENTS_COUNT_IMPLEMENTATION(...) \
		GET_NTH_ARGUMENT(__VA_ARGS__)

#else

	// реализация макроса получения количества параметров макроса
	// (специально для Microsoft Visual C++, в котором есть ошибка в препроцессоре)
#	define GET_ARGUMENTS_COUNT_IMPLEMENTATION(...) \
		COMBINE(GET_NTH_ARGUMENT, (__VA_ARGS__))

#endif


// вернуть 11-ый параметр макроса
#define GET_NTH_ARGUMENT(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, n, ...) \
	n

// получить количество параметров макроса (не определяет нулевое число параметров)
#define GET_ARGUMENTS_COUNT_WITHOUT_ZERO(...) \
	GET_ARGUMENTS_COUNT_IMPLEMENTATION(__VA_ARGS__, SEQUENCE_N())

// получить количество параметров макроса
#define GET_ARGUMENTS_COUNT(...) \
	GET_ARGUMENTS_COUNT_WITHOUT_ZERO(SEQUENCE_N __VA_ARGS__())


#endif // ARGUMENTS_COUNT_MACRO_H