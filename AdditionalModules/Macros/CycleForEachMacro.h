#ifndef CYCLE_FOR_EACH_MACRO_H

#define CYCLE_FOR_EACH_MACRO_H


#include "ArgumentsCountMacro.h"
#include "ControlMacros.h"
#include "ListMacros.h"
#include "StringMacros.h"


#if defined(MS_VC_COMPILER)
	// реализация цикла "foreach" для компилятора Microsoft Visual C++

	// выполнить следующую итерацию цикла "foreach"
#	define CYCLE_FOR_EACH_NEXT_ITERATION(iterationNumber, executeIteration, ...) \
		EXECUTE_IF \
			( \
				GET_ARGUMENTS_COUNT(__VA_ARGS__), \
				CONCATENATE(CYCLE_FOR_EACH_ITERATION_, iterationNumber), \
				executeIteration, \
				__VA_ARGS__ \
			)

// выполнить итерацию цикла "foreach"
#	define CYCLE_FOR_EACH_ITERATION(iterationNumber, executeIteration, param, ...) \
		executeIteration(param, iterationNumber) \
		CYCLE_FOR_EACH_NEXT_ITERATION(ARITHMETIC_INC(iterationNumber), executeIteration, __VA_ARGS__)

	// выполнить итерацию цикла "foreach" с номером 0
#	define CYCLE_FOR_EACH_ITERATION_0(executeIteration, param, ...) \
		CYCLE_FOR_EACH_ITERATION(0, executeIteration, param, __VA_ARGS__)

	// выполнить итерацию цикла "foreach" с номером 1
#	define CYCLE_FOR_EACH_ITERATION_1(executeIteration, param, ...) \
		CYCLE_FOR_EACH_ITERATION(1, executeIteration, param, __VA_ARGS__)

	// выполнить итерацию цикла "foreach" с номером 2
#	define CYCLE_FOR_EACH_ITERATION_2(executeIteration, param, ...) \
		CYCLE_FOR_EACH_ITERATION(2, executeIteration, param, __VA_ARGS__)

	// выполнить итерацию цикла "foreach" с номером 3
#	define CYCLE_FOR_EACH_ITERATION_3(executeIteration, param, ...) \
		CYCLE_FOR_EACH_ITERATION(3, executeIteration, param, __VA_ARGS__)

	// выполнить итерацию цикла "foreach" с номером 4
#	define CYCLE_FOR_EACH_ITERATION_4(executeIteration, param, ...) \
		CYCLE_FOR_EACH_ITERATION(4, executeIteration, param, __VA_ARGS__)

	// выполнить итерацию цикла "foreach" с номером 5
#	define CYCLE_FOR_EACH_ITERATION_5(executeIteration, param, ...) \
		CYCLE_FOR_EACH_ITERATION(5, executeIteration, param, __VA_ARGS__)

	// выполнить итерацию цикла "foreach" с номером 6
#	define CYCLE_FOR_EACH_ITERATION_6(executeIteration, param, ...) \
		CYCLE_FOR_EACH_ITERATION(6, executeIteration, param, __VA_ARGS__)

	// выполнить итерацию цикла "foreach" с номером 7
#	define CYCLE_FOR_EACH_ITERATION_7(executeIteration, param, ...) \
		CYCLE_FOR_EACH_ITERATION(7, executeIteration, param, __VA_ARGS__)

	// выполнить итерацию цикла "foreach" с номером 8
#	define CYCLE_FOR_EACH_ITERATION_8(executeIteration, param, ...) \
		CYCLE_FOR_EACH_ITERATION(8, executeIteration, param, __VA_ARGS__)

	// выполнить итерацию цикла "foreach" с номером 9
#	define CYCLE_FOR_EACH_ITERATION_9(executeIteration, param, ...) \
		CYCLE_FOR_EACH_ITERATION(9, executeIteration, param, __VA_ARGS__)

	// выполнить итерацию цикла "foreach" с номером 10
#	define CYCLE_FOR_EACH_ITERATION_10(executeIteration, param, ...) \
		executeIteration(param)

	// выполнить цикл "foreach"
#	define CYCLE_FOR_EACH(executeIteration, ...) \
		CYCLE_FOR_EACH_NEXT_ITERATION(0, executeIteration, __VA_ARGS__)


	// выполнить специальную версию цикла "foreach", в которой отличаются действия,
	// выполняемые во время первой и последующих итераций
#	define CYCLE_FOR_EACH_SPECIAL(executeIteration0, executeIteration, ...) \
		CYCLE_FOR_EACH(executeIteration0, GET_HEAD_OF_LIST(__VA_ARGS__)) \
		CYCLE_FOR_EACH_NEXT_ITERATION(1, executeIteration, GET_TAIL_OF_LIST(__VA_ARGS__))


#elif defined(GNU_GCC_COMPILER)
	// реализация цикла "foreach" для компилятора GNU GCC C/C++

	// выполнить итерацию цикла "foreach" с номером 0
#	define CYCLE_FOR_EACH_ITERATION_0(iterationNumber, executeIteration)

	// выполнить итерацию цикла "foreach" с номером 1
#	define CYCLE_FOR_EACH_ITERATION_1(iterationNumber, executeIteration, param) \
		executeIteration(param, iterationNumber) \
		CYCLE_FOR_EACH_ITERATION_0(ARITHMETIC_INC(iterationNumber), executeIteration)

	// выполнить итерацию цикла "foreach" с номером 2
#	define CYCLE_FOR_EACH_ITERATION_2(iterationNumber, executeIteration, param, ...) \
		executeIteration(param, iterationNumber) \
		CYCLE_FOR_EACH_ITERATION_1(ARITHMETIC_INC(iterationNumber), executeIteration, __VA_ARGS__)

	// выполнить итерацию цикла "foreach" с номером 3
#	define CYCLE_FOR_EACH_ITERATION_3(iterationNumber, executeIteration, param, ...) \
		executeIteration(param, iterationNumber) \
		CYCLE_FOR_EACH_ITERATION_2(ARITHMETIC_INC(iterationNumber), executeIteration, __VA_ARGS__)

	// выполнить итерацию цикла "foreach" с номером 4
#	define CYCLE_FOR_EACH_ITERATION_4(iterationNumber, executeIteration, param, ...) \
		executeIteration(param, iterationNumber) \
		CYCLE_FOR_EACH_ITERATION_3(ARITHMETIC_INC(iterationNumber), executeIteration, __VA_ARGS__)

	// выполнить итерацию цикла "foreach" с номером 5
#	define CYCLE_FOR_EACH_ITERATION_5(iterationNumber, executeIteration, param, ...) \
		executeIteration(param, iterationNumber) \
		CYCLE_FOR_EACH_ITERATION_4(ARITHMETIC_INC(iterationNumber), executeIteration, __VA_ARGS__)

	// выполнить итерацию цикла "foreach" с номером 6
#	define CYCLE_FOR_EACH_ITERATION_6(iterationNumber, executeIteration, param, ...) \
		executeIteration(param, iterationNumber) \
		CYCLE_FOR_EACH_ITERATION_5(ARITHMETIC_INC(iterationNumber), executeIteration, __VA_ARGS__)

	// выполнить итерацию цикла "foreach" с номером 7
#	define CYCLE_FOR_EACH_ITERATION_7(iterationNumber, executeIteration, param, ...) \
		executeIteration(param, iterationNumber) \
		CYCLE_FOR_EACH_ITERATION_6(ARITHMETIC_INC(iterationNumber), executeIteration, __VA_ARGS__)

	// выполнить итерацию цикла "foreach" с номером 8
#	define CYCLE_FOR_EACH_ITERATION_8(iterationNumber, executeIteration, param, ...) \
		executeIteration(param, iterationNumber) \
		CYCLE_FOR_EACH_ITERATION_7(ARITHMETIC_INC(iterationNumber), executeIteration, __VA_ARGS__)

	// выполнить итерацию цикла "foreach" с номером 9
#	define CYCLE_FOR_EACH_ITERATION_9(iterationNumber, executeIteration, param, ...) \
		executeIteration(param, iterationNumber) \
		CYCLE_FOR_EACH_ITERATION_8(ARITHMETIC_INC(iterationNumber), executeIteration, __VA_ARGS__)

	// выполнить итерацию цикла "foreach" с номером 10
#	define CYCLE_FOR_EACH_ITERATION_10(iterationNumber, executeIteration, param, ...) \
		executeIteration(param, iterationNumber) \
		CYCLE_FOR_EACH_ITERATION_9(ARITHMETIC_INC(iterationNumber), executeIteration, __VA_ARGS__)


	// реализация цикла "foreach"
#	define CYCLE_FOR_EACH_IMPLEMENTATION(startWith, ...) \
		CONCATENATE(CYCLE_FOR_EACH_ITERATION_, ARITHMETIC_DEC(GET_ARGUMENTS_COUNT(__VA_ARGS__))) \
			(startWith, __VA_ARGS__)

	// выполнить цикл "foreach"
#	define CYCLE_FOR_EACH(...) \
		CYCLE_FOR_EACH_IMPLEMENTATION(0, __VA_ARGS__)


#	define CYCLE_FOR_EACH_SPECIAL_HEAD(executeIteration0, ...) \
		EXECUTE_IF(GET_ARGUMENTS_COUNT(__VA_ARGS__), CYCLE_FOR_EACH, executeIteration0, HEAD_OF_LIST(__VA_ARGS__))

#	define CYCLE_FOR_EACH_SPECIAL_TAIL(executeIteration, ...) \
		EXECUTE_IF \
			( \
				ARITHMETIC_DEC(GET_ARGUMENTS_COUNT(__VA_ARGS__)), \
				CYCLE_FOR_EACH_IMPLEMENTATION, 1, executeIteration, TAIL_OF_LIST(__VA_ARGS__) \
			)

	// выполнить специальную версию цикла "foreach", в которой отличаются действия,
	// выполняемые во время первой и последующих итераций
#	define CYCLE_FOR_EACH_SPECIAL(executeIteration0, executeIteration, ...) \
		CYCLE_FOR_EACH_SPECIAL_HEAD(executeIteration0, __VA_ARGS__) \
		CYCLE_FOR_EACH_SPECIAL_TAIL(executeIteration, __VA_ARGS__)

#else

	// Реализация цикла "foreach" доступна только для компиляторов Microsoft Visual C++
	// и GNU GCC C/C++. Для других компиляторов требуется проверка и доработка кода.
#	error "foreach" cycle implementation is available only for Microsoft Visual C++ \
and GNU GCC C/C++. Manual check and code modification are needed for other compilers.

#endif


#endif // CYCLE_FOR_EACH_MACRO_H