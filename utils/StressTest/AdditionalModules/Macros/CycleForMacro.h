#ifndef CYCLE_FOR_MACRO_H

#define CYCLE_FOR_MACRO_H


#include "ArithmeticMacros.h"
#include "StringMacros.h"


// выполнить итерацию цикла "for" с номером 0
#define CYCLE_FOR_ITERATION_0(executeIteration, iterationNumber)

// выполнить итерацию цикла "for" с номером 1
#define CYCLE_FOR_ITERATION_1(executeIteration, iterationNumber) \
	executeIteration(iterationNumber)

// выполнить итерацию цикла "for" с номером 2
#define CYCLE_FOR_ITERATION_2(executeIteration, iterationNumber) \
	executeIteration(iterationNumber) \
	CYCLE_FOR_ITERATION_1(executeIteration, ARITHMETIC_INC(iterationNumber))

// выполнить итерацию цикла "for" с номером 3
#define CYCLE_FOR_ITERATION_3(executeIteration, iterationNumber) \
	executeIteration(iterationNumber) \
	CYCLE_FOR_ITERATION_2(executeIteration, ARITHMETIC_INC(iterationNumber))

// выполнить итерацию цикла "for" с номером 4
#define CYCLE_FOR_ITERATION_4(executeIteration, iterationNumber) \
	executeIteration(iterationNumber) \
	CYCLE_FOR_ITERATION_3(executeIteration, ARITHMETIC_INC(iterationNumber))

// выполнить итерацию цикла "for" с номером 5
#define CYCLE_FOR_ITERATION_5(executeIteration, iterationNumber) \
	executeIteration(iterationNumber) \
	CYCLE_FOR_ITERATION_4(executeIteration, ARITHMETIC_INC(iterationNumber))

// выполнить итерацию цикла "for" с номером 6
#define CYCLE_FOR_ITERATION_6(executeIteration, iterationNumber) \
	executeIteration(iterationNumber) \
	CYCLE_FOR_ITERATION_5(executeIteration, ARITHMETIC_INC(iterationNumber))

// выполнить итерацию цикла "for" с номером 7
#define CYCLE_FOR_ITERATION_7(executeIteration, iterationNumber) \
	executeIteration(iterationNumber) \
	CYCLE_FOR_ITERATION_6(executeIteration, ARITHMETIC_INC(iterationNumber))

// выполнить итерацию цикла "for" с номером 8
#define CYCLE_FOR_ITERATION_8(executeIteration, iterationNumber) \
	executeIteration(iterationNumber) \
	CYCLE_FOR_ITERATION_7(executeIteration, ARITHMETIC_INC(iterationNumber))

// выполнить итерацию цикла "for" с номером 9
#define CYCLE_FOR_ITERATION_9(executeIteration, iterationNumber) \
	executeIteration(iterationNumber) \
	CYCLE_FOR_ITERATION_8(executeIteration, ARITHMETIC_INC(iterationNumber))

// выполнить итерацию цикла "for" с номером 10
#define CYCLE_FOR_ITERATION_10(executeIteration, iterationNumber) \
	executeIteration(iterationNumber) \
	CYCLE_FOR_ITERATION_9(executeIteration, ARITHMETIC_INC(iterationNumber))

// выполнить цикл "for"
#define CYCLE_FOR(executeIteration, startWith, iterationsCount) \
	CONCATENATE(CYCLE_FOR_ITERATION_, iterationsCount)(executeIteration, startWith)


// выполнить специальную версию цикла "for", в которой отличаются действия,
// выполняемые во время первой и последующих итераций
#define CYCLE_FOR_SPECIAL(executeIteration0, executeIteration, startWith, iterationsCount) \
	CYCLE_FOR(executeIteration0, startWith, MAKE_BOOLEAN(iterationsCount)) \
	CYCLE_FOR(executeIteration, ARITHMETIC_INC(startWith), ARITHMETIC_DEC(iterationsCount))


#endif // CYCLE_FOR_MACRO_H