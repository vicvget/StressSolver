#ifndef ARITHMETIC_MACROS_H

#define ARITHMETIC_MACROS_H


#include "StringMacros.h"


// увеличить число 0 на единицу
#define ARITHMETIC_INC_0 1

// увеличить число 1 на единицу
#define ARITHMETIC_INC_1 2

// увеличить число 2 на единицу
#define ARITHMETIC_INC_2 3

// увеличить число 3 на единицу
#define ARITHMETIC_INC_3 4

// увеличить число 4 на единицу
#define ARITHMETIC_INC_4 5

// увеличить число 5 на единицу
#define ARITHMETIC_INC_5 6

// увеличить число 6 на единицу
#define ARITHMETIC_INC_6 7

// увеличить число 7 на единицу
#define ARITHMETIC_INC_7 8

// увеличить число 8 на единицу
#define ARITHMETIC_INC_8 9

// увеличить число 9 на единицу
#define ARITHMETIC_INC_9 10

// увеличить число 10 на единицу
#define ARITHMETIC_INC_10 11

// увеличить число на единицу
#define ARITHMETIC_INC(value) \
	CONCATENATE(ARITHMETIC_INC_, value)


// уменьшить число 0 на единицу (выставлено для корректной работы ряда макросов)
#define ARITHMETIC_DEC_0 0

// уменьшить число 1 на единицу
#define ARITHMETIC_DEC_1 0

// уменьшить число 2 на единицу
#define ARITHMETIC_DEC_2 1

// уменьшить число 3 на единицу
#define ARITHMETIC_DEC_3 2

// уменьшить число 4 на единицу
#define ARITHMETIC_DEC_4 3

// уменьшить число 5 на единицу
#define ARITHMETIC_DEC_5 4

// уменьшить число 6 на единицу
#define ARITHMETIC_DEC_6 5

// уменьшить число 7 на единицу
#define ARITHMETIC_DEC_7 6

// уменьшить число 8 на единицу
#define ARITHMETIC_DEC_8 7

// уменьшить число 9 на единицу
#define ARITHMETIC_DEC_9 8

// уменьшить число 10 на единицу
#define ARITHMETIC_DEC_10 9

// уменьшить число на единицу
#define ARITHMETIC_DEC(value) \
	CONCATENATE(ARITHMETIC_DEC_, value)


#endif // ARITHMETIC_MACROS_H