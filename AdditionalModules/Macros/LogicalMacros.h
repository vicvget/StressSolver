#ifndef LOGICAL_MACROS_H

#define LOGICAL_MACROS_H


#include "StringMacros.h"


// превратить число 0 в значение логического типа (0/1)
#define MAKE_BOOLEAN_0 0

// превратить число 1 в значение логического типа (0/1)
#define MAKE_BOOLEAN_1 1

// превратить число 2 в значение логического типа (0/1)
#define MAKE_BOOLEAN_2 1

// превратить число 3 в значение логического типа (0/1)
#define MAKE_BOOLEAN_3 1

// превратить число 4 в значение логического типа (0/1)
#define MAKE_BOOLEAN_4 1

// превратить число 5 в значение логического типа (0/1)
#define MAKE_BOOLEAN_5 1

// превратить число 6 в значение логического типа (0/1)
#define MAKE_BOOLEAN_6 1

// превратить число 7 в значение логического типа (0/1)
#define MAKE_BOOLEAN_7 1

// превратить число 8 в значение логического типа (0/1)
#define MAKE_BOOLEAN_8 1

// превратить число 9 в значение логического типа (0/1)
#define MAKE_BOOLEAN_9 1

// превратить число 10 в значение логического типа (0/1)
#define MAKE_BOOLEAN_10 1

// превратить число в значение логического типа (0/1)
#define MAKE_BOOLEAN(value) \
	CONCATENATE(MAKE_BOOLEAN_, value)


// логическое отрицание значения 0
#define LOGICAL_NOT_0 1

// логическое отрицание значения 1
#define LOGICAL_NOT_1 0

// реализация операции логического отрицания
#define LOGICAL_NOT_IMPLEMETATION(value) \
	CONCATENATE(LOGICAL_NOT_, value)

// логическое отрицание
#define LOGICAL_NOT(value) \
	LOGICAL_NOT_IMPLEMETATION(MAKE_BOOLEAN(value))


#endif // LOGICAL_MACROS_H