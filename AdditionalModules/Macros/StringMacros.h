#ifndef STRING_MACROS_H

#define STRING_MACROS_H


// запятая
#define COMMA \
	,


// реализация макроса объединения двух аргументов
#define CONCATENATE_IMPLEMENTATION(a, b) \
	a ## b

// объединить две аргумента
#define CONCATENATE(a, b) \
	CONCATENATE_IMPLEMENTATION(a, b)


// превратить аргумент макроса в строковый литерал (уровень вложенности 1)
#define STRINGIZE1(x) \
	#x

// превратить аргумент макроса в строковый литерал (уровень вложенности 2)
#define STRINGIZE2(x) \
	STRINGIZE1(x)

// превратить аргумент макроса в строковый литерал (уровень вложенности 3)
#define STRINGIZE3(x) \
	STRINGIZE2(x)

// превратить аргумент макроса в строковый литерал (уровень вложенности 4)
#define STRINGIZE4(x) \
	STRINGIZE3(x)

// превратить аргумент макроса в строковый литерал (уровень вложенности 5)
#define STRINGIZE5(x) \
	STRINGIZE4(x)

// превратить аргумент макроса в строковый литерал (уровень вложенности 6)
#define STRINGIZE6(x) \
	STRINGIZE5(x)

// превратить аргумент макроса в строковый литерал (уровень вложенности 7)
#define STRINGIZE7(x) \
	STRINGIZE6(x)

// превратить аргумент макроса в строковый литерал (уровень вложенности 8)
#define STRINGIZE8(x) \
	STRINGIZE7(x)

// превратить аргумент макроса в строковый литерал (уровень вложенности 9)
#define STRINGIZE9(x) \
	STRINGIZE8(x)

// превратить аргумент макроса в строковый литерал (уровень вложенности 10)
#define STRINGIZE10(x) \
	STRINGIZE9(x)

// превратить аргумент макроса в строковый литерал (используется 11 уровней вложенности)
#define STRINGIZE(x) \
	STRINGIZE10(x)


#endif // STRING_MACROS_H