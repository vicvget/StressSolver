#ifndef ARGUMENTS_LIST_MACROS_H

#define ARGUMENTS_LIST_MACROS_H


#include "ArgumentsCountMacro.h"
#include "StringMacros.h"
#include "CycleForMacro.h"
#include "CycleForEachMacro.h"


// TODO: отделить макросы, отвечающие за формирование списка полей класса
// в отдельный заголовочный файл (или переименовать соответствующим образом этот)


// сформировать наименование переменной
#define VARIABLE_NAME(number) \
	CONCATENATE(variable, number)


// сформировать наименование первого аргумента функции
#define ARGUMENT_NAME(number) \
	VARIABLE_NAME(number)


// сформировать наименование одного из последующих аргументов функции (с запятой в начале)
#define NEXT_ARGUMENT_NAME(number) \
	ADD_COMMA_SINGLE(ARGUMENT_NAME(number))


// сформировать объявление первого аргумента функции
#define ARGUMENT_DECLARATION(type, number) \
	type ARGUMENT_NAME(number)


// сформировать объявление одного из последующих аргументов функции (с запятой в начале)
#define NEXT_ARGUMENT_DECLARATION(type, number) \
	ADD_COMMA_SINGLE(ARGUMENT_DECLARATION(type, number))


// сформировать наименование поля класса
#define FIELD_DECLARATION(type, number) \
	type VARIABLE_NAME(number);


// сформировать список аргументов функции при ее вызове для заданного числа параметров
#define CREATE_ARGUMENTS_LIST_FOR_NUMBER(number) \
	CYCLE_FOR_SPECIAL(ARGUMENT_NAME, NEXT_ARGUMENT_NAME, 0, number)


// сформировать список аргументов функции при ее вызове для заданного списка типов
#define CREATE_ARGUMENTS_LIST(...) \
	CREATE_ARGUMENTS_LIST_FOR_NUMBER(GET_ARGUMENTS_COUNT(__VA_ARGS__))


// сформировать список параметров функции при ее объявлении для заданного списка типов
#define CREATE_PARAMETERS_LIST(...) \
	CYCLE_FOR_EACH_SPECIAL(ARGUMENT_DECLARATION, NEXT_ARGUMENT_DECLARATION, __VA_ARGS__)


// сформировать список полей класса для заданного списка типов
#define CREATE_FIELD_LIST(...) \
	CYCLE_FOR_EACH(FIELD_DECLARATION, __VA_ARGS__)


#endif // ARGUMENTS_LIST_MACROS_H