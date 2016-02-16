#ifndef LIST_MACROS_H

#define LIST_MACROS_H


#include "ArgumentsCountMacro.h"
#include "ControlMacros.h"


// откусить голову списка, состоящего по крайней мере из одного элемента
#define HEAD_OF_LIST(head, ...) \
	head


// откусить голову списка
#define GET_HEAD_OF_LIST(...) \
	EXECUTE_IF(GET_ARGUMENTS_COUNT(__VA_ARGS__), HEAD_OF_LIST, __VA_ARGS__)


// откусить хвост списка, состоящего по крайней мере из одного элемента
#define TAIL_OF_LIST(head, ...) \
	__VA_ARGS__


// откусить хвост списка
#define GET_TAIL_OF_LIST(...) \
	EXECUTE_IF(GET_ARGUMENTS_COUNT(__VA_ARGS__), TAIL_OF_LIST, __VA_ARGS__)


#endif // LIST_MACROS_H