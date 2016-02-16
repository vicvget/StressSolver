#ifndef CONTROL_MACROS_H

#define CONTROL_MACROS_H


#include "AdditionalMacros.h"
#include "LogicalMacros.h"
#include "StringMacros.h"


// выполнить действие, соответствующее ветке "ложь"
#define STATEMENT_IF_0(trueAction)

// выполнить действие, соответствующее ветке "истина"
#define STATEMENT_IF_1(trueAction) \
	trueAction

// реализация операции "выполнить действие по условию"
#define STATEMENT_IF_IMPLEMENTATION(condition, trueAction) \
	CONCATENATE(STATETEMENT_IF_, condition)(trueAction)

// операция "выполнить действие по условию"
#define STATEMENT_IF(condition, trueAction) \
	STATETEMENT_IF_IMPLEMENTATION(MAKE_BOOLEAN(condition), trueAction)


// выполнить действие, соответствующее ветке "ложь"
#define STATEMENT_IF_ELSE_0(trueAction, falseAction) \
	falseAction

// выполнить действие, соответствующее ветке "истина"
#define STATEMENT_IF_ELSE_1(trueAction, falseAction) \
	trueAction

// реализация операции "ветвление"
#define STATEMENT_IF_ELSE_IMPLEMENTATION(condition, trueAction, falseAction) \
	CONCATENATE(STATEMENT_IF_ELSE_, condition)(trueAction, falseAction)

// операция "ветвление"
#define STATEMENT_IF_ELSE(condition, trueAction, falseAction) \
	STATEMENT_IF_ELSE_IMPLEMENTATION(MAKE_BOOLEAN(condition), trueAction, falseAction)


// не выполнять действие (условие ложно)
#define EXECUTE_IF_0(action, ...)


#if !defined(MS_VC_COMPILER)

	// выполнить действие (условие истинно)
#	define EXECUTE_IF_1(action, ...) \
		action(__VA_ARGS__)

// реализация операции выполнения действия по условию
#	define EXECUTE_IF_IMPLEMENTATION(condition, ...) \
		CONCATENATE(EXECUTE_IF_, condition)(__VA_ARGS__)

// выполнить действие по условию
#	define EXECUTE_IF(condition, ...) \
		EXECUTE_IF_IMPLEMENTATION(MAKE_BOOLEAN(condition), __VA_ARGS__)

#else

	// выполнить действие (условие истинно)
#	define EXECUTE_IF_1(action, ...) \
		COMBINE(action, (__VA_ARGS__))

// реализация операции выполнения действия по условию
#	define EXECUTE_IF_IMPLEMETATION(condition, action, ...) \
		CONCATENATE(EXECUTE_IF_, condition)(action, __VA_ARGS__)

// выполнить действие по условию
#	define EXECUTE_IF(condition, action, ...) \
		EXECUTE_IF_IMPLEMETATION(MAKE_BOOLEAN(condition), action, __VA_ARGS__)

#endif


#endif // CONTROL_MACROS_H