#ifndef EXUMERATIONS_H

#define EXUMERATIONS_H


#include <type_traits>


// определение, является ли тип Type строго типизированным перечислением
// (строго типизировнное перечисление - это перечисление, элементы которого неявно не приводятся к int)
template
	<
		typename Type // проверяемый тип
	>
using IsScopedEnumeration = std::integral_constant
	<
		bool,
		std::is_enum<Type>::value && !std::is_convertible<Type, int>::value
	>;


// явное приведение значения переменной типа "перечисление" к underlying типу этого перечисления (compile-time)
#define TO_INTEGRAL_TYPE(enumElement) \
	static_cast<std::underlying_type_t<decltype(enumElement)>>(enumElement)


/**
* Явное приведение значения переменной типа "перечисление" к underlying типу этого перечисления (run-time)
* @param e - переменная типа "перечисление" (Enumeration)
* @return значение underlying типа для перечисления Enumeration, соответствующее значению переменной e
*/
template
	<
		typename Enumeration // тип перечисления
	>
// запрещаем использование этой функции для типов, не являющихся перечислениями
std::enable_if_t
	<
		IsScopedEnumeration<Enumeration>::value,
		std::underlying_type_t<Enumeration>
	>
ToIntegralType
	(
		Enumeration e
	)
{
	return TO_INTEGRAL_TYPE(e);
}


/**
* Оператор "побитовое 'или'" для переменных типа "перечисление"
* @param e1 - переменная типа "перечисление" (Enumeration)
* @param e2 - переменная типа "перечисление" (Enumeration)
* @return результат "побитового 'или'" для значений underlying типа,
* соответствующих значениям переменных e1 и e2 и преобразованный к типу Enumeration
*/
template
	<
		typename Enumeration // тип перечисления
	>
// запрещаем использование этого оператора для типов, не являющихся перечислениями
std::enable_if_t
	<
		IsScopedEnumeration<Enumeration>::value,
		Enumeration
	>
operator |
	(
		Enumeration e1,
		Enumeration e2
	)
{
	return static_cast<Enumeration>(ToIntegralType(e1) | ToIntegralType(e2));
}


/**
* Оператор "побитовое 'или' с присвоением" для переменных типа "перечисление"
* @param e - переменная типа "перечисление", которой будет присвоен результат операции (Enumeration)
* @param otherE - переменная типа "перечисление" (Enumeration)
* @return ссылка на результат "побитового 'или'" для значений underlying типа,
* соответствующих значениям переменных e и otherE и преобразованный к типу Enumeration,
* сохраненный в переменной e
*/
template
	<
		typename Enumeration
	>
// запрещаем использование этого оператора для типов, не являющихся перечислениями
std::enable_if_t
	<
		IsScopedEnumeration<Enumeration>::value,
		Enumeration&
	>
operator |=
	(
		Enumeration& e,
		Enumeration otherE
	)
{
	e = e | otherE;

	return e;
}


/**
* Оператор "побитовое 'и'" для переменных типа "перечисление"
* @param e1 - переменная типа "перечисление" (Enumeration)
* @param e2 - переменная типа "перечисление" (Enumeration)
* @return результат "побитового 'и'" для значений underlying типа,
* соответствующих значениям переменных e1 и e2 и преобразованный к типу Enumeration
*/
template
	<
		typename Enumeration // тип перечисления
	>
// запрещаем использование этого оператора для типов, не являющихся перечислениями
std::enable_if_t
	<
		IsScopedEnumeration<Enumeration>::value,
		Enumeration
	>
operator &
	(
		Enumeration e1,
		Enumeration e2
	)
{
	return static_cast<Enumeration>(ToIntegralType(e1) & ToIntegralType(e2));
}


/**
* Оператор "побитовое 'и' с присвоением" для переменных типа "перечисление"
* @param e - переменная типа "перечисление", которой будет присвоен результат операции (Enumeration)
* @param otherE - переменная типа "перечисление" (Enumeration)
* @return ссылка на результат "побитового 'и'" для значений underlying типа,
* соответствующих значениям переменных e и otherE и преобразованный к типу Enumeration,
* сохраненный в переменной e
*/
template
	<
		typename Enumeration
	>
// запрещаем использование этого оператора для типов, не являющихся перечислениями
std::enable_if_t
	<
		IsScopedEnumeration<Enumeration>::value,
		Enumeration&
	>
operator &=
	(
		Enumeration& e,
		Enumeration otherE
	)
{
	e = e & otherE;

	return e;
}


/**
* Оператор "сдвиг влево" для переменных типа "перечисление"
* @param e - переменная типа "перечисление" (Enumeration)
* @param shift - переменная underlying типа для типа "перечисление" (Enumeration)?
* равная числу битов, на которое производится сдвиг
* @return результат "сдвига влево" для значения underlying типа,
* соответствующего значению переменной e и преобразованный к типу Enumeration
*/
template
	<
		typename Enumeration // тип перечисления
	>
// запрещаем использование этого оператора для типов, не являющихся перечислениями
std::enable_if_t
	<
		IsScopedEnumeration<Enumeration>::value,
		Enumeration
	>
operator <<
	(
		Enumeration e,
		std::underlying_type_t<Enumeration> shift
	)
{
	return static_cast<Enumeration>(ToIntegralType(e) << shift);
}


/**
* Оператор "сдвиг влево с присвоением" для переменной типа "перечисление"
* @param e - переменная типа "перечисление", которой будет присвоен результат операции (Enumeration)
* @param shift - переменная underlying типа для типа "перечисление" (Enumeration),
* равная числу битов, на которое производится сдвиг
* @return ссылка на результат "сдвига влево" для значений underlying типа,
* соответствующего значению переменной e и преобразованный к типу Enumeration,
* сохраненный в переменной e
*/
template
	<
		typename Enumeration
	>
// запрещаем использование этого оператора для типов, не являющихся перечислениями
std::enable_if_t
	<
		IsScopedEnumeration<Enumeration>::value,
		Enumeration&
	>
operator <<=
	(
		Enumeration& e,
		std::underlying_type_t<Enumeration> shift
	)
{
	e = e << shift;

	return e;
}


#endif // EXUMERATIONS_H
