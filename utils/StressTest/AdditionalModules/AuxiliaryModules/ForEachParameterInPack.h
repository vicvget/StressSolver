#ifndef FOR_EACH_PARAMETER_IN_PACK_H

#define FOR_EACH_PARAMETER_IN_PACK_H


/**
* Вспомогательная структура, необходимая для выполнения некоторой операции
* над variadic template parameter pack
*/
struct Do
{

	// Конструкторы и деструктор

	/**
	* Конструктор (принимает любое число параметров)
	*/
	Do(...)
	{
		// ничего не надо делать
	}

};


/**
* Функция, которая возвращает значение, переданное на вход,
* и одновременно аккумулирует некоторый тип
* @param value входное значение, которое будет передано на выход
* @return входное значение
*/
template
	<
		typename TypeToAccumulate, // тип, который необходимо аккумулировать
		typename Type // тип значения, переданного на вход
	>
auto Accumulate
	(
		Type&& value
	)
	->	decltype(std::forward<Type>(value))
{
	return std::forward<Type>(value);
}


#endif // FOR_EACH_PARAMETER_IN_PACK_H