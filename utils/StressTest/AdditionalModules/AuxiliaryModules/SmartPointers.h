#ifndef SMART_POINTERS_H

#define SMART_POINTERS_H


#include <memory>
#include <functional>


// пространство имен для вспомогательного кода, обслуживающего основной код по созданию умных указателей
namespace SmartPointersImplementation
{

	// тип используемого пользовательского deleter'а
	template
		<
			class Type // класс объекта, который будет хранить std::unique_ptr
		>
	using DeleterType = std::function<void(Type*)>;


	/**
	* Создать пользовательский deleter типа DeleterType
	* @param deleter - пользовательский deleter, функция-член класса Type
	* @return созданный пользовательский deleter типа DeleterType
	*/
	template
		<
			class Type, // класс объекта, который будет хранить std::unique_ptr
			typename DeleterPrototype, // прототип функции-члена класса Type, которая выступает в роли deleter'а
			typename Deleter = DeleterType<Type> // тип пользовательского deleter'а
		>
	Deleter MakeUserDeleter
		(
			DeleterPrototype Type::* deleter
		)
	{
		return
			[deleter]
				(
					Type* object
				)
			{
				(object->*deleter)();
				delete object;
			};
	}

}


// тип используемого умного указателя с пользовательским deleter'ом
template
	<
		class Type // класс объекта, который будет хранить у4мный указатель
	>
using UniquePointerWithDeleter =
	std::unique_ptr<Type, SmartPointersImplementation::DeleterType<Type>>;


/**
* Создать экземпляр std::unique_ptr с пользовательским deleter'ом
* @param deleter - пользовательский deleter
* @param arguments - аргументы, передаваемые конструктору создаваемого объекта
* @return созданный экземпляр std::unique_ptr с пользовательским deleter'ом
*/
template
	<
		typename Type, // тип объекта, который будет хранить std::unique_ptr
		typename Deleter, // тип пользовательского deleter'а
		typename... Arguments // типы аргументов, передаваемых конструктору создаваемого объекта
	>
std::unique_ptr<Type, Deleter> MakeUniqueWithDeleter
	(
		Deleter&& deleter,
		Arguments&&... arguments
	)
{
	return {new Type(std::forward<Arguments>(arguments)...), std::forward<Deleter>(deleter)};
}

/**
* Создать экземпляр std::unique_ptr с пользовательским deleter'ом,
* который является функцией-членом класса хранимого объекта
* @param deleter - пользовательский deleter, функция-член класса Type
* @param arguments - аргументы, передаваемые конструктору создаваемого объекта
* @return созданный экземпляр std::unique_ptr с пользовательским deleter'ом
*/
template
	<
		class Type, // класс объекта, который будет хранить std::unique_ptr
		typename DeleterPrototype, // прототип функции-члена класса Type, которая выступает в роли deleter'а
		typename... Arguments // типы аргументов, передаваемых конструктору создаваемого объекта
	>
UniquePointerWithDeleter<Type> MakeUniqueWithMemberDeleter
	(
		DeleterPrototype Type::* deleter,
		Arguments&&... arguments
	)
{
	return MakeUniqueWithDeleter<Type>
		(
			SmartPointersImplementation::MakeUserDeleter(deleter),
			std::forward<Arguments>(arguments)...
		);
}

/**
* Создать экземпляр std::unique_ptr с пользовательским deleter'ом,
* который является функцией-членом класса хранимого объекта
* @param object - объект, который должен хранить создаваемый экземпляр std::unique_ptr
* @param deleter - пользовательский deleter, функция-член класса Type
* @return созданный экземпляр std::unique_ptr с пользовательским deleter'ом
*/
template
	<
		class Type, // класс объекта, который будет хранить std::unique_ptr
		typename DeleterPrototype // прототип функции-члена класса Type, которая выступает в роли deleter'а
	>
UniquePointerWithDeleter<Type> MakeUniqueWithMemberDeleter
	(
		Type* object,
		DeleterPrototype Type::* deleter
	)
{
	return {object, SmartPointersImplementation::MakeUserDeleter(deleter)};
}


// Макрос функции создания объекта класса std::unique_ptr<Type>,
// в котором в роли deleter'а выступает функция-член Type::Function
#define DEFINE_CREATE_OBJECT_MEMBER_FUNCTION(Type, Function) \
 \
template \
	< \
		typename... Arguments \
	> \
static \
auto Create \
	( \
		Arguments&&... arguments \
	) \
	-> decltype(MakeUniqueWithMemberDeleter<Type>(&Type::Function, std::forward<Arguments>(arguments)...)) \
{ \
	return MakeUniqueWithMemberDeleter<Type>(&Type::Function, std::forward<Arguments>(arguments)...); \
}


#endif // SMART_POINTERS_H