#ifndef SMART_POINTERS_H

#define SMART_POINTERS_H


#include <memory>
#include <functional>


// ������������ ���� ��� ���������������� ����, �������������� �������� ��� �� �������� ����� ����������
namespace SmartPointersImplementation
{

	// ��� ������������� ����������������� deleter'�
	template
		<
			class Type // ����� �������, ������� ����� ������� std::unique_ptr
		>
	using DeleterType = std::function<void(Type*)>;


	/**
	* ������� ���������������� deleter ���� DeleterType
	* @param deleter - ���������������� deleter, �������-���� ������ Type
	* @return ��������� ���������������� deleter ���� DeleterType
	*/
	template
		<
			class Type, // ����� �������, ������� ����� ������� std::unique_ptr
			typename DeleterPrototype, // �������� �������-����� ������ Type, ������� ��������� � ���� deleter'�
			typename Deleter = DeleterType<Type> // ��� ����������������� deleter'�
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


// ��� ������������� ������ ��������� � ���������������� deleter'��
template
	<
		class Type // ����� �������, ������� ����� ������� �4���� ���������
	>
using UniquePointerWithDeleter =
	std::unique_ptr<Type, SmartPointersImplementation::DeleterType<Type>>;


/**
* ������� ��������� std::unique_ptr � ���������������� deleter'��
* @param deleter - ���������������� deleter
* @param arguments - ���������, ������������ ������������ ������������ �������
* @return ��������� ��������� std::unique_ptr � ���������������� deleter'��
*/
template
	<
		typename Type, // ��� �������, ������� ����� ������� std::unique_ptr
		typename Deleter, // ��� ����������������� deleter'�
		typename... Arguments // ���� ����������, ������������ ������������ ������������ �������
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
* ������� ��������� std::unique_ptr � ���������������� deleter'��,
* ������� �������� ��������-������ ������ ��������� �������
* @param deleter - ���������������� deleter, �������-���� ������ Type
* @param arguments - ���������, ������������ ������������ ������������ �������
* @return ��������� ��������� std::unique_ptr � ���������������� deleter'��
*/
template
	<
		class Type, // ����� �������, ������� ����� ������� std::unique_ptr
		typename DeleterPrototype, // �������� �������-����� ������ Type, ������� ��������� � ���� deleter'�
		typename... Arguments // ���� ����������, ������������ ������������ ������������ �������
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
* ������� ��������� std::unique_ptr � ���������������� deleter'��,
* ������� �������� ��������-������ ������ ��������� �������
* @param object - ������, ������� ������ ������� ����������� ��������� std::unique_ptr
* @param deleter - ���������������� deleter, �������-���� ������ Type
* @return ��������� ��������� std::unique_ptr � ���������������� deleter'��
*/
template
	<
		class Type, // ����� �������, ������� ����� ������� std::unique_ptr
		typename DeleterPrototype // �������� �������-����� ������ Type, ������� ��������� � ���� deleter'�
	>
UniquePointerWithDeleter<Type> MakeUniqueWithMemberDeleter
	(
		Type* object,
		DeleterPrototype Type::* deleter
	)
{
	return {object, SmartPointersImplementation::MakeUserDeleter(deleter)};
}


// ������ ������� �������� ������� ������ std::unique_ptr<Type>,
// � ������� � ���� deleter'� ��������� �������-���� Type::Function
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