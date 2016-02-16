#ifndef FOR_EACH_PARAMETER_IN_PACK_H

#define FOR_EACH_PARAMETER_IN_PACK_H


/**
* ��������������� ���������, ����������� ��� ���������� ��������� ��������
* ��� variadic template parameter pack
*/
struct Do
{

	// ������������ � ����������

	/**
	* ����������� (��������� ����� ����� ����������)
	*/
	Do(...)
	{
		// ������ �� ���� ������
	}

};


/**
* �������, ������� ���������� ��������, ���������� �� ����,
* � ������������ ������������ ��������� ���
* @param value ������� ��������, ������� ����� �������� �� �����
* @return ������� ��������
*/
template
	<
		typename TypeToAccumulate, // ���, ������� ���������� ��������������
		typename Type // ��� ��������, ����������� �� ����
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