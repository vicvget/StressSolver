#ifndef INTEGER_SEQUENCE_H

#define INTEGER_SEQUENCE_H


/**
* ��������� �����, ����������� �������� �������� ������������� ���
* � ������������������ ������������� ��������
*/
template
	<
		typename Integer, // ������������� ���
		Integer... Numbers // ������������������ �������� ����� ����
	>
struct IntegerSequence
{
};


/**
* ������������������ �������� (������� ����� ��� size_t)
*/
template
	<
		size_t... Numbers // ������������������ ��������
	>
using IndexSequence = IntegerSequence<size_t, Numbers...>;


/**
* ���������� ���������� ������������ ������������� ������������������
*/
template
	<
		typename Integer, // ������������� ���
		Integer CurrentNumber, // ������� ����� �����, ����������� �������������������
		Integer EndNumber, // �������� ����� �����, ������� ����� ���� ���������� �������������������
		Integer... Numbers // ��������� ������������������ �������� ����� ����
	>
struct IntegerSequenceGenerator_
	:
		IntegerSequenceGenerator_<Integer, CurrentNumber + 1, EndNumber, Numbers..., CurrentNumber + 1>
{
};


/**
* ���������� ���������� ������������ ������������� ������������������
* (������������� ��������� ���������)
*/
template
	<
		typename Integer, // ������������� ���
		Integer EndNumber, // �������� ����� �����, ������� ����� ���� ���������� �������������������
		Integer... Numbers // ��������� ������������������ �������� ����� ����
	>
struct IntegerSequenceGenerator_<Integer, EndNumber, EndNumber, Numbers...>
{

	// ��� ������������� ������������������, ��������������� �����������
	using type = IntegerSequence<Integer, Numbers...>;

};


/**
* ��������� ������������ ������������� ������������������
*/
template
	<
		typename Integer, // ������������� ���
		Integer BeginNumber, // ����� �����, � �������� ������ �������� ��������� ������������������
		Integer EndNumber // �������� ����� �����, ������� ����� ���� ���������� �������������������
	>
struct IntegerSequenceGenerator
	:
		IntegerSequenceGenerator_<Integer, BeginNumber, EndNumber, BeginNumber>
{
};


/**
* ������������ ������������� ������������������
*/
template
	<
		typename Integer, // ������������� ���
		Integer BeginNumber, // ����� �����, � �������� ������ �������� ��������� ������������������
		Integer EndNumber // �������� ����� �����, ������� ����� ���� ���������� �������������������
	>
using IncreasingIntegerSequence = typename IntegerSequenceGenerator<Integer, BeginNumber, EndNumber>::type;


/**
* ������������ ������������������ �������� (������� ����� ��� size_t, �� ��������� ���������� � 0)
*/
template
	<
		size_t Number // ���������� ����� � ������������������
	>
using IncreasingIndexSequence = IncreasingIntegerSequence<size_t, 0, Number - 1>;


#endif // INTEGER_SEQUENCE_H