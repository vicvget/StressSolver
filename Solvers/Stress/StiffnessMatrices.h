#ifndef STIFFNESS_MATRICES_H

#define STIFFNESS_MATRICES_H

#include "StiffnessMatrix.h"
#include "StiffnessMatrixPackedInCoordinateFormat.h"
#include "StiffnessMatrixPackedInCSCFormat.h"
#include "../../AdditionalModules/Macros/ArgumentsListMacros.h"

#include <tuple>


/**
* ����� ������� ������ ���������
*/
template
	<
		class... Matrices // ������ ������ ���������
	>
class StiffnessMatricesTuple
{
public:

	// ����

	// ������� ��� ���� ������������� ������� - std::tuple
	using Tuple = std::tuple<Matrices...>;


	// ������������ � ����������

	/**
	* �����������
	* @param nodesCount ���������� ����� ��������� ������������� ����
	*/
	explicit
	StiffnessMatricesTuple
		(
			size_t nodesCount
		);

	/**
	* ����������� �����������
	* @param otherStiffnessMatricesTuple - ������ ������ ������ ���������, �� �������� ������������ �����������
	*/
	StiffnessMatricesTuple
		(
			StiffnessMatricesTuple&& otherStiffnessMatricesTuple
		);


	// ���������

	/**
	* �������� ������ ������ ���������
	* @return ������ ������ ���������
	*/
	const Tuple& GetTuple() const;

	/**
	* �������� ������ ������ ���������
	* @return ������ ������ ���������
	*/
	Tuple& GetTuple();


	// ������ ������� ������ ���������

	// ���������� ������ ������� ������ ���������
#	define STIFFNESS_MATRICES_FUNCTION_DECLARATION(Function, ReturnType, ...) \
	 \
	ReturnType Function \
		( \
			CREATE_PARAMETERS_LIST(__VA_ARGS__) \
		);

#	define STIFFNESS_MATRICES_FUNCTION STIFFNESS_MATRICES_FUNCTION_DECLARATION

#	include "StiffnessMatricesFunctions.h"

#	undef STIFFNESS_MATRICES_FUNCTION

#	undef STIFFNESS_MATRICES_FUNCTION_DECLARATION


private:

	// ������ ������ ���������
	Tuple _tuple;

};


// ������������ ������ ������ ���������
#define USED_STIFFNESS_MATRICES \
	/* StiffnessMatrix, */ \
	/* StiffnessMatrixPackedInCoordinateFormat, */ \
	StiffnessMatrixPackedInCSCFormat


// ����� ������������� StiffnessMatrices

// ��� ����� StiffnessMatrix � StiffnessMatrixPackedInCoordinateFormat
extern
template class StiffnessMatricesTuple<USED_STIFFNESS_MATRICES>;


// ������� ��� StiffnessMatricesTuple, ������������������� ������������� �������� ������ ���������
using StiffnessMatrices = StiffnessMatricesTuple<USED_STIFFNESS_MATRICES>;


#endif // STIFFNESS_MATRICES_H