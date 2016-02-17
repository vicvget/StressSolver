#ifndef StressStrainFortranIterativeSolverH

#define StressStrainFortranIterativeSolverH


#include "StressStrainFortranSolver.h"
#include "StiffnessMatrices.h"
#include "StiffnessRHS.h"

#include <cmath>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>


// ������������ �����, � ������� ��������� ������� ���������, �� ���������
const std::string StiffnessMatrixFileName = "StiffnessMatrix";

// ������������ �����, � ������� ��������� ������ ������ ������ ���� �� ������ ������� ���������, �� ���������
const std::string StiffnessRHSFileName = "StiffnessRHS";

// �����, �������������� � ��������� ���������� ����������� ��� ������������ ������� ���������
const double StiffnessEpsilon = 1e-5;


/**
* ��������������� ���������, �������������� ��� ������������� ����
*/
template
	<
		typename Type // ���������������� ���
	>
struct Identity
{
};


/** �������������� ��� �� ��������
*/
class StressStrainFortranIterativeSolver
	:public StressStrainFortranSolver
{
public:

	// ������� ������ � ��������� �����������
	StressStrainFortranIterativeSolver
		(
			double* params,
			int* links,
			int nLinks,
			double *nodes,
			int nNodes,
			double gridStep,
			double timeStep,
			int numThreads,
			int stride
		);

#pragma region overriden

	virtual
	void Solve
		(
			const int nIteratons
		);

	/**
	* ������ ������ ������ ������ �����-�����
	*/
	virtual
	void Solve1();

	/**
	* ������ ������ ������ ������ �����-�����
	*/
	virtual
	void Solve2();

	/**
	* ������ ������� ������ ������ �����-�����
	*/
	virtual
	void Solve3();

	/**
	* ������ ��������� ������ ������ �����-�����
	*/
	virtual
	void Solve4();

	/**
	* ������ ����� ������ ������ �����-�����
	*/
	virtual
	void Solve5();

	/** �������� ��������
	* @param data - ������ ��� ������ �������� ��� ���������� ���������
	*/
	virtual
	void GetDisplacement
		(
			float* data
		);

	/** �������� ���������� �� ������ ������ ���������
	* @param data - ������ ��� ������ ���������� ��� ���������� ���������
	*/
	virtual
	void GetStressesByFirstTheoryOfStrength
		(
			float* data
		);

	/** �������� ���������� �� von Mises
	* @param data - ������ ��� ������ ���������� ��� ���������� ���������
	*/
	virtual
	void GetStressesByVonMises
		(
			float* data
		);

#pragma endregion


//protected:

	// ������������ ��� ������� ��������� (��� ������ ������ ���������)
	using StiffnessMatrixType = StiffnessMatrixPackedInCSCFormat;


	// ����� �������� ���������� �����
	int _iterationNumber;


	void pravsubfl();

	void FindStressStrainMatrix
		(
			double (&matrix)[6][6]
		);

	void ApplyBoundary();

	void ApplyMass();

	/**
	* ������������ ���� �� ������ ������� ���������
	* @param stiffnessMatrixFileName - ����� ������������ ������ ��� ������ ����������� ������� ���������
	* � ���� ��������
	* @param stiffnessMatrixFileFormats - ������ �������� ������, � ������� ����� ��������
	* ����������� ������� ���������
	* @param stiffnessRHSFileName - ����� ������������ ������ ��� ������ ������������ �������
	* ������ ������ ���� �� ������ ������� ���������
	* @param writeStiffnessRHSVectorToTextFileFlag - ����, ������� �� ���������� ����������� ������
	* ������ ������ ���� �� ������ ������� ��������� � ��������� ����
	* @param stiffnessEpsilon - �����, �������������� � ��������� ���������� �����������
	* ��� ������������ ������� ���������
	*/
	void CreateSystemOfLinearEquationsForStiffness
		(
			const string& stiffnessMatrixFileName = StiffnessMatrixFileName,
			StiffnessMatrixFileFormats stiffnessMatrixFileFormats = //TO_INTEGRAL_TYPE(StiffnessMatrixFileFormat::BinaryFile),
				TO_INTEGRAL_TYPE(StiffnessMatrixFileFormat::BinaryFile | StiffnessMatrixFileFormat::TextFile),
			const string& stiffnessRHSFileName = StiffnessRHSFileName,
			//bool writeStiffnessRHSVectorToTextFileFlag = false,
			bool writeStiffnessRHSVectorToTextFileFlag = true,
			double stiffnessEpsilon = StiffnessEpsilon
		);

	/**
	* ������������ ������ ������� ������ ������������ �������� ������� ��������� ��������� �������������
	* @return ����������� ������ ������� ������ ������������ �������� ������� ��������� ��������� �������������
	*/
	IsSealedFlagsList FormIsSealedFlagsList() const;

	/**
	* ������������ ������� ���������
	* @param isSealedFlagsList - ������ ������� ������ ������������ �������� �������
	* ��������� ��������� �������������
	* @param stiffnessMatrixFileName - ����� ������������ ������ ��� ������ ����������� ������� ���������
	* � ���� ��������
	* @param stiffnessMatrixFileFormats - ������ �������� ������, � ������� ����� ��������
	* ����������� ������� ���������
	* @param stiffnessEpsilon - �����, �������������� � ��������� ���������� �����������
	* ��� ������������ ������� ���������
	*/
	void CreateStiffnessMatrix
		(
			const IsSealedFlagsList& isSealedFlagsList,
			const string& stiffnessMatrixFileName,
			StiffnessMatrixFileFormats stiffnessMatrixFileFormats,
			double stiffnessEpsilon
		);

	/**
	* ������������ ������� ���������
	* @param isSealedFlagsList - ������ ������� ������ ������������ �������� �������
	* ��������� ��������� �������������
	* @param stiffnessEpsilon - �����, �������������� � ��������� ���������� �����������
	* ��� ������������ ������� ���������
	* @return �������������� ������� ���������
	*/
	StiffnessMatrixType CreateStiffnessMatrix
		(
			const IsSealedFlagsList& isSealedFlagsList,
			double stiffnessEpsilon
		);

	/**
	* ������������ ������� ���������, ����������� � CSC (compressed sparse column) �������
	* @param isSealedFlagsList - ������ ������� ������ ������������ �������� �������
	* ��������� ��������� �������������
	* @param stiffnessEpsilon - �����, �������������� � ��������� ���������� �����������
	* ��� ������������ ������� ���������
	* @param ��������������� ��������, ����������� ��� ��������� ���� UsedStiffnessMatrixType
	* @return �������������� ������� ���������, ����������� � CSC (compressed sparse column) �������
	*/
	template
		<
			class UsedStiffnessMatrixType = StiffnessMatrixType // ������������ ��� ������� ���������
		>
	std::enable_if_t<std::is_same<UsedStiffnessMatrixType, StiffnessMatrixPackedInCSCFormat>::value, StiffnessMatrixPackedInCSCFormat>
	CreateStiffnessMatrixPackedInCSCFormat
		(
			const IsSealedFlagsList& isSealedFlagsList,
			double stiffnessEpsilon,
			Identity<UsedStiffnessMatrixType> = {}
		);

	/**
	* ������������ ������� ���������, ����������� � CSC (compressed sparse column) �������
	* @param isSealedFlagsList - ������ ������� ������ ������������ �������� �������
	* ��������� ��������� �������������
	* @param stiffnessEpsilon - �����, �������������� � ��������� ���������� �����������
	* ��� ������������ ������� ���������
	* @param ��������������� ��������, ����������� ��� ��������� ���� UsedStiffnessMatrixType
	* @return �������������� ������� ���������, ����������� � CSC (compressed sparse column) �������
	*/
	template
		<
			class UsedStiffnessMatrixType = StiffnessMatrixType // ������������ ��� ������� ���������
		>
	std::enable_if_t<std::is_same<UsedStiffnessMatrixType, StiffnessMatrices>::value, StiffnessMatrixPackedInCSCFormat>
	CreateStiffnessMatrixPackedInCSCFormat
		(
			const IsSealedFlagsList& isSealedFlagsList,
			double stiffnessEpsilon,
			Identity<UsedStiffnessMatrixType> = {}
		);

	/**
	* ��������� ������� ���������
	* @param stiffnessMatrix - ������� ���������, ������� ����������� ����������
	* @param isSealedFlagsList - ������ ������� ������ ������������ �������� �������
	* ��������� ��������� �������������
	* @param stiffnessEpsilon - �����, �������������� � ��������� ���������� �����������
	* ��� ������������ ������� ���������
	*/
	void FindStiffnessMatrix
		(
			StiffnessMatrixType& stiffnessMatrix,
			const IsSealedFlagsList& isSealedFlagsList,
			double stiffnessEpsilon
		);

	/**
	* ��������� ������� ��������� ���������� ���� �������� ��������� � �������� � ������ ��������
	* @param stiffnessMatrix - ������� ���������, ������� ����������� ����������
	* @param nodeId - ������ ������� ��������
	* @param isSealedFlagsList - ������ ������� ������ ������������ �������� �������
	* ��������� ��������� �������������
	* @param stiffnessEpsilon - �����, �������������� � ��������� ���������� �����������
	* ��� ������������ ������� ���������
	*/
	void FindStiffnessMatrixForAdjElements
		(
			StiffnessMatrixType& stiffnessMatrix,
			int nodeId,
			const IsSealedFlagsList& isSealedFlagsList,
			double stiffnessEpsilon
		);

	/**
	* ��������� ������� ��������� ���������� ��������� �������� � �������� � ������ ��������
	* @param stiffnessMatrix - ������� ���������, ������� ����������� ����������
	* @param nodeId - ������ ������� ��������
	* @param adjNodeId - ������ ��������� ��������
	* @param isSealedFlagsList - ������ ������� ������ ������������ �������� �������
	* ��������� ��������� �������������
	* @param stiffnessEpsilon - �����, �������������� � ��������� ���������� �����������
	* ��� ������������ ������� ���������
	*/
	void FindStiffnessMatrixForElement
		(
			StiffnessMatrixType& stiffnessMatrix,
			int nodeId,
			int adjNodeId,
			const IsSealedFlagsList& isSealedFlagsList,
			double stiffnessEpsilon
		);

	/**
	* ��������� ������� ��������� ���������� ��������, ��������� ������� ������� �������� ����������
	* @param stiffnessMatrix - ������� ���������, ������� ����������� ����������
	* @param nodeId - ������ ��������, ��������� ������� ������� �������� ����������
	* @param isSealedFlagsList - ������ ������� ������ ������������ �������� �������
	* ��������� ��������� �������������
	*/
	void FindStiffnessMatrixForSealedElement
		(
			StiffnessMatrixType& stiffnessMatrix,
			int nodeId,
			const IsSealedFlagsList& isSealedFlagsList
		);

	/**
	* ��������� ������� ��������� ���������� ��������� �������� � �������� � ������ ��������
	* � ������ �������� �������
	* @param stiffnessMatrix - ������� ���������, ������� ����������� ����������
	* @param nodeId - ������ ������� ��������
	* @param dof - ������ ������ ������� �������
	* @param adjNodeId - ������ ��������� ��������
	* @param isSealedFlagsList - ������ ������� ������ ������������ �������� �������
	* ��������� ��������� �������������
	* @param stiffnessEpsilon - �����, �������������� � ��������� ���������� �����������
	* ��� ������������ ������� ���������
	*/
	void FindStiffnessMatrixForDof
		(
			StiffnessMatrixType& stiffnessMatrix,
			int nodeId,
			int dof,
			int adjNodeId,
			const IsSealedFlagsList& isSealedFlagsList,
			double stiffnessEpsilon
		);

	/**
	* �������� �������� ������� ��������� ��� �������� � ������ ��������, ������� �������� �����������
	* �� ����������� changedDof ���� � �������� � �������� changedNodeId
	* @param changedNodeId - ������ ��������, � �������� �������������� ����
	* @param changedDof - ������ ������� �������, �� ������� �������������� ����
	* @param nodeId - ������ ������� ��������
	* @param stiffnessEpsilon - �����, �������������� � ��������� ���������� �����������
	* ��� ������������ ������� ���������
	*/
	void UpdateForcesForNode
		(
			int changedNodeId,
			int changedDof,
			int nodeId,
			double stiffnessEpsilon
		);

	/**
	* ������������ ������ ������ ������ ���� �� ������ ������� ���������
	* @param isSealedFlagsList - ������ ������� ������ ������������ �������� �������
	* ��������� ��������� �������������
	* @param stiffnessRHSFileName - ����� ������������ ������ ��� ������ ������������ �������
	* ������ ������ ���� �� ������ ������� ���������
	* @param writeStiffnessRHSVectorToTextFileFlag - ����, ������� �� ���������� ����������� ������
	* ������ ������ ���� �� ������ ������� ��������� � ��������� ����
	*/
	void CreateStiffnessRHSVector
		(
			const IsSealedFlagsList& isSealedFlagsList,
			const string& stiffnessRHSFileName,
			bool writeStiffnessRHSVectorToTextFileFlag
		)	const;

	/**
	* ������������ ������ ������ ������ ���� �� ������ ������� ���������
	* @param isSealedFlagsList - ������ ������� ������ ������������ �������� �������
	* ��������� ��������� �������������
	* @return �������������� ������ ������ ������ ���� �� ������ ������� ���������, ����������� ����������
	*/
	StiffnessRHSVector CreateStiffnessRHSVector
		(
			const IsSealedFlagsList& isSealedFlagsList
		)	const;

	/**
	* ��������� ���������� ������ ������ ������ ���� �� ������ ������� ���������
	* @param stiffnessRHSVector - ������ ������ ������ ���� �� ������ ������� ���������, ����������� ����������
	*/
	void FindStiffnessRHSVector
		(
			StiffnessRHSVector& stiffnessRHSVector
		)	const;

};

/**
* ������������ ������� ���������, ����������� � CSC (compressed sparse column) �������
* @param isSealedFlagsList - ������ ������� ������ ������������ �������� �������
* ��������� ��������� �������������
* @param stiffnessEpsilon - �����, �������������� � ��������� ���������� �����������
* ��� ������������ ������� ���������
* @param ��������������� ��������, ����������� ��� ��������� ���� UsedStiffnessMatrixType
* @return �������������� ������� ���������, ����������� � CSC (compressed sparse column) �������
*/
template
	<
		class UsedStiffnessMatrixType
	>
std::enable_if_t<std::is_same<UsedStiffnessMatrixType, StiffnessMatrixPackedInCSCFormat>::value, StiffnessMatrixPackedInCSCFormat>
StressStrainFortranIterativeSolver::CreateStiffnessMatrixPackedInCSCFormat
	(
		const IsSealedFlagsList& isSealedFlagsList,
		double stiffnessEpsilon,
		Identity<UsedStiffnessMatrixType>
	)
{
	return CreateStiffnessMatrix(isSealedFlagsList, stiffnessEpsilon);
}

/**
* ������������ ������� ���������, ����������� � CSC (compressed sparse column) �������
* @param isSealedFlagsList - ������ ������� ������ ������������ �������� �������
* ��������� ��������� �������������
* @param stiffnessEpsilon - �����, �������������� � ��������� ���������� �����������
* ��� ������������ ������� ���������
* @param ��������������� ��������, ����������� ��� ��������� ���� UsedStiffnessMatrixType
* @return �������������� ������� ���������, ����������� � CSC (compressed sparse column) �������
*/
template
	<
		class UsedStiffnessMatrixType
	>
std::enable_if_t<std::is_same<UsedStiffnessMatrixType, StiffnessMatrices>::value, StiffnessMatrixPackedInCSCFormat>
StressStrainFortranIterativeSolver::CreateStiffnessMatrixPackedInCSCFormat
	(
		const IsSealedFlagsList& isSealedFlagsList,
		double stiffnessEpsilon,
		Identity<UsedStiffnessMatrixType>
	)
{
	return std::get<0>(CreateStiffnessMatrix(isSealedFlagsList, stiffnessEpsilon).GetTuple());
}

#endif