#include "StressStrainSolverExports.h"
#include "StressStrainCppIterativeSolver.h"

#ifdef USE_KNC
#include "StressStrainCppIterativeSolverKNC.h"
#else
#include "StressStrainCppIterativeSolverAVX.h"
#include "StressStrainCppIterativeSolverFMA.h"
#endif
//#include "StressStrainCppIterativeSolverSSE.h"
//#include "../../Fcore/Exceptions/fcExceptions.h"

#include <memory>


//StressStrainSolver* _hsolver = NULL;


namespace Stress
{

	/** �������������
	* @param params - ���������
	* @param links - �����
	* @params nElements - ���������� ���������
	* @params elements - ���������� ��������� (������ �����)
	*
	* @return ��������� �� ��������
	*/
	DLL_FUNCTION
		void* Init
		(
		double* params,
		int* links,
		int nLinks,
		double *gridElements,
		int nElements,
		double gridStep,
		double timeStep,
		int procType,
		int numThreads,
		bool useFloat,
		int solverType
		)
	{	
		StressStrainSolver* hsolver = NULL;
		switch (solverType)
		{
		case 0:
			hsolver = new StressStrainCppIterativeSolver
				(
				params,
				links,
				nLinks,
				gridElements,
				nElements,
				gridStep,
				timeStep,
				numThreads,
				4
				);
			break;
#ifndef USE_KNC
		case 1:
			hsolver = new StressStrainCppIterativeSolverAVX
				(
				params,
				links,
				nLinks,
				gridElements,
				nElements,
				gridStep,
				timeStep,
				numThreads,
				4
				);
			break;
		case 2:
			hsolver = new StressStrainCppIterativeSolverFMA
			(
				params,
				links,
				nLinks,
				gridElements,
				nElements,
				gridStep,
				timeStep,
				numThreads,
				4
			);
			break;
#else
		case 3:
			hsolver = new StressStrainCppIterativeSolverKNC
			(
				params,
				links,
				nLinks,
				gridElements,
				nElements,
				gridStep,
				timeStep,
				numThreads,
				4
			);
			break;
#endif
		default:
			std::cout << "Unsupported solver type " << std::endl;
		}
		
		return (void*)hsolver;
	}

	/** ���������� ico
	* @param hsolver - ���������� ��������
	*/
	DLL_FUNCTION void InitIco
	(
		void* hsolver,
		const char* fileName,
		bool readIco,
		bool writeIco
	)
	{
		((StressStrainSolver*)hsolver)->InitIco(fileName, readIco, writeIco, 0);
	}

	DLL_FUNCTION
		bool ReadIco
		(
			void* hsolver,
			const char* fileName
		)
	{
		return ((StressStrainSolver*)hsolver)->ReadIco(fileName);
	}

	DLL_FUNCTION
		void WriteIco
		(
			void* hsolver,
			const char* fileName
		)
	{
		((StressStrainSolver*)hsolver)->WriteIco(fileName);
	}


	/** ������� ��������
	* @param hsolver - ���������� ��������
	*/
	DLL_FUNCTION
		void ReleaseMemory
		(
		void* &hsolver
		)
	{
		if (hsolver)
		{
			StressStrainSolver* stressSolver = (StressStrainSolver*)hsolver;
			delete stressSolver;
			hsolver = NULL;
		}
	}

	/** �������� ����� ��������� �� ������
	* @param hsolver - ���������� ��������
	*/
	DLL_FUNCTION
		float* GetMemoryPointer
		(
		const void* hsolver
		)
	{
		return ((StressStrainSolver*)hsolver)->GetMemoryPointer();
	}

	/** �������� ������ ������
	* @param hsolver - ���������� ��������
	*/
	DLL_FUNCTION
		int GetMemorySize
		(
		const void* hsolver
		)
	{
		return ((StressStrainSolver*)hsolver)->GetMemorySize();
	}

	/** �������� ����� ��������� �� ������
	* @param hsolver - ���������� ��������
	*/
	DLL_FUNCTION
		void* GetAccelerations(const void* hsolver)
	{
		return (void*)((StressStrainSolver*)hsolver)->GetDataInternal(DT_Accelerations);
	}

	/** ���������� ������
	* @param hsolver - ���������� ��������
	* @param links - ������ 6*_nElements � ���������� ������ (1 - ���� �����, 0 - ��� �����)
	*/
	DLL_FUNCTION
		void AddLinks
		(
		const void* hsolver,
		int* links
		)
	{
		((StressStrainSolver*)hsolver)->AddLinks(links);
	}

	/** �������� ��������� �������
	* @param hsolver - ���������� ��������
	* @param boundaryNodesIndices - ������� ��������� �����
	* @param numberOfBoundaryNodes - ���������� ��������� ����� � ��������� �������
	* @param numberOfBoundaryNodesInFullBoundary - ���������� ��������� ����� � ������ �������
	* @param bcKind - ��� ��������� ������ (=1 - ����, ���������� ����������� �� ������, 
	* �������� �������� �� 6 ��������� - 3 �������������� ���������� ������� ���� 
	* �� ������������ ���� � ��� ������� ������������ ����)
	* @param bcParams - ���������� ��������� �������
	* =1: 6 ���������� - ���������� ���� Fx,Fy,Fz,Mx,My,Mz
	*/
	DLL_FUNCTION
		void AddPartialBoundary
		(
		const void *hsolver,
		int* boundaryNodesIndices, 
		int numberOfBoundaryNodes,
		int numberOfBoundaryNodesInFullBoundary,
		int bcKind,
		double* bcParams
		)
	{
		((StressStrainSolver*)hsolver)->AddPartialBoundary
			(
			boundaryNodesIndices, 
			numberOfBoundaryNodes,
			numberOfBoundaryNodesInFullBoundary,
			bcKind,
			bcParams
			);
	}

	/** �������� ��������� �������
	* @param hsolver - ���������� ��������
	* @param boundaryNodesIndices - ������� ��������� �����
	* @param numberOfBoundaryNodes - ���������� ��������� �����
	* @param bcKind - ��� ��������� ������ (=1 - ����, ���������� ����������� �� ������, 
	* �������� �������� �� 6 ��������� - 3 �������������� ���������� ������� ���� 
	* �� ������������ ���� � ��� ������� ������������ ����)
	* @param bcParams - ���������� ��������� �������
	* =1: 6 ���������� - ���������� ���� Fx,Fy,Fz,Mx,My,Mz
	*/
	DLL_FUNCTION
		void AddBoundary
		(
		const void *hsolver,
		int* boundaryNodesIndices, 
		int numberOfBoundaryNodes,
		int bcKind,
		double* bcParams
		)
	{
		((StressStrainSolver*)hsolver)->AddBoundary
			(
			boundaryNodesIndices, 
			numberOfBoundaryNodes,
			bcKind,
			bcParams
			);
	}

	/** �������� ������� ���������� �������
	* @param hsolver - ���������� ��������
	* @param bcNumber - ����� ���������� ������� (� ������� ����������),
	* @param bcParamNumber - ����� ��������� � ���������� ��������� �������
	* @param bcParamValue - ����� �������� ���������
	*/
	DLL_FUNCTION
		void ChangeBoundaryParam
		(
		const void *hsolver,
		int bcNumber,
		int bcParamNumber,
		double bcParamValue
		)
	{
		((StressStrainSolver*)hsolver)->ChangeBoundaryParam
			(
			bcNumber,
			bcParamNumber,
			bcParamValue
			);	
	}

	DLL_FUNCTION
		double GetBoundaryParam
		(
		const void *hsolver,
		int bcNumber,
		int bcParamNumber
		)
	{
		return ((StressStrainSolver*)hsolver)->GetBoundaryParam
			(
			bcNumber,
			bcParamNumber
			);	
	}

	/**
	* �������� �����
	* @param hsolver - ���������� ��������
	* @param scaleFactor - ���������� �����������
	*/
	DLL_FUNCTION
	void UpdateBuffer
		(
			const void* hsolver,
			double scaleFactor
		)
	{
		((StressStrainSolver*)hsolver)->UpdateBuffer(scaleFactor);
	}


	/** ������ ��������� �������
	* @param nIterations - ���������� ��������
	*/
	DLL_FUNCTION
		void InitialSolve
		(
		void* hsolver
		)
	{
		((StressStrainSolver*)hsolver)->InitialSolve();
	}

	/** ������
	* @param nIterations - ���������� ��������
	*/
	DLL_FUNCTION
		void Solve
		(
			void* hsolver,
			int nIterations
		)
	{
		((StressStrainSolver*)hsolver)->Solve(nIterations);
	}

	/** ������ ������ ������ ������ �����-�����
	* @param hsolver - ���������� ��������
	*/
	DLL_FUNCTION
		void Solve1
		(
			void* hsolver
		)
	{
		((StressStrainSolver*)hsolver)->Solve1();
	}

	/** ������ ������ ������ ������ �����-�����
	* @param hsolver - ���������� ��������
	*/
	DLL_FUNCTION
		void Solve2
		(
			void* hsolver
		)
	{
		((StressStrainSolver*)hsolver)->Solve2();
	}

	/** ������ ������� ������ ������ �����-�����
	* @param hsolver - ���������� ��������
	*/
	DLL_FUNCTION
		void Solve3
		(
			void* hsolver
		)
	{
		((StressStrainSolver*)hsolver)->Solve3();
	}

	/** ������ ��������� ������ ������ �����-�����
	* @param hsolver - ���������� ��������
	*/
	DLL_FUNCTION
		void Solve4
		(
			void* hsolver
		)
	{
		((StressStrainSolver*)hsolver)->Solve4();
	}

	/** ������ ����� ������ ������ �����-�����
	* @param hsolver - ���������� ��������
	*/
	DLL_FUNCTION
		void Solve5
		(
			void* hsolver
		)
	{
		((StressStrainSolver*)hsolver)->Solve5();
	}

	DLL_FUNCTION
		double GetData
		(
		const void* hsolver,
		int nNode,
		int dir,
		int type
		)
	{

		return ((StressStrainSolver*)hsolver)->GetData
			(
			nNode,
			dir,
			type
			);
	}

	//DLL_FUNCTION
	//void SolveSystemOfLinearEquationsForStiffness()
	//{
	//	try
	//	{
	//		StiffnessMatrixPackedInCSCFormat stiffnessMatrixCSC;

	//		stiffnessMatrixCSC.ReadFromBinaryFile(StiffnessMatrixFileName + ".bin");

	//		StiffnessRHSVector stiffnessRHSVector = ReadStiffnessRHSVectorFromBinaryFile(StiffnessRHSFileName + ".bin");
	//		StiffnessRHSVector stiffnessResultsVector =
	//			StressStrainFortranStaticSolver2::SolveSystemOfLinearEquationsForStiffness
	//				(
	//					stiffnessMatrixCSC,
	//					stiffnessRHSVector
	//				);

	//		WriteStiffnessRHSVectorToTextFile("Results.txt", stiffnessResultsVector);
	//	}
	//	catch (const exceptions::CoreException& exception)
	//	{
	//		std::cout << exception.ToString() << std::endl;
	//	}
	//}

}