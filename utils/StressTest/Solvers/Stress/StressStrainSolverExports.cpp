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

	/** Инициализация
	* @param params - параметры
	* @param links - связи
	* @params nElements - количество элементов
	* @params elements - координаты элементов (тройки чисел)
	*
	* @return указатель на решатель
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

	/** Установить ico
	* @param hsolver - дескриптор решателя
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


	/** Удалить решатель
	* @param hsolver - дескриптор решателя
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

	/** Получить адрес указателя на память
	* @param hsolver - дескриптор решателя
	*/
	DLL_FUNCTION
		float* GetMemoryPointer
		(
		const void* hsolver
		)
	{
		return ((StressStrainSolver*)hsolver)->GetMemoryPointer();
	}

	/** Получить размер памятм
	* @param hsolver - дескриптор решателя
	*/
	DLL_FUNCTION
		int GetMemorySize
		(
		const void* hsolver
		)
	{
		return ((StressStrainSolver*)hsolver)->GetMemorySize();
	}

	/** Получить адрес указателя на память
	* @param hsolver - дескриптор решателя
	*/
	DLL_FUNCTION
		void* GetAccelerations(const void* hsolver)
	{
		return (void*)((StressStrainSolver*)hsolver)->GetDataInternal(DT_Accelerations);
	}

	/** Добавление связей
	* @param hsolver - дескриптор решателя
	* @param links - массив 6*_nElements с признаками связей (1 - есть связь, 0 - нет связи)
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

	/** Добавить граничное условие
	* @param hsolver - дескриптор решателя
	* @param boundaryNodesIndices - индексы граничных точек
	* @param numberOfBoundaryNodes - количество граничных точек в частичной границе
	* @param numberOfBoundaryNodesInFullBoundary - количество граничных точек в полной границе
	* @param bcKind - тип граничных словий (=1 - сила, равномерно разделенная по точкам, 
	* задается вектором из 6 компонент - 3 поступательных компоненты вектора силы 
	* по координатным осям и три момента относительно осей)
	* @param bcParams - переменные граничных условий
	* =1: 6 параметров - компоненты силы Fx,Fy,Fz,Mx,My,Mz
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

	/** Добавить граничное условие
	* @param hsolver - дескриптор решателя
	* @param boundaryNodesIndices - индексы граничных точек
	* @param numberOfBoundaryNodes - количество граничных точек
	* @param bcKind - тип граничных словий (=1 - сила, равномерно разделенная по точкам, 
	* задается вектором из 6 компонент - 3 поступательных компоненты вектора силы 
	* по координатным осям и три момента относительно осей)
	* @param bcParams - переменные граничных условий
	* =1: 6 параметров - компоненты силы Fx,Fy,Fz,Mx,My,Mz
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

	/** Изменить парметр граничного условия
	* @param hsolver - дескриптор решателя
	* @param bcNumber - номер граничного условия (в порядке добавления),
	* @param bcParamNumber - номер параметра в параметрах граничных условий
	* @param bcParamValue - новое значение параметра
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
	* Обновить буфер
	* @param hsolver - дескриптор решателя
	* @param scaleFactor - масштабный коэффициент
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


	/** Расчет начальноо времени
	* @param nIterations - количество итераций
	*/
	DLL_FUNCTION
		void InitialSolve
		(
		void* hsolver
		)
	{
		((StressStrainSolver*)hsolver)->InitialSolve();
	}

	/** Расчет
	* @param nIterations - количество итераций
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

	/** Расчет первой стадии метода Рунге-Кутты
	* @param hsolver - дескриптор решателя
	*/
	DLL_FUNCTION
		void Solve1
		(
			void* hsolver
		)
	{
		((StressStrainSolver*)hsolver)->Solve1();
	}

	/** Расчет второй стадии метода Рунге-Кутты
	* @param hsolver - дескриптор решателя
	*/
	DLL_FUNCTION
		void Solve2
		(
			void* hsolver
		)
	{
		((StressStrainSolver*)hsolver)->Solve2();
	}

	/** Расчет третьей стадии метода Рунге-Кутты
	* @param hsolver - дескриптор решателя
	*/
	DLL_FUNCTION
		void Solve3
		(
			void* hsolver
		)
	{
		((StressStrainSolver*)hsolver)->Solve3();
	}

	/** Расчет четвертой стадии метода Рунге-Кутты
	* @param hsolver - дескриптор решателя
	*/
	DLL_FUNCTION
		void Solve4
		(
			void* hsolver
		)
	{
		((StressStrainSolver*)hsolver)->Solve4();
	}

	/** Расчет пятой стадии метода Рунге-Кутты
	* @param hsolver - дескриптор решателя
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