#include "StressStrainFortranStaticSolver2.h"

#include "PardisoSolver.h"
#include "../../Fcore/Exceptions/fcExceptions.h"


StressStrainFortranStaticSolver2::StressStrainFortranStaticSolver2
	(
		double* params, 
		int* links, 
		int nLinks, 
		double *nodes, 
		int nNodes, 
		double gridStep, 
		double timeStep,
		int numThreads
	)
	:
		StressStrainFortranIterativeSolver
			(
				params,
				links, 
				nLinks, 
				nodes, 
				nNodes, 
				gridStep,
				timeStep,
				numThreads,
				3
			)
{
}

// virtual
void StressStrainFortranStaticSolver2::Solve
	(
		const int nIterations
	)
{
	try
	{
		StiffnessMatrixPackedInCSCFormat stiffnessMatrixCSC;
		StiffnessRHSVector stiffnessRHSVector;

		CreateSystemOfLinearEquationsForStiffness(stiffnessMatrixCSC, stiffnessRHSVector);

		StiffnessRHSVector stiffnessResultsVector =
			SolveSystemOfLinearEquationsForStiffness(stiffnessMatrixCSC, stiffnessRHSVector);

		UseSolutionOfSystemOfLinearEquationsForStiffness(stiffnessResultsVector);
	}
	catch (const exceptions::CoreException& exception)
	{
		std::cout << exception.ToString() << std::endl;
	}
}

/**
* Решить СЛАУ на основе матрицы жесткости
* @param stiffnessMatrixCSC - матрица жесткости в формате CSC
* @param stiffnessRHSVector - вектор правых частей СЛАУ на основе матрицы жесткости
* @return рассчитанный вектор результатов СЛАУ на основе матрицы жесткости
*/
StiffnessRHSVector StressStrainFortranStaticSolver2::SolveSystemOfLinearEquationsForStiffness
	(
		const StiffnessMatrixPackedInCSCFormat& stiffnessMatrixCSC,
		const StiffnessRHSVector& stiffnessRHSVector
	)
{
	if (stiffnessMatrixCSC.GetMatrixSize() != stiffnessRHSVector.size())
	{
		exceptions::ThrowMessage("Dimensions of stiffness matrix and RHS vector are not equal!");
	}

	StiffnessRHSVector stiffnessResultsVector(stiffnessRHSVector.size());

	std::cout << "Init solution" << std::endl;

#ifdef USE_MKL
	const size_t matrixSizePlus1 = stiffnessMatrixCSC.GetMatrixSize() + 1;
	auto iaCSC = std::make_unique<int[]>(matrixSizePlus1);
	auto jaCSC = std::make_unique<int[]>(stiffnessMatrixCSC.GetMatrixNonZeroElementsCount());

	stiffnessMatrixCSC.CopyColumnFirstElements(iaCSC.get());
	stiffnessMatrixCSC.CopyRowIds(jaCSC.get());

	PardisoSolver _pardisoSolver;

	if
		(
			!_pardisoSolver.Init
				(
					stiffnessMatrixCSC.GetMatrixSize(),
					iaCSC.get(),
					jaCSC.get(),
					const_cast<double*>(stiffnessMatrixCSC.GetElements())
				)
		)
	{
		exceptions::ThrowMessage("PARDISO INIT ERROR!");
	}
	std::cout << "Solution" << std::endl;
	if (!_pardisoSolver.Solve(const_cast<double*>(stiffnessRHSVector.data()), stiffnessResultsVector.data()))
	{
		exceptions::ThrowMessage("PARDISO SOLVE ERROR!");
	}
	std::cout << "End of solution" << std::endl;
#else
	exceptions::ThrowMessage("MKL IS NOT USED (#define USE_MKL)");
#endif

	return stiffnessResultsVector;
}

/**
* Сформировать СЛАУ на основе матрицы жесткости
* @param stiffnessMatrixCSC - сформированная матрица жесткости в формате CSC
* @param stiffnessRHSVector - сформированный вектор правых частей СЛАУ на основе матрицы жесткости
* @param stiffnessEpsilon - число, использующееся в численном нахождении производной
* при формировании матрицы жесткости
*/
// static
void StressStrainFortranStaticSolver2::CreateSystemOfLinearEquationsForStiffness
	(
		StiffnessMatrixPackedInCSCFormat& stiffnessMatrixCSC,
		StiffnessRHSVector& stiffnessRHSVector,
		double stiffnessEpsilon
	)
{
	const IsSealedFlagsList isSealedFlagsList = FormIsSealedFlagsList();

	stiffnessRHSVector = CreateStiffnessRHSVector(isSealedFlagsList);
	stiffnessMatrixCSC = CreateStiffnessMatrixPackedInCSCFormat(isSealedFlagsList, stiffnessEpsilon);
}

/**
* Использовать решение СЛАУ на основе матрицы жесткости для заполнения соответствующих массивов решателя
* @param stiffnessResultsVector - вектор результатов СЛАУ на основе матрицы жесткости
*/
void StressStrainFortranStaticSolver2::UseSolutionOfSystemOfLinearEquationsForStiffness
	(
		const StiffnessRHSVector& stiffnessResultsVector
	)
{
	std::copy_n(stiffnessResultsVector.data(), stiffnessResultsVector.size(), _stress);

	size_t index{};

	for (auto resultsElement : stiffnessResultsVector)
	{
		_dataInternal[index] += resultsElement;
		++index;
	}
}