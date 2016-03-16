#pragma once


#include "CommonSolvers.h"

#include <string>


namespace SpecialSolversTest
{

	using namespace SpecialSolvers;


	// Вспомогательные функции создания сетки
	
	// Создает тестовую сетку nx x ny x nz с шагом step
	// количество узлов = nx x ny x nz
	// количество связей nLinks
	void CreateTestGrid
		(
			double*& nodes,
			int*& links,
			int& nLinks,
			int nx,
			int ny,
			int nz,
			double step,
			const std::string& fileName
		);

	void CreateTestGrid
		(
			const GridParams& gridParams,
			const std::string& rlcGridFileName
		);

}