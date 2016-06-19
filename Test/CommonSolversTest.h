#pragma once


#include "CommonSolvers.h"

#include <string>


namespace SpecialSolversTest
{

	using namespace SpecialSolvers;


	// ��������������� ������� �������� �����
	
	// ������� �������� ����� nx x ny x nz � ����� step
	// ���������� ����� = nx x ny x nz
	// ���������� ������ nLinks
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