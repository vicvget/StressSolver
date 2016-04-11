#pragma once

#include "StressSolver.h"


namespace SpecialSolversTest
{

	enum EDOF
	{
		dof_x = 0,
		dof_y,
		dof_z,
		dof_rx,
		dof_ry,
		dof_rz
	};

	enum EFACE
	{
		face_left = 0,
		face_right,
		face_top,
		face_bottom,
		face_front,
		face_back
	};

	namespace StressStrainStuff
	{

		using namespace SpecialSolvers;
		using namespace SpecialSolvers::StressStrainStuff;


		StressStrainSolver MakeSolver(
			const GridParams& gridParams,
			const SpecialParams& specialParams,
			const IntegrationParams& integrationParams,
			const std::string& fileRlc,
			EFACE sealedFace,
			EFACE forcedFace,
			double force,
			EDOF dof,
			const int solverType);

		StressStrainSolver MakeSolver
		(
			const GridParams& gridParams,
			const SpecialParams& specialParams,
			const IntegrationParams& integrationParams,
			const std::string& fileRlc,
			const int solverType
		);

		void SetSealedForceBc(
			StressStrainSolver hStressSolver,
			EFACE faceSealed,
			EFACE faceForced,
			double force,
			EDOF dof,
			const GridParams& gridParams);

		void SetSealedForceForceBc(
			StressStrainSolver hStressSolver,
			EFACE faceSealed,
			EFACE faceForced,
			double force,
			EDOF dof,
			EDOF dof2,
			const GridParams& gridParams);

		void SetElementForceBc(
			StressStrainSolver hStressSolver,
			double torque,
			size_t nodeId,
			EDOF dof);

		void Test();
		void Test3x1x1(int solverType = 0);
		void Test1x3x10(int solverType = 0);
		void Test3x3x10(int solverType = 0);
		void TestSolveSystemOfLinearEquationsForStiffness();

	}

}