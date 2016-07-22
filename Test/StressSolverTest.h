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
		enum ECode;

		using namespace SpecialSolvers;
		using namespace SpecialSolvers::StressStrainStuff;


		SolverHandler MakeSolver(
			const GridParams& gridParams,
			const SpecialParams& specialParams,
			const IntegrationParams& integrationParams,
			const std::string& solverUid,
			EFACE sealedFace,
			EFACE forcedFace,
			double force,
			EDOF dof,
			const int solverType);

		SolverHandler MakePlateSolver(
			const GridParams& gridParams,
			const SpecialParams& specialParams,
			const IntegrationParams& integrationParams,
			const std::string& solverUid,
			double force,
			EDOF dof,
			const int solverType);

		SolverHandler MakeSolver
		(
			const GridParams& gridParams,
			const SpecialParams& specialParams,
			const IntegrationParams& integrationParams,
			const std::string& solverUid,
			const int solverType
		);

		void SetSealedForceBc(
			SolverHandler hStressSolver,
			EFACE faceSealed,
			EFACE faceForced,
			double force,
			EDOF dof,
			const GridParams& gridParams);

		void SetSealedForceForceBc(
			SolverHandler hStressSolver,
			EFACE faceSealed,
			EFACE faceForced,
			double force,
			EDOF dof,
			EDOF dof2,
			const GridParams& gridParams);

		void SetElementForceBc(
			SolverHandler hStressSolver,
			double torque,
			size_t nodeId,
			EDOF dof);

		void OverrideStiffness(
			SolverHandler hStressSolver,
			double elasticFactorLinear,
			double elasticFactorAngular,
			double dampingFactorLinear,
			double dampingFactorAngular,
			double stiffnessScale);

		void SetSealedForcePlateBc(
			SolverHandler hStressSolver,
			size_t side,
			double force,
			EDOF dof);


		void Test();
		void Test1x1x3(int solverType, ECode code);
		void Test10x1x3(int solverType, ECode code);
		void Test10x3x1(int solverType, ECode code);
		void Test10x3x3(int solverType, ECode code);
		void Test1x10x10(int solverType, ECode code);
		void Test1x11x11(int solverType, ECode code);

		void Test1x2x10(int solverType, ECode code);
		void Test2x1x10(int solverType, ECode code);
		void Test3x3x10(int solverType, ECode code);
		void Test2x2x5(int solverType,  ECode code);
		void TestSolveSystemOfLinearEquationsForStiffness();

	}

}