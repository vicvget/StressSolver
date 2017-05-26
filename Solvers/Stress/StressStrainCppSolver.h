#pragma once

#include "StressStrainSolver.h"
#include "BoundaryParams.h"
//#include "RotationSolver.h"

#include <vector>
//#include <string>
#include "FTimer.h"

using std::string;
using std::vector;

namespace Stress
{
	/** �������� ��� ������� �������� ������ ��� �� ���������� ������������� �����
	*/
	class StressStrainCppSolver: public StressStrainSolver
	{
	public:

		void PrintTime();
		// ������� ������ � ��������� �����������
		StressStrainCppSolver
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

		virtual
			~StressStrainCppSolver();


#pragma region overriden

		virtual
			void AddLinks
			(
			const int* links
			);

		virtual
			void AddBoundary
			(
			int* boundaryNodesIndices,
			int numberOfBoundaryNodes,
			int bcKind,
			double* bcParams
			);

		void AddPartialBoundary
			(
			int* boundaryNodesIndicesInPart,
			int numberOfNodesInPart,
			int numberOfNodes,
			int bcKind,
			double* bcParams
			);

		virtual
			void ChangeBoundaryParam
			(
			const int bcNumber,
			const int bcParamNumber,
			const double bcParamValue
			);

		virtual
			double GetBoundaryParam
			(
			const int bcNumber,
			const int bcParamNumber
			);
		virtual
			void UpdateBuffer();
		virtual
			float UpdateBufferWithOutput();
#pragma endregion

		//protected:

		// ������� ������ ��������
		bool _isFirstSolution;

		// ������� ������ ��������
		vector<BoundaryParams> _boundaryParamsSet;

		// ���-�������, ������������ ����� �� solve
		double *_initX;
		double *_initDX;
		double *_hDDX1;
		double *_hDDX2;
		double *_hDDX3;

		// ���-�������
		double* _varX;		// �������� �� �����������
		double* _varDX;		// �������� �� �����������
		double* _varDDX;	// �������� �� 2 �����������

		double* _radiusVectors;
		int* _linkedElements;

		size_t _nVariables;				// ���������� �����������
		size_t _nIteration;				// ����� ��������
		int _numThreads;				// ����� ������������ ������� OpenMP

		//double _time;					// �����
		double _timeStep;				// ��� �� �������
		double _timeStep2;				// ��� �� ������� /2
		double _timeStep4;				// ��� �� ������� /4

		double _gridStep;				// ��� �����
		double _gridStep2;				// ��� ����� � ��������
		double _gridStep3;				// ��� ����� � ����

		double _elasticModulus;			// ������ ��������� (������������ ������ ��� ���������� ����������)

		double _dampingFactorAngular;	// ����������� ����������� �������� �������������
		double _dampingFactorLinear;	// ����������� ����������� ��������� �������������
		double _elasticFactorLinear;	// ����������� ����������� �������� ���������
		double _elasticFactorAngular;	// ����������� ����������� ������� ���������
		bool _isStiffnessOverriden;		// ������� ���������� �������� ��������� � �������������
		bool _isInertiaOverriden;		// ������� ���������� �������� ����� � �������� �������

		// ���������� ������������ ��� ���������� �� X,Y,Z
		double _stressScalingFactors[3];
		double _puassonFactor;


		//double _density;				// ��������� ���������
		double _cellMass;				// ����� ��������
		double _cellInertia;			// ������ ��������
		double _stiffScale;				// ���������������
		double _poissonRatio;			// ����������� ��������
		double _lameMatrix[6][6];		// ������� �������� �� ����������� � �����������

		FTimer _testTimer;			// ������ �������

		// ���������� ������� �� ���������� ����
		void FindStressStrainMatrix();

		// ���������� ������-������ ��������� ������� ����� 
		// ��� ���������� ����� � ����������� �� ����� ��������
		// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
		// @param side - ����� �����
		// @returns ������-������
		double* GetRadiusVector(size_t side) const;
		int GetLinkedElement(size_t elementId, size_t dof) const;

		void OverrideStiffness(
			double elasticModulus,
			double shearModulus,
			double dampingFactorLinear,
			double dampingFactorAngular,
			double stiffnessScale);

		void OverrideInertia(
			double mass,
			double inertia);
		
		static void CrossProduct(double* v1, double* v2, double* res);
		
		void OverrideScalingFactors(
			double stressScalingFactorX,
			double stressScalingFactorY,
			double stressScalingFactorZ);

		virtual void CalculateForces();
		void CalculateForcesOmp();
		virtual void ApplyBoundary();
		virtual void ApplyMass();
	};
};