#ifdef OPENCL
#include <CL\cl.h>
#endif // OPENCL



#pragma once

#include "StressStrainCppSolver.h"
#include "StressStrainCppIterativeSolver.h"
#include <cmath>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <immintrin.h>
#pragma once

using std::string;
using std::vector;

namespace Stress
{
	struct double3
	{
		double x, y, z, w;
		
	};

	double3 operator-(double3 a,double3 b);
	double3 operator+(double3 a, double3 b);
	double3 operator-(double3 a);

	void mulPrivateMM(double* a, double* b, double * res);

	void mulPrivateMV(double * a, double* b, double* res);

	typedef unsigned int uint;
#ifdef OPENCL
	class StressStrainCppIterativeSolverOpenCL
		:
		public StressStrainCppIterativeSolver
	{
	public:

		// ������� ������ � ��������� �����������
		StressStrainCppIterativeSolverOpenCL
		(
			double* params,
			int* links,
			int nLinks,
			double *coordinates,
			int nElements,
			double gridStep,
			double timeStep,
			int numThreads,
			int stride
		);

		virtual
			~StressStrainCppIterativeSolverOpenCL();

#pragma region overriden

		virtual
			void InitialSolve
			(
			);


		//OpenCL buffers
		/*	global double * _dataInternal,
	global double * _dataRotationMatrix,
	constant size_t * _linkedElements,
	global double* radiusVectors,
	global double* boundaryForces,
	global size_t* boundaryFixed,
	global double* _elementStressFactorCache, //
	global double* _stressScalingFactor,
	global double* _stress,
	global double* _coordinates,
	double _elasticModulus,			// ������ ��������� (������������ ������ ��� ���������� ����������)
	double _dampingFactorAngular,	// ����������� ����������� �������� �������������
	double _dampingFactorLinear,	// ����������� ����������� ��������� �������������
	double _elasticFactorLinear,	// ����������� ����������� �������� ���������
	double _elasticFactorAngular	// ����������� ����������� ������� ���������*/

        __declspec(align(32))cl_mem DataInternalBuffer;
		__declspec(align(32))cl_mem DataRotationMatrixBuffer;
		__declspec(align(32))cl_mem LinkedElemenetsBuffer;
		__declspec(align(32))cl_mem RadiusVectorsBuffer;
		__declspec(align(32))cl_mem BoundaryForcesBuffer;
		__declspec(align(32))cl_mem BoundaryFixedBuffer;
		__declspec(align(32))cl_mem ElementStressFactorCacheBuffer;
		__declspec(align(32))cl_mem StressScalingFactorBuffer;
		__declspec(align(32))cl_mem StressBuffer;
		__declspec(align(32))cl_mem CoordinatesBuffer;

		__declspec(align(32))cl_double* DataRotationMatrix;
		__declspec(align(32))cl_double* Stress;
		__declspec(align(32))cl_double* DataInternal;

		void InitOpenCLContext();

		void InitBuffers();

		void UploadPartialOpenCLBuffers();

		void UploadAllOpenCLBuffers();

		void DownloadOpenCLBuffers();

		virtual	void Solve(const int nIteratons);
		void clCalculate();

		void clInitKernel();
		virtual void SolveFull(const int nIteratons);
		cl_context clContext;
		cl_device_id deviceId;
		cl_program program;
		cl_kernel kernel;
		cl_command_queue command_queue;
		cl_device_info device_info;
		size_t max_all_work_items;

		cl_double* allBoundaryForces;
		cl_uint* FreeDegrees;
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
		virtual	void GetStressesByVonMises(float* data);
		virtual
			void CalculateStrains(size_t side, double* shiftStrains, double* velocityStrains, size_t nodeId1, size_t nodeId2) const;
		virtual
			void CalculateStrainsUa(size_t side, double* shiftStrains, double* velocityStrains, size_t nodeId1, size_t nodeId2) const;
		/** �������� ���������� �� von Mises
		* @param data - ������ ��� ������ ���������� ��� ���������� ���������
		*/
		virtual void GetStressesX(float* data);
		virtual void GetStressesY(float* data);
		virtual void GetStressesZ(float* data);

		virtual void GetStressesXY(float* data);
		virtual void GetStressesXZ(float* data);
		virtual void GetStressesYZ(float* data);

#pragma endregion

		double df[12]; // debug

	protected:
		// ����� �������� ���������� �����
		int _iterationNumber;

		void ApplyBoundary();

		void ApplyMass();

		void DebugCalculateForces();

};
#endif //OPENCL
}