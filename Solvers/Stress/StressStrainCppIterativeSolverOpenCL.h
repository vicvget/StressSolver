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

		// создает объект с заданными параметрами
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
	double _elasticModulus,			// модуль упругости (используется только для вычисления напряжений)
	double _dampingFactorAngular,	// приведенный коэффициент углового демпфирования
	double _dampingFactorLinear,	// приведенный коэффициент линейного демпфирования
	double _elasticFactorLinear,	// приведенный коэффициент линейной жесткости
	double _elasticFactorAngular	// приведенный коэффициент угловой жесткости*/

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
		* Расчет первой стадии метода Рунге-Кутты
		*/
		virtual
			void Solve1();

		/**
		* Расчет второй стадии метода Рунге-Кутты
		*/
		virtual
			void Solve2();

		/**
		* Расчет третьей стадии метода Рунге-Кутты
		*/
		virtual
			void Solve3();

		/**
		* Расчет четвертой стадии метода Рунге-Кутты
		*/
		virtual
			void Solve4();

		/**
		* Расчет пятой стадии метода Рунге-Кутты
		*/
		virtual
			void Solve5();

		/** Получить смещения
		* @param data - массив для записи смещений как скалярного параметра
		*/
		virtual
			void GetDisplacement
			(
				float* data
			);

		/** Получить напряжения по первой теории прочности
		* @param data - массив для записи напряжений как скалярного параметра
		*/
		virtual
			void GetStressesByFirstTheoryOfStrength
			(
				float* data
			);

		/** Получить напряжения по von Mises
		* @param data - массив для записи напряжений как скалярного параметра
		*/
		virtual	void GetStressesByVonMises(float* data);
		virtual
			void CalculateStrains(size_t side, double* shiftStrains, double* velocityStrains, size_t nodeId1, size_t nodeId2) const;
		virtual
			void CalculateStrainsUa(size_t side, double* shiftStrains, double* velocityStrains, size_t nodeId1, size_t nodeId2) const;
		/** Получить напряжения по von Mises
		* @param data - массив для записи напряжений как скалярного параметра
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
		// номер итерации расчетного цикла
		int _iterationNumber;

		void ApplyBoundary();

		void ApplyMass();

		void DebugCalculateForces();

};
#endif //OPENCL
}