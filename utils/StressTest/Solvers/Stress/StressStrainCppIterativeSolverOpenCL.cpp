#include <CL\cl.h>
#include "StressStrainCppIterativeSolverOpenCL.h"
#include "StressStrainCppSolver.h"
#include "../../AdditionalModules/fmath/Vector3.h"
#include "Common.h"
#include "../../AdditionalModules/fmath/Matrix3x3.h"
#include "../../AdditionalModules/fmath/Matrix3x4.h"

#define MAX_SOURCE_SIZE 50000

//#define USE_EXP \\ set this flag when targeting experimental platforms. Otherwise comment it out

#ifdef NON_BDW
#define PLATFROM_NAME "Intel(R) OpenCL"
#else
#define PLATFROM_NAME "Experimental OpenCL 2.1 CPU Only Platform"
#endif

using namespace MathHelpers;
namespace Stress
{

	Stress::StressStrainCppIterativeSolverOpenCL::StressStrainCppIterativeSolverOpenCL
	(
		double * params,
		int * links,
		int nLinks,
		double * coordinates,
		int nElements,
		double gridStep,
		double timeStep,
		int numThreads,
		int stride
	) :
		StressStrainCppIterativeSolver
		(
			params,
			links,
			nLinks,
			coordinates,
			nElements,
			gridStep,
			timeStep,
			numThreads,
			stride
		)
	{
		SetUid("OpenCL_Solver");
		std::cout << "OpenCL SOLVER" << std::endl << std::flush;
		allBoundaryForces = 0;
		FreeDegrees = 0;
	}

	Stress::StressStrainCppIterativeSolverOpenCL::~StressStrainCppIterativeSolverOpenCL()
	{

	}

	void Stress::StressStrainCppIterativeSolverOpenCL::InitialSolve()
	{
		_rotationSolver->InitialSolve();
		InitOpenCLContext();
		clInitKernel();
		InitBuffers();
		StressStrainCppIterativeSolver::CalculateForces();
		UploadAllOpenCLBuffers();
	}

	void Stress::StressStrainCppIterativeSolverOpenCL::InitOpenCLContext()
	{		
		cl_int code = 0;
		cl_platform_id platformId;
		cl_uint numPlatforms = 0;
		cl_uint numDevices = 0;
#ifndef USE_EXP

		code = clGetPlatformIDs(1,&platformId,&numPlatforms);
		if (code != 0)
		{
			fprintf(stderr, "Failed to get plantform ids.\n");
			exit(code);
		}
#else

		cl_platform_id * platforms = NULL;
		char vendor_name[128] = { 0 };
		cl_uint num_platforms = 0;

		// get number of available platforms
		cl_int err = clGetPlatformIDs(0, NULL, &num_platforms);
		if (CL_SUCCESS != err)
		{
			exit(err);
		}
		platforms = (cl_platform_id*)malloc(sizeof(platformId)* num_platforms);
		if (NULL == platforms)
		{
			// handle error
		}
		err = clGetPlatformIDs(num_platforms, platforms, NULL);
		if (CL_SUCCESS != err)
		{
			// handle error
		}
		// return the OpenCL 2.1 development environment platform
		for (cl_uint ui = 0; ui< num_platforms; ++ui)
		{
			err = clGetPlatformInfo(platforms[ui],
				CL_PLATFORM_NAME,
				128 * sizeof(char),
				vendor_name,
				NULL);
			if (CL_SUCCESS != err)
			{
				// handle error
			}
			if (vendor_name != NULL)
			{
				if (!strcmp(vendor_name, PLATFROM_NAME))
				{
					platformId = platforms[ui];
					break;
				}
			}
		}
#endif
		code = clGetDeviceIDs(platformId, CL_DEVICE_TYPE_GPU, 1, &deviceId, &numDevices);
		if (code != 0)
		{
			fprintf(stderr, "Failed to get device ids. %d\n",code);
			exit(code);
		}
		clContext = clCreateContext(NULL, 1, &deviceId, NULL, NULL, &code);
		if (code != 0)
		{
			fprintf(stderr, "Failed to create cotext. %d\n", code);
			exit(code);
		}
	}

	void Stress::StressStrainCppIterativeSolverOpenCL::InitBuffers()
	{
		int ret = 0;
		DataInternalBuffer = clCreateBuffer(clContext, CL_MEM_READ_WRITE, 6 * sizeof(cl_double)*vecStride*_nElements, NULL, &ret);
		if (ret != 0)
		{
			fprintf(stderr, "Failed to create buffers.\n");
			exit(ret);
		}
		DataRotationMatrixBuffer = clCreateBuffer(clContext, CL_MEM_READ_WRITE, matStride * sizeof(cl_double)*_nElements, NULL, &ret);
		if (ret != 0)
		{
			fprintf(stderr, "Failed to create buffers.\n");
			exit(ret);
		}

		LinkedElemenetsBuffer = clCreateBuffer(clContext, CL_MEM_READ_WRITE, sizeof(cl_int)*_nElements * 6, NULL, &ret);
		if (ret != 0)
		{
			fprintf(stderr, "Failed to create buffers.\n");
			exit(ret);
		}
		RadiusVectorsBuffer = clCreateBuffer(clContext, CL_MEM_READ_WRITE, vecStride * 6 * sizeof(cl_double), NULL, &ret);
		if (ret != 0)
		{
			fprintf(stderr, "Failed to create buffers.\n");
			exit(ret);
		}
		BoundaryForcesBuffer = clCreateBuffer(clContext, CL_MEM_READ_WRITE, vecStride2 * sizeof(cl_double) * _nElements, NULL, &ret);
		if (ret != 0)
		{
			fprintf(stderr, "Failed to create buffers.\n");
			exit(ret);
		}
		BoundaryFixedBuffer = clCreateBuffer(clContext, CL_MEM_READ_WRITE, 6 * sizeof(cl_int)*_nElements, NULL, &ret);
		if (ret != 0)
		{
			fprintf(stderr, "Failed to create buffers.\n");
			exit(ret);
		}
		ElementStressFactorCacheBuffer = clCreateBuffer(clContext, CL_MEM_READ_WRITE, vecStride * sizeof(cl_double)*_nElements, NULL, &ret);
		if (ret != 0)
		{
			fprintf(stderr, "Failed to create buffers.\n");
			exit(ret);
		}
		StressScalingFactorBuffer = clCreateBuffer(clContext, CL_MEM_READ_WRITE, vecStride * sizeof(cl_double)*_nElements, NULL, &ret);
		if (ret != 0)
		{
			fprintf(stderr, "Failed to create buffers.\n");
			exit(ret);
		}
		StressBuffer = clCreateBuffer(clContext, CL_MEM_READ_WRITE, 2 * vecStride * sizeof(cl_double)*_nElements, NULL, &ret);
		if (ret != 0)
		{
			fprintf(stderr, "Failed to create buffers.\n");
			exit(ret);
		}
		CoordinatesBuffer = clCreateBuffer(clContext, CL_MEM_READ_WRITE, 3 * sizeof(cl_double)*_nElements, NULL, &ret);
		if (ret != 0)
		{
			fprintf(stderr, "Failed to create buffers.\n");
			exit(ret);
		}
	}
	void Stress::StressStrainCppIterativeSolverOpenCL::UploadPartialOpenCLBuffers()
	{
		int ret = 0;

		DataInternal = (cl_double*)_dataInternal;
		ret = clEnqueueWriteBuffer(command_queue, DataInternalBuffer, CL_FALSE, 0, 6 * sizeof(cl_double) * vecStride * _nElements, DataInternal, 0, NULL, NULL);
		if (ret != 0)
		{
			fprintf(stderr, "Failed to write buffers.\n");
			exit(ret);
		}
		ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&DataInternalBuffer);
		if (ret != 0)
		{
			fprintf(stderr, "Failed with args setup.\n");
			exit(ret);
		}
		DataRotationMatrix = (cl_double*)_dataRotationMtx;
		if (ret != 0)
		{
			fprintf(stderr, "Failed to create buffers.\n");
			exit(ret);
		}
		ret = clEnqueueWriteBuffer(command_queue, DataRotationMatrixBuffer, CL_FALSE, 0, matStride * sizeof(cl_double)*_nElements, DataRotationMatrix, 0, NULL, NULL);
		if (ret != 0)
		{
			fprintf(stderr, "Failed to write buffers.\n");
			exit(ret);
		}
		ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&DataRotationMatrixBuffer);
		if (ret != 0)
		{
			fprintf(stderr, "Failed with args setup.\n");
			exit(ret);
		}
		Stress = _stress;
		ret = clEnqueueWriteBuffer(command_queue, StressBuffer, CL_TRUE, 0, 2 * vecStride * sizeof(cl_double)*_nElements, Stress, 0, NULL, NULL);
		ret = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *)&StressBuffer);

	}

	void Stress::StressStrainCppIterativeSolverOpenCL::UploadAllOpenCLBuffers()
	{
		//в девайс
		UploadPartialOpenCLBuffers();

		int ret = 0;
		

		cl_uint* LinkedElements = (cl_uint*)_linkedElements;
		ret = clEnqueueWriteBuffer(command_queue, LinkedElemenetsBuffer, CL_FALSE, 0, sizeof(cl_int)*_nElements * 6, LinkedElements, 0, NULL, NULL);
		ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&LinkedElemenetsBuffer);
		
		cl_double* RadiusVectors = (cl_double*)_radiusVectors;
		ret = clEnqueueWriteBuffer(command_queue, RadiusVectorsBuffer, CL_FALSE, 0, vecStride * 6* sizeof(cl_double), RadiusVectors, 0, NULL, NULL);
		ret = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&RadiusVectorsBuffer);

		if(allBoundaryForces==0)
			allBoundaryForces = (double*)malloc(vecStride2*sizeof(double)*_nElements);

		//1 степень свободы отпущена 0 закреплена
		if(FreeDegrees==0)
			FreeDegrees = (cl_uint*)malloc(6*sizeof(cl_uint)*_nElements);


		for (int i = 0; i < _nElements; i++)
		{
			allBoundaryForces[i*vecStride2] = 0;
			allBoundaryForces[i*vecStride2+1] = 0;
			allBoundaryForces[i*vecStride2+2] = 0;
			allBoundaryForces[i*vecStride2+3] = 0;
			allBoundaryForces[i*vecStride2+4] = 0;
			allBoundaryForces[i*vecStride2 + 5] = 0;
			allBoundaryForces[i*vecStride2 + 6] = 0;
			allBoundaryForces[i*vecStride2 + 7] = 0;

			FreeDegrees[i * 6] = 1; FreeDegrees[i * 6+3] = 1;
			FreeDegrees[i * 6+1] = 1; FreeDegrees[i * 6+4] = 1;
			FreeDegrees[i * 6+2] = 1; FreeDegrees[i * 6+5] = 1;
			for (int c = 0; c < _boundaryParamsSet.size(); c++)
			{
				for (int n = 0; n < _boundaryParamsSet.at(c).GetNodesCount(); n++)
				{
					if (_boundaryParamsSet.at(c).GetNode(n)-1==i)
					{
						//Force
						if (_boundaryParamsSet.at(c).GetKind() == 3)
						{
							for (int p = 0; p < 3; p++)
							{
								double param = _boundaryParamsSet.at(c).GetParam(p);
								allBoundaryForces[i*vecStride2 + p] += param / _boundaryParamsSet.at(c).GetNodesCountInFullBoundary();
								param = _boundaryParamsSet.at(c).GetParam(p + 3);
								allBoundaryForces[i*vecStride2 + vecStride + p] += param / _boundaryParamsSet.at(c).GetNodesCountInFullBoundary();
							}
						}//sealed
						else if (_boundaryParamsSet.at(c).GetKind() == 4)
						{
							for (int p = 0; p < 3; p++)
							{
								if (_boundaryParamsSet.at(c).GetParam(p) < 0)
								{
									FreeDegrees[i*6 + p] = 0;
									FreeDegrees[i*6 + p + 3] = 0;
								}
								else
								{
									FreeDegrees[i*6 + p] = 1;
									FreeDegrees[i*6 + p + 3] = 1;
								}
							}
						}
					}
				}
			}
		}

		ret = clEnqueueWriteBuffer(command_queue, BoundaryForcesBuffer, CL_FALSE, 0, vecStride2 * sizeof(cl_double) * _nElements, allBoundaryForces, 0, NULL, NULL);
		if (ret != 0)
		{
			fprintf(stderr, "Failed to write buffers.\n");
			exit(ret);
		}
		ret = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&BoundaryForcesBuffer);
		if (ret != 0)
		{
			fprintf(stderr, "Failed with args setup.\n");
			exit(ret);
		}
		ret = clEnqueueWriteBuffer(command_queue, BoundaryFixedBuffer, CL_FALSE, 0, 6*sizeof(cl_int)*_nElements,FreeDegrees, 0, NULL, NULL);
		ret = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&BoundaryFixedBuffer);
		
		cl_double* ElementStressFactorCache = _elementStressFactorCache;
		ret = clEnqueueWriteBuffer(command_queue, ElementStressFactorCacheBuffer, CL_FALSE, 0, vecStride * sizeof(cl_double)*_nElements, ElementStressFactorCache, 0, NULL, NULL);
		ret = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&ElementStressFactorCacheBuffer);

		cl_double* StressScalingFactor = _stressScalingFactors;
		ret = clEnqueueWriteBuffer(command_queue, StressScalingFactorBuffer, CL_FALSE, 0, vecStride * sizeof(cl_double), StressScalingFactor, 0, NULL, NULL);
		ret = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *)&StressScalingFactorBuffer);

		cl_double* Coordinates = _coordinates;
		ret = clEnqueueWriteBuffer(command_queue, CoordinatesBuffer, CL_TRUE, 0, 3 * sizeof(cl_double)*_nElements, Coordinates, 0, NULL, NULL);
		ret = clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *)&CoordinatesBuffer);
		/*	double _elasticModulus,			// модуль упругости (используется только для вычисления напряжений)
		double _dampingFactorAngular,	// приведенный коэффициент углового демпфирования
		double _dampingFactorLinear,	// приведенный коэффициент линейного демпфирования
		double _elasticFactorLinear,	// приведенный коэффициент линейной жесткости
		double _elasticFactorAngular	// приведенный коэффициент угловой жесткости*/
		ret = clSetKernelArg(kernel, 10, sizeof(cl_double), (void *)&_elasticModulus);
		ret = clSetKernelArg(kernel, 11, sizeof(cl_double), (void *)&_dampingFactorAngular);
		ret = clSetKernelArg(kernel, 12, sizeof(cl_double), (void *)&_dampingFactorLinear);
		ret = clSetKernelArg(kernel, 13, sizeof(cl_double), (void *)&_elasticFactorLinear);
		ret = clSetKernelArg(kernel, 14, sizeof(cl_double), (void *)&_elasticFactorAngular);
		ret = clSetKernelArg(kernel, 15, sizeof(cl_int),(void*)&_nElements);
		ret = clSetKernelArg(kernel, 16, sizeof(cl_double), (void*)&_timeStep);
	}

	void Stress::StressStrainCppIterativeSolverOpenCL::DownloadOpenCLBuffers()
	{
		//из GPU
		cl_int ret = 0;
		ret = clEnqueueReadBuffer(command_queue, DataInternalBuffer, CL_TRUE, 0, 6 * sizeof(cl_double)*vecStride*_nElements, DataInternal, 0, NULL, NULL);

		ret = clEnqueueReadBuffer(command_queue, DataRotationMatrixBuffer, CL_TRUE, 0, matStride * sizeof(cl_double)*_nElements, DataRotationMatrix, 0, NULL, NULL);

		ret = clEnqueueReadBuffer(command_queue, StressBuffer, CL_TRUE,0, 2 * vecStride * sizeof(cl_double)*_nElements, Stress,0,NULL,NULL);

		ret = clEnqueueReadBuffer(command_queue, CoordinatesBuffer,CL_TRUE,0,3*_nElements*sizeof(cl_double),_coordinates,0,NULL,NULL);
	}

	void Stress::StressStrainCppIterativeSolverOpenCL::Solve(const int nIterations)
	{

		//	Solve(nIterations);
		_iterationNumber = 0;
		_testTimer.Start(0);
		while (_iterationNumber != nIterations && _rotationSolver->IsValid())
		{
			_iterationNumber++;
			_nIteration++;
			clCalculate();
		}
		_testTimer.Stop(0);
		DownloadOpenCLBuffers();

	}

	void Stress::StressStrainCppIterativeSolverOpenCL::clCalculate()
	{
		cl_int ret = 0;
		for (size_t i = 0; i < _nElements; i += max_all_work_items)
		{
			size_t offset = i;
			size_t count = 0;
			if (offset+max_all_work_items>_nElements)
			{
				count = offset + max_all_work_items - _nElements;
			}
			else
			{
				count = max_all_work_items;
			}
			if ((ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, (size_t*)&offset, (size_t*)&count, NULL, 0, NULL, NULL)) != 0)
			{
				fprintf(stderr, "Failed to execute kernels.\n");
				exit(ret);
			}
		}
	}




	void StressStrainCppIterativeSolverOpenCL::ApplyBoundary()
	{
		//работает с типами из frm_provider
		//BCT_stresstrainBoundaryForce = 3,
		//BCT_stresstrainBoundarySealing = 4,

		double* accelerationsPointer = GetElementAcceleration(0);
		vector<BoundaryParams>::iterator it = _boundaryParamsSet.begin();

		while (it != _boundaryParamsSet.end())
		{
			switch (it->GetKind())
			{
			case 3:
				it->ApplyForceBoundary(accelerationsPointer);
				break;

			case 4:
				it->ApplySealedBoundary(accelerationsPointer);
				break;

			default:
				break;
			}
			it++;
		}
	}

	void StressStrainCppIterativeSolverOpenCL::ApplyMass()
	{
#ifdef _DEBUG
		//_controlfp(0, EM_ZERODIVIDE);
		//_control87(~_EM_ZERODIVIDE, _MCW_EM);
#endif

		//double* accelerations = GetElementAcceleration(0); // debug
		//std::cout << "Num threads = " << _numThreads << std::endl;
#ifdef OMP_SOLVE
#pragma omp parallel for num_threads(_numThreads)
#endif
		for (int elementId = 0; elementId < _nElements; elementId++)
		{
			MakeVec3(GetElementAcceleration(elementId)) /= _cellMass;
			MakeVec3(GetElementAccelerationAngular(elementId)) /= _cellInertia;
		}
	}


	void Stress::StressStrainCppIterativeSolverOpenCL::clInitKernel()
	{
		cl_program program;
		cl_int ret = 0;
		FILE *fp;
		const char fileName[] = "../../Solvers/Stress/OpenCLSolverKernels.cl";
		size_t source_size;
		char *source_str;

		try {
			fp = fopen(fileName, "r");
			if (!fp) {
				fprintf(stderr, "Failed to load kernel.\n");
				exit(ret);
			}
			source_str = (char *)malloc(MAX_SOURCE_SIZE);
			ZeroMemory(source_str, sizeof(char)*MAX_SOURCE_SIZE);
			source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
			fclose(fp);
		}
		catch (int a) {
			printf("%d", a);
		}

		program = clCreateProgramWithSource(clContext, 1, (const char **)&source_str, &source_size, &ret);
		if (ret != 0)
		{
			fprintf(stderr, "Failed to create program with source.\n");
			exit(ret);
		}
		puts("Building OpenCLSolverKernels.cl...");
		char args[200] = "-cl-std=CL2.0";
		ret = clBuildProgram(program, 1, &deviceId, args, NULL, NULL);
		if (ret != 0)
		{
			fprintf(stderr, "Failed to build kernel program.\n");
			if (ret == CL_BUILD_PROGRAM_FAILURE) {
				//Определение размера лога
				size_t log_size;
				clGetProgramBuildInfo(program, deviceId, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
				//выделение памяти под лог
				char *log = (char *)malloc(log_size);
				//считывание лога
				clGetProgramBuildInfo(program, deviceId, CL_PROGRAM_BUILD_LOG, log_size, log, NULL);

				fprintf(stderr, "%s\n", log);
			}
			exit(ret);
		}
		puts("Done");
		kernel = clCreateKernel(program, "clCalculate", &ret);
		if (ret != 0)
		{
			fprintf(stderr, "Failed to create kernel\n");
			exit(ret);
		}
		command_queue = clCreateCommandQueue(clContext, deviceId, 0, &ret);
		if (ret != 0)
		{
			fprintf(stderr, "Failed to create command queue.\n");
			exit(ret);
		}
		device_info = CL_DEVICE_MAX_WORK_ITEM_SIZES;
		size_t parameter[3] = { 0 };
		size_t fact_size = 0;
		ret = clGetDeviceInfo(deviceId, device_info, sizeof(size_t)*3, (void*)parameter,&fact_size);
		if (ret != 0)
		{
			fprintf(stderr, "Failed to getting gevice info.\n");
			exit(ret);
		}
		max_all_work_items = _nElements;
	}

	void Stress::StressStrainCppIterativeSolverOpenCL::SolveFull(const int nIterations)
	{
		Solve(nIterations);
	}

	void Stress::StressStrainCppIterativeSolverOpenCL::Solve1()
	{
	}

	void Stress::StressStrainCppIterativeSolverOpenCL::Solve2()
	{
	}

	void Stress::StressStrainCppIterativeSolverOpenCL::Solve3()
	{
	}

	void Stress::StressStrainCppIterativeSolverOpenCL::Solve4()
	{
	}

	void Stress::StressStrainCppIterativeSolverOpenCL::Solve5()
	{
	}

	void Stress::StressStrainCppIterativeSolverOpenCL::GetDisplacement(float * data)
	{

		this->StressStrainCppIterativeSolver::GetDisplacement(data);
	}

	void Stress::StressStrainCppIterativeSolverOpenCL::GetStressesByFirstTheoryOfStrength(float * data)
	{
		this->StressStrainCppIterativeSolver::GetStressesByFirstTheoryOfStrength(data);
	}

	void Stress::StressStrainCppIterativeSolverOpenCL::GetStressesByVonMises(float * data)
	{
		this->StressStrainCppIterativeSolver::GetStressesByVonMises(data);
	}

	void Stress::StressStrainCppIterativeSolverOpenCL::CalculateStrains(size_t side, double * shiftStrains, double * velocityStrains, size_t nodeId1, size_t nodeId2) const
	{
	}

	void Stress::StressStrainCppIterativeSolverOpenCL::CalculateStrainsUa(size_t side, double * shiftStrains, double * velocityStrains, size_t nodeId1, size_t nodeId2) const
	{
	}

	void Stress::StressStrainCppIterativeSolverOpenCL::GetStressesX(float * data)
	{
		this->StressStrainCppIterativeSolver::GetStressesX(data);
	}

	void Stress::StressStrainCppIterativeSolverOpenCL::GetStressesY(float * data)
	{
		this->StressStrainCppIterativeSolver::GetStressesY(data);
	}

	void Stress::StressStrainCppIterativeSolverOpenCL::GetStressesZ(float * data)
	{
		this->StressStrainCppIterativeSolver::GetStressesZ(data);
	}

	void Stress::StressStrainCppIterativeSolverOpenCL::GetStressesXY(float * data)
	{
		this->StressStrainCppIterativeSolver::GetStressesXY(data);
	}

	void Stress::StressStrainCppIterativeSolverOpenCL::GetStressesXZ(float * data)
	{
		this->StressStrainCppIterativeSolver::GetStressesXZ(data);
	}

	void Stress::StressStrainCppIterativeSolverOpenCL::GetStressesYZ(float * data)
	{
		this->StressStrainCppIterativeSolver::GetStressesYZ(data);
	}

	double3 Stress::operator-(double3 a, double3 b)
	{
		double3 result;
		result.x = a.x - b.x;
		result.y = a.y - b.y;
		result.z = a.z - b.z;
		return result;
	}

	double3 Stress::operator+(double3 a, double3 b)
	{
		double3 result;
		result.x = a.x + b.x;
		result.y = a.y + b.y;
		result.z = a.z + b.z;
		return result;
	}

	double3 Stress::operator-(double3 a)
	{
		double3 result;
		result.x = - a.x;
		result.y = - a.y;
		result.z = - a.z;
		return result;
	}
}
