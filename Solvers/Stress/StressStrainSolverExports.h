#pragma once

#ifndef STANDALONE
#define STANDALONE
#endif

#ifdef STANDALONE
#define EXTERN
#else
#define EXTERN extern "C"
#endif


#ifndef STANDALONE
#ifdef _USRDLL
#define DLL_FUNCTION EXTERN __declspec(dllexport)
#else
#define DLL_FUNCTION EXTERN __declspec(dllimport)
#endif
#else
#define DLL_FUNCTION
#endif
#include <string>


namespace Stress
{

	/** �������������
	* @param params - ���������
	* @param links - �����
	* @param nLinks - ���������� ������
	* @param nodes - ���������� ����� (������ �����)
	* @param nNodes - ���������� �����
	* @param solverType - ��� �������� 0-iterative, 1-statics
	* @return ��������� �� ��������
	*/
	DLL_FUNCTION
	void* Init
		(
			const std::string& solverUid,
			double* params,
			int* links,
			int nLinks,
			double *nodes,
			int nNodes,
			double gridStep,
			double timeStep,
			int procType = 0,
			int numThreads = 0,
			bool useFloat = false,
			int solverType = 0
		);

	DLL_FUNCTION
	void InitIco
		(
			void* hsolver,
			const char* fileName,
			bool readIco,
			bool writeIco
		);

	DLL_FUNCTION
	bool ReadIco
		(
			void* hsolver,
			const char* fileName
		);

	DLL_FUNCTION
	void WriteIco
		(
			void* hsolver,
			const char* fileName
		);

	/** ������������� ��� ������
	* @param params - ���������
	* @param nNodes - ���������� �����
	* @param nodes - ���������� ����� (������ �����)
	*
	* @return ��������� �� ��������
	*/
	/*
	DLL_FUNCTION
	void* Init
	(
	double* params,
	float *nodes,
	int nNodes,
	double gridStep,
	double timeStep,
	int procType = 0,
	bool useFloat = false
	);
	*/

	/** ������� ��������
	* @param hsolver - ���������� ��������
	*/
	DLL_FUNCTION
	void ReleaseMemory
		(
			void* &hsolver
		);

	/** �������� ����� ��������� �� ������
	* @param hsolver - ���������� ��������
	*/
	DLL_FUNCTION
	float* GetMemoryPointer
		(
			const void* hsolver
		);

	/** �������� ������ ������
	* @param hsolver - ���������� ��������
	*/
	DLL_FUNCTION
	int GetMemorySize
		(
			const void* hsolver
		);

	/** �������� ����� ��������� �� ������
	* @param hsolver - ���������� ��������
	*/
	DLL_FUNCTION
	void* GetAccelerations(const void* hsolver);

	/** ���������� ������
	* @param hsolver - ���������� ��������
	* @param links - ������ 6*_nElements � ���������� ������ (1 - ���� �����, 0 - ��� �����)
	*/
	DLL_FUNCTION
	void AddLinks
		(
			const void* hsolver,
			int* links
		);

	/** �������� ��������� �������
	* @param hsolver - ���������� ��������
	* @param boundaryNodesIndices - ������� ��������� �����
	* @param numberOfBoundaryNodes - ���������� ��������� �����
	* @param bcKind - ��� ��������� ������ (=1 - ����, ���������� ����������� �� ������, 
	* �������� �������� �� 6 ��������� - 3 �������������� ���������� ������� ���� 
	* �� ������������ ���� � ��� ������� ������������ ����)
	* @param bcParams - ���������� ��������� �������
	* =1: 6 ���������� - ���������� ���� Fx,Fy,Fz,Mx,My,Mz
	*/
	DLL_FUNCTION
	void AddBoundary
		(
			const void *hsolver,
			int* boundaryNodesIndices, 
			int numberOfBoundaryNodes,
			int bcKind,
			double* bcParams
		);

	/** �������� ������� ���������� �������
	* @param hsolver - ���������� ��������
	* @param bcNumber - ����� ���������� ������� (� ������� ����������),
	* @param bcParamNumber - ����� ��������� � ���������� ��������� �������
	* @param bcParamValue - ����� �������� ���������
	*/
	DLL_FUNCTION
	void ChangeBoundaryParam
		(
			const void *hsolver,
			int bcNumber,
			int bcParamNumber,
			double bcParamValue
		);

	/** �������� �������� �� ������ ���������� ������� � �� �������
	* @param hsolver - ���������� ��������
	* @param bcNumber - ����� ���������� �������
	* @param bcParamNumber - ����� ���������
	*/
	DLL_FUNCTION
	double GetBoundaryParam
		(
			const void *hsolver,
			int bcNumber,
			int bcParamNumber
		);

	/**
	* �������� �����
	* @param hsolver - ���������� ��������
	* @param scaleFactor - ���������� �����������
	*/
	DLL_FUNCTION
	float UpdateBuffer
		(
			const void* hsolver,
			double scaleFactor = 1.0
		);

	DLL_FUNCTION
		void InitialSolve
		(
		void* hsolver
		);

	/** ������
	* @param hsolver - ���������� ��������
	* @param nIterations - ���������� ��������
	*/
	DLL_FUNCTION
	void Solve
		(
			void* hsolver,
			int nIterations
		);

	/**
	* ������ ������ ������ ������ �����-�����
	* @param hsolver - ���������� ��������
	*/
	DLL_FUNCTION
	void Solve1
		(
			void* hsolver
		);

	/**
	* ������ ������ ������ ������ �����-�����
	* @param hsolver - ���������� ��������
	*/
	DLL_FUNCTION
	void Solve2
		(
			void* hsolver
		);

	/**
	* ������ ������� ������ ������ �����-�����
	* @param hsolver - ���������� ��������
	*/
	DLL_FUNCTION
	void Solve3
		(
			void* hsolver
		);

	/**
	* ������ ��������� ������ ������ �����-�����
	* @param hsolver - ���������� ��������
	*/
	DLL_FUNCTION
	void Solve4
		(
			void* hsolver
		);

	/**
	* ������ ����� ������ ������ �����-�����
	* @param hsolver - ���������� ��������
	*/
	DLL_FUNCTION
	void Solve5
		(
			void* hsolver
		);

	DLL_FUNCTION
	double GetData
		(
			const void* hsolver,
			int nNode,
			int dir,
			int type
		);

	//DLL_FUNCTION
	//void SolveSystemOfLinearEquationsForStiffness();

}

