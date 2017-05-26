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
	/** Решатель НДС методом динамики систем тел на регулярной ортогональной сетке
	*/
	class StressStrainCppSolver: public StressStrainSolver
	{
	public:

		void PrintTime();
		// создает объект с заданными параметрами
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

		// признак первой итерации
		bool _isFirstSolution;

		// признак первой итерации
		vector<BoundaryParams> _boundaryParamsSet;

		// кэш-массивы, используемые тольк ов solve
		double *_initX;
		double *_initDX;
		double *_hDDX1;
		double *_hDDX2;
		double *_hDDX3;

		// кэш-массивы
		double* _varX;		// смещения до неизвестных
		double* _varDX;		// смещения до производных
		double* _varDDX;	// смещения до 2 производных

		double* _radiusVectors;
		int* _linkedElements;

		size_t _nVariables;				// количество неизвестных
		size_t _nIteration;				// номер итерации
		int _numThreads;				// число параллельных потоков OpenMP

		//double _time;					// время
		double _timeStep;				// шаг по времени
		double _timeStep2;				// шаг по времени /2
		double _timeStep4;				// шаг по времени /4

		double _gridStep;				// шаг сетки
		double _gridStep2;				// шаг сетки в квадрате
		double _gridStep3;				// шаг сетки в кубе

		double _elasticModulus;			// модуль упругости (используется только для вычисления напряжений)

		double _dampingFactorAngular;	// приведенный коэффициент углового демпфирования
		double _dampingFactorLinear;	// приведенный коэффициент линейного демпфирования
		double _elasticFactorLinear;	// приведенный коэффициент линейной жесткости
		double _elasticFactorAngular;	// приведенный коэффициент угловой жесткости
		bool _isStiffnessOverriden;		// признак перегрузки реальной жесткости и демпфирования
		bool _isInertiaOverriden;		// признак перегрузки реальной массы и моментов инерции

		// Масштабные коэффициенты для напряжений по X,Y,Z
		double _stressScalingFactors[3];
		double _puassonFactor;


		//double _density;				// плотность материала
		double _cellMass;				// масса элемента
		double _cellInertia;			// момент элемента
		double _stiffScale;				// масштабирование
		double _poissonRatio;			// коэффициент Пуассона
		double _lameMatrix[6][6];		// матрица перехода от деформациям к напряжениям

		FTimer _testTimer;			// замеры времени

		// построение матрицы по параметрам Ламе
		void FindStressStrainMatrix();

		// Возвращает радиус-вектор крепления упругой связи 
		// для регулярной сетки в зависимости от грани элемента
		// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
		// @param side - номер грани
		// @returns радиус-вектор
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