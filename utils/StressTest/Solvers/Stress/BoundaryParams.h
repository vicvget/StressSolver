#ifndef BoundaryParamsH
#define BoundaryParamsH


#include "AuxiliaryStressStuff.h"

#include <string>


class BoundaryParams
{
public:

	BoundaryParams
		(
			const BoundaryParams& bp
		);
	
	BoundaryParams
		(
			const int kind,
			const double* params,
			const int* nodes,
			const int nodesCount,
			size_t stride = 3
		);

	~BoundaryParams();

	int GetKind() const;

	void GetParams
		(
			const double*& params,
			int& paramsCount
		)	const;

	double GetParam
		(
			int paramIndex
		)	const;

	void GetNodes
		(
			const int*& nodes,
			int& nodesCount
		)	const;

	int GetNodesCount() const;

	int GetNode
		(
			int nodeIndex
		)	const;

	void SetParam
		(
			double param,
			int paramIndex
		);


	void ApplyForceBoundary
		(
			double* data
		)	const;

	void CorrectForceBoundary
		(
			double* accelerations,
			const double* velocities,
			const double* shifts,
			const float* nodes,
			double elasticModulus,
			double dampingFactor
		)	const;

	void ApplySealedBoundary
		(
			double* data
		)	const;

	void FormIsSealedFlagsList
		(
			IsSealedFlagsList& isSealedFlagsList
		)	const;

	void SetNodesCountInFullBoundary
		(
			int nodesCountInFullBoundary
		);

	int GetNodesCountInFullBoundary() const;


private:

	int _kind; //BCT_StressStrainBoundaryForce = 3, BCT_StressStrainBoundarySealing = 4
	double* _params;
	int* _nodes;
	int _nodesCount;
	int _nodesCountInFullBoundary;

	size_t _stride; // расстояние между векторами неизвестных
	size_t _varStride; // расстояние между векторами неизвестных


	int GetParamsCount() const;

	static
	void CheckIndex
		(
			int index,
			int count,
			const std::string& message
		);

};

#endif