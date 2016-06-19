#include "BoundaryParams.h"
#include <cstring>

using std::string;


BoundaryParams::BoundaryParams
(
const BoundaryParams& bp
)
{
	_varStride = bp._varStride;
	_stride = bp._stride;
	_kind = bp._kind;
	_nodesCount = bp._nodesCount;
	_nodesCountInFullBoundary = bp._nodesCountInFullBoundary;

	int nParams = 6;

	if (_kind == 3)
	{
		nParams++;
	}
	_params = new double[nParams];
	_nodes = new int[_nodesCount];
	memcpy(_params, bp._params, sizeof(double) * nParams);
	memcpy(_nodes, bp._nodes, sizeof(int) * _nodesCount);
}

BoundaryParams::BoundaryParams
(
const int kind,
const double* params,
const int* nodes,
const int nodesCount,
size_t stride
)
:
_kind(kind),
_nodesCount(nodesCount),
_nodesCountInFullBoundary(nodesCount),
_stride(stride),
_varStride(stride * 2)
{
	_params = NULL;
	_nodes = NULL;

	int nParams = 6;

	if (kind == 3)
	{
		nParams++;
	}

	if (kind == 3 || kind == 4)
	{
		_params = new double[nParams];
		_nodes = new int[_nodesCount];
		memcpy(_params, params, sizeof(double) * nParams);
		memcpy(_nodes, nodes, sizeof(int) * _nodesCount);
	}
}

BoundaryParams::~BoundaryParams()
{
	if (_params != NULL)
	{
		delete[] _params;
		_params = NULL;
	}
	if (_nodes != NULL)
	{
		delete[] _nodes;
		_nodes = NULL;
	}
}

int BoundaryParams::GetKind() const
{
	return _kind;
}

void BoundaryParams::GetParams
(
const double*& params,
int& paramsCount
)	const
{
	params = _params;
	paramsCount = GetParamsCount();
}

double BoundaryParams::GetParam
(
int paramIndex
)	const
{
	CheckIndex
		(
		paramIndex,
		GetParamsCount(),
		"Incorrect paramIndex value: paramIndex < 0 or paramIndex >= _paramsCount (BoundaryParams)"
		);

	return _params[paramIndex];
}

void BoundaryParams::GetNodes
(
const int*& nodes,
int& nodesCount
)	const
{
	nodes = _nodes;
	nodesCount = _nodesCount;
}

int BoundaryParams::GetNode
(
int nodeIndex
)	const
{
	return _nodes[nodeIndex];
}

int BoundaryParams::GetNodesCount() const
{
	return _nodesCount;
}

void BoundaryParams::SetParam
(
double param,
int paramIndex
)
{
	CheckIndex
		(
		paramIndex,
		GetParamsCount(),
		"Incorrect paramIndex value: paramIndex < 0 or paramIndex >= _paramsCount (BoundaryParams)"
		);
	_params[paramIndex] = param;
}

void BoundaryParams::ApplyForceBoundary
(
double* data
)	const
{
	for (int j = 0; j < _nodesCount; j++)
	{
		for (int i = 0; i < 3; i++)
		{
			int nodeId = _nodes[j] - 1;
			int nodeOffset = static_cast<int>(_varStride)* nodeId;

			data[nodeOffset + i] += _params[i] / _nodesCountInFullBoundary; // равномерное распределение
			data[nodeOffset + i + _stride] += _params[i + 3] / _nodesCountInFullBoundary; // равномерное распределение
		}
	}
}

void BoundaryParams::CorrectForceBoundary
(
double* accelerations,
const double* velocities,
const double* shifts,
const float* nodes,
double elasticModulus,
double dampingFactor
)	const
{
	/// TODO: ЗАЧЕМ ДЕЛИТЬ НА 10 ????
	double conserv = elasticModulus / 10;
	double dissip = dampingFactor / 10;


	for (int j = 0; j < _nodesCount; j++)
	{
		int nodeId = _nodes[j] - 1;
		int nodeOffset = static_cast<int>(_varStride)* nodeId;

		for (int i = 0; i < 3; i++)
		{
			double nodePosisiton = nodes[3 * nodeId + i];

			accelerations[nodeOffset + i] -=
				conserv * (shifts[nodeOffset + i] - nodePosisiton) +
				dissip * velocities[nodeOffset + i];
			accelerations[nodeOffset + i + _stride] -=
				conserv * (shifts[nodeOffset + i + _stride] - nodePosisiton) +
				dissip * velocities[nodeOffset + i + _stride];
		}
	}
}

void BoundaryParams::ApplySealedBoundary
(
double* data
)	const
{
	for (int j = 0; j < _nodesCount; j++)
	{
		int nodeOffset = static_cast<int>(_varStride)* (_nodes[j] - 1);

		for (int i = 0; i < 3; i++)
		{
			if (_params[i] < 0)
			{
				data[nodeOffset + i] = 0.0;
			}
			if (_params[i + 3] < 0)
			{
				data[nodeOffset + i + _stride] = 0.0;
			}
		}
	}
}

void BoundaryParams::FormIsSealedFlagsList
(
IsSealedFlagsList& isSealedFlagsList
)	const
{
	for (int boundaryNodeIndex = 0; boundaryNodeIndex < _nodesCount; ++boundaryNodeIndex)
	{
		int nodeIndex = _nodes[boundaryNodeIndex] - 1;
		IsSealedFlags& isSealedFlags = isSealedFlagsList[nodeIndex];

		for (int dofIndex = 0; dofIndex < 6; ++dofIndex)
		{
			if (_params[dofIndex] < 0)
			{
				isSealedFlags[dofIndex] = true;
			}
		}
	}
}

void BoundaryParams::SetNodesCountInFullBoundary(int nodesCountInFullBoundary)
{
	_nodesCountInFullBoundary = nodesCountInFullBoundary;
}

int BoundaryParams::GetNodesCountInFullBoundary() const
{
	return _nodesCountInFullBoundary;
}

int BoundaryParams::GetParamsCount() const
{
	int nParams = 6;

	if (_kind == 3)
	{
		nParams++;
	}

	return nParams;
}

// static
void BoundaryParams::CheckIndex
(
int index,
int count,
const string& message
)
{
	if
		(
		(index < 0) ||
		(index >= count)
		)
	{
		throw message;
	}
}