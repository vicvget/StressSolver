#include "RegularMemory.h"
#include <math.h>

StressStrainSolver::StressStrainSolver
	(
		int const nNodes
	)
	:
		_dataSize(nNodes * 4),
		_nNodes(nNodes)
{
	_data = new float[_dataSize];
	_dataInternal = new double[nNodes * 6 * 3];
	_stress = new double[nNodes * 6];
	memset(_stress, 0, _nNodes * 6 * sizeof(double));
}

StressStrainSolver::~StressStrainSolver()
{
	if (_data != NULL)
	{
		delete [] _data;
		_data = NULL;
	}
	if (_dataInternal != NULL)
	{
		delete [] _dataInternal;
		_dataInternal = NULL;
	}
	if (_stress != NULL)
	{
		delete[] _stress;
		_stress = NULL;
	}
}

float* StressStrainSolver::GetMemoryPointer() const
{
	return _data;
}

int StressStrainSolver::GetMemorySize() const
{
	return _dataSize * sizeof(float);
}

double* StressStrainSolver::GetDataInternal(DataType dataType) const
{
	return _dataInternal+_nNodes * 6 * (int)dataType;
}

double StressStrainSolver::GetData
	(
		int nNode,
		int dir,
		int type
	)
{
	//std::cout << _dataInternal[_nNodes * 6 * type + nNode + dir] << std::endl;

	return _dataInternal[_nNodes * 6 * type + nNode * 6 + dir];
}

/** Получить скалярный параметр
* @param data - массив для записи скалярного параметра
*/
// virtual
void StressStrainSolver::GetScalarParameter
	(
		float* data
	)
{
	GetStressesByVonMises(data);
}