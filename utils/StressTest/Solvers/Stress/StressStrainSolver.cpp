#include "StressStrainSolver.h"
#include <cmath>
#include <cstdlib>
#include <fstream>

using std::ifstream;
using std::ofstream;

StressStrainSolver::StressStrainSolver
	(
		int const nNodes
	)
	:
		_dataSize(nNodes * 4),
		_nNodes(nNodes),
		_readIco(false),
		_writeIco(false)
{
	_data = new float[_dataSize];

	size_t dataInternalSize		= sizeof(double)*nNodes * vecStride2 * 3;
	size_t dataRotationMtxSize	= sizeof(double)*nNodes * matStride;
	size_t stressVectorSize		= sizeof(double)*nNodes * vecStride2;

	// массивы для внутреннего хранения данных
	_dataInternal		= (double*)_aligned_malloc(dataInternalSize, alignment); // поступательные и вращательные 3 - перемещения, скорости, ускорения
	_dataRotationMtx	= (double*)_aligned_malloc(dataRotationMtxSize, alignment); // TODO: выравнивание матриц поворота
	_stress				= (double*)_aligned_malloc(stressVectorSize, alignment);
	
	// обнуление массивов
	memset(_stress, 0, stressVectorSize);
	memset(_dataRotationMtx, 0, dataRotationMtxSize);
	memset(_dataInternal, 0, dataInternalSize);

	// единичные матрицы поворота
	for(size_t i = 0; i < _nNodes; i++)
	{
		_dataRotationMtx[i * matStride] = 1.;
		_dataRotationMtx[i * matStride + vecStride + 1] = 1.;
		_dataRotationMtx[i * matStride + vecStride2 + 2] = 1.;
	}
}

void StressStrainSolver::SetZeroVelocities()
{
	memset(_dataInternal + _nNodes * vecStride2, 0, sizeof(double)*_nNodes * vecStride2);
}


StressStrainSolver::~StressStrainSolver()
{
	if (_data != NULL)
	{
		delete [] _data;
		_data = NULL;
	}

	_aligned_free(_dataInternal);
	_aligned_free(_dataRotationMtx);
	_aligned_free(_stress);
}

void StressStrainSolver::InitIco ( const string& fileName,
		bool readIco,
		bool writeIco,
		int nWriteIteration)
{
	_fileName = fileName;
	_readIco = readIco;
	_writeIco = writeIco;
	_nWriteIteration = nWriteIteration;
}

bool StressStrainSolver::ReadIco(const char* fileName)
{
	ifstream ifs(fileName, std::ios_base::binary);
	int nNodes;
	double* dataPointer = GetDataInternal(DataType::DT_Shifts);
	double* mtxPointer = GetDataRotaionMtx();

	if(ifs.is_open())
	{
		ifs.read((char*)&nNodes, sizeof(int));
		ifs.read((char*)dataPointer, nNodes*sizeof(double)*6);
		ifs.read((char*)mtxPointer, nNodes*sizeof(double)*9);
		ifs.close();
		return true;
	}
	return false;
}

void StressStrainSolver::WriteIco(const char* fileName) const
{
	ofstream ofs(fileName, std::ios_base::binary);
	double* dataPointer = GetDataInternal(DataType::DT_Shifts);
	double* mtxPointer = GetDataRotaionMtx();

	if(ofs.is_open())
	{
		ofs.write((char*)&_nNodes, sizeof(int));
		ofs.write((char*)dataPointer, _nNodes*sizeof(double)*6);
		ofs.write((char*)mtxPointer, _nNodes*sizeof(double)*9);
		ofs.close();
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
	return _dataInternal + _nNodes * vecStride2 * (int)dataType;
}

double* StressStrainSolver::GetDataRotaionMtx() const
{
	return _dataRotationMtx;
}


double StressStrainSolver::GetData
	(
		int nNode,
		int dir,
		int type
	)
{
	//std::cout << _dataInternal[_nElements * 6 * type + nNode + dir] << std::endl;

	return _dataInternal[_nNodes * vecStride2 * type + nNode * vecStride2 + dir];
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
	//GetStressesByVonMises(data);
	GetStressesByFirstTheoryOfStrength(data);
}