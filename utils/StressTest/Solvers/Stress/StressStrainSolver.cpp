#include "StressStrainSolver.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include "Common.h"
#include "../../AdditionalModules/fmath/Vector3.h"

using std::ifstream;
using std::ofstream;

StressStrainSolver::StressStrainSolver
	(
		int const nNodes,
		int stride = 4
	)
	:
		_dataSize(nNodes * stride),
		_nElements(nNodes),
		_readIco(false),
		_writeIco(false),
		vecStride(stride),
		vecStride2(stride*2),
		matStride(stride*3)
{

const size_t vecStride = stride;				// смещение векторов
const size_t vecStride2 = vecStride * 2;
const size_t matStride = vecStride * 3; // смещения матриц поворота (3x3, 3x4)

	_data = new float[_dataSize];

	size_t dataInternalSize		= sizeof(double)*nNodes * vecStride2 * 3;
	size_t dataRotationMtxSize	= sizeof(double)*nNodes * matStride;
	size_t stressVectorSize		= sizeof(double)*nNodes * vecStride2;

	// Additional alignment for KNC
	size_t 
	delta = dataInternalSize % 8;
	if (delta) dataInternalSize += (8 - delta);
	delta = dataRotationMtxSize % 8;
	if (delta) dataRotationMtxSize += (8 - delta);
	delta = stressVectorSize % 8;
	if (delta) stressVectorSize += (8 - delta);

	// массивы для внутреннего хранения данных
	//TODO: refactor to 3 arrays
	_dataInternal		= (double*)aligned_alloc(dataInternalSize, ALIGNMENT); // поступательные и вращательные 3 - перемещения, скорости, ускорения
	_dataRotationMtx	= (double*)aligned_alloc(dataRotationMtxSize, ALIGNMENT); // TODO: выравнивание матриц поворота
	_stress				= (double*)aligned_alloc(stressVectorSize, ALIGNMENT);
	_buffer				= (double*)aligned_alloc(8*5, ALIGNMENT);

	// обнуление массивов
	memset(_stress, 0, stressVectorSize);
	memset(_dataRotationMtx, 0, dataRotationMtxSize);
	memset(_dataInternal, 0, dataInternalSize);

	// единичные матрицы поворота
	for(size_t i = 0; i < _nElements; i++)
	{
		_dataRotationMtx[i * matStride] = 1.;
		_dataRotationMtx[i * matStride + vecStride + 1] = 1.;
		_dataRotationMtx[i * matStride + vecStride2 + 2] = 1.;
	}
}

void StressStrainSolver::SetZeroVelocities()
{
	memset(_dataInternal + _nElements * vecStride2, 0, sizeof(double)*_nElements * vecStride2);
}


StressStrainSolver::~StressStrainSolver()
{
	if (_data != NULL)
	{
		delete [] _data;
		_data = NULL;
	}

	aligned_free(_dataInternal);
	aligned_free(_dataRotationMtx);
	aligned_free(_stress);
	aligned_free(_buffer);
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
		ofs.write((char*)&_nElements, sizeof(int));
		ofs.write((char*)dataPointer, _nElements*sizeof(double)*6);
		ofs.write((char*)mtxPointer, _nElements*sizeof(double)*9);
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
	return _dataInternal + _nElements * vecStride2 * (int)dataType;
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

	return _dataInternal[_nElements * vecStride2 * type + nNode * vecStride2 + dir];
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
	//GetStressesByFirstTheoryOfStrength(data);
}


double* StressStrainSolver::GetElementGridCoordinates(size_t elementId) const
{
	return _elements + elementId * 3;
}

//virtual 
double* StressStrainSolver::GetElementStress(size_t elementId) const
{
	return _stress + (elementId * vecStride2);
}

double* StressStrainSolver::GetElementStressAngular(size_t elementId) const
{
	return _stress + (elementId * vecStride2 + vecStride);
}


//virtual 
double* StressStrainSolver::GetElementShift(size_t elementId) const
{
	return _dataInternal + (elementId * vecStride2);
}

//virtual 
double StressStrainSolver::GetElementDisplacement(size_t elementId) const
{
	return (MathHelpers::MakeVec3(GetElementShift(elementId)) - MathHelpers::MakeVec3(GetElementGridCoordinates(elementId))).Magnitude();
}


//virtual 
double* StressStrainSolver::GetElementVelocity(size_t elementId) const
{
	return _dataInternal + (_nElements * vecStride2 + elementId * vecStride2);
}

//virtual 
double* StressStrainSolver::GetElementAcceleration(size_t elementId) const
{
	return _dataInternal + (_nElements * vecStride2 * 2 + elementId * vecStride2);
}

//virtual 
double* StressStrainSolver::GetElementShiftAngular(size_t elementId) const
{
	return _dataInternal + (elementId * vecStride2 + vecStride);
}

//virtual 
double* StressStrainSolver::GetElementVelocityAngular(size_t elementId) const
{
	return _dataInternal + (_nElements * vecStride2 + elementId * vecStride2 + vecStride);
}

//virtual 
double* StressStrainSolver::GetElementAccelerationAngular(size_t elementId) const
{
	return _dataInternal + (_nElements * vecStride2 * 2 + elementId * vecStride2 + vecStride);
}

//virtual 
double* StressStrainSolver::GetRotationMatrix(size_t elementId) const
{
	return _dataRotationMtx + (elementId * matStride);
}

void StressStrainSolver::SetUid(const string& uid)
{
	_uid = uid;
}
