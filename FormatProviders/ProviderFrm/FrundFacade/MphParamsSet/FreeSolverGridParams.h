#ifndef FreeSolverGridParamsH
#define FreeSolverGridParamsH

#ifdef USE_PROVIDER_FRM_OLD
#	include "../../../../Fcore/FrundFacade/MphParamsSet/ParamsSetBase.h"
#	include "../../../../Fcore/FrundFacade/MphParamsSet/SolverGridParams.h"
#	include "../../../../Fcore/FrundFacade/MphParamsSet/BcMapper.h"
#	include "../../../../Fcore/FrundFacade/BasicTypes.h"
#	include "../../../../Fcore/FrundFacade/FGeometryPoint.h"
#else
#	include "ParamsSetBase.h"
#	include "SolverGridParams.h"
#	include "BcMapper.h"
#	include "../BasicTypes.h"
#	include "../FGeometryPoint.h"
#endif

#include <vector>
#include <ios>
using std::vector;
using std::string;

/** Набор параметров сетки,
* общих для всех свободных междисциплинарных решателей
*
* @author Gromov Eugene
*/
class FreeSolverGridParams: public ParamsSetBase
{
private:
	// шаг сетки - уточнить, нужен или нет
	double _gridStep;
	// у свободного решателя только один маппер
	BcMapper *_bcMapper;

public:
	FreeSolverGridParams();
	FreeSolverGridParams(const string &name);
	FreeSolverGridParams(SolverGridParams &solverGridParams);
	FreeSolverGridParams(const FreeSolverGridParams &freeSolverGridParams);
	~FreeSolverGridParams() { /*delete _bcMapper;*/ };

	BcMapper* GetModifiableMapper() { return _bcMapper; }

#pragma region accessors

	double GridStep() const {return _gridStep;}
	void GridStep(double val){_gridStep = val;}
	BcMapper GetBcMapper() const { return *_bcMapper; }
	void SetBcMapper(BcMapper val) { delete _bcMapper; _bcMapper = new BcMapper(val); }

#pragma endregion

	void SaveBinary(ofstream &ofs);
	void LoadBinary(ifstream &ifs);

#pragma region overriden	

	void Init();
	void Save(ofstream& ofs) const;
	void Load(ifstream& ifs);

#pragma endregion
};

#endif