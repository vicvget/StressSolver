#ifndef SOLVER_PARAMS_BASE_H

#define SOLVER_PARAMS_BASE_H


#include "ParamsSetBase.h"
#include "SolverSpecialParams.h"

/** Набор параметров, общих для всех междисциплинарных решателей
*
* @author Getmanskiy Victor
*/
class SolverParamsBase: public ParamsSetBase
{
protected:
	// тип решателя
	SolverTypes _solverType;	
	// номер набора параметров интегрирования из заголовка
	size_t _intParamsId;
	// номер набора параметров сетки из заголовка
	size_t _gridParamsId;
	// номер набора специальных параметров решателя
	size_t _specialParamsId;

	
	bool _readIco;  // признак чтения ico
	bool _writeIco; // признак записи ico

	// признак активации решателя
	bool _isEnabled;

	// прозрачность тела
	double _solidTransparency;
	
public:
	void Disable() {_isEnabled = false;}
	void Enable() {_isEnabled = true;}

	SolverParamsBase(string &name, SolverTypes solverType);
	SolverParamsBase(SolverTypes solverType);
	~SolverParamsBase() {};

#pragma region properties
	bool IsReadIco() const {return _readIco;}
	bool IsWriteIco() const {return _writeIco;}
	void SetIsReadIco(bool readIco) {_readIco=readIco;}
	void SetIsWriteIco(bool writeIco) {_writeIco=writeIco;}
	SolverTypes GetSolverType() const {return _solverType;}
	void SetSolverType(SolverTypes solverType) {_solverType = solverType;}
	size_t GetIntParamsId() const {return _intParamsId;}
	void SetIntParamsId(int intParamsId) {_intParamsId = intParamsId;}
	size_t GetGridParamsId() const {return _gridParamsId;}
	void SetGridParamsId(int gridParamsId) {_gridParamsId = gridParamsId;}
	bool GetIsEnabled() const {return _isEnabled;}
	void SetIsEnabled(bool isEnabled) {_isEnabled = isEnabled;}
	void SetSpecialParamsId(int specialParamsId) {_specialParamsId = specialParamsId;}
	size_t GetSpecialParamsId() const {return _specialParamsId;}
	double SolidTransparency() const { return _solidTransparency; }
	void SolidTransparency(double val) { _solidTransparency = val; }
#pragma endregion 

// TODO:
#pragma region overriden
	virtual void Init();
	virtual void Save(ofstream& ofs) const;
	virtual void Load(ifstream& ifs);
#pragma endregion
};


#endif // SOLVER_PARAMS_BASE_H