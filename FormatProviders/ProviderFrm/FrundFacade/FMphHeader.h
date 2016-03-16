#ifndef FMPH_HEADER_H

#define FMPH_HEADER_H


#include "BasicTypes.h"
#include "FVariablesGroupsSet.h"
#include "FElement.h"
#include "../../../fcore/fcore.h"
#include "MphParamsSet/BcParams.h"
#include "MphParamsSet/SolverIntParams.h"
#include "MphParamsSet/SolverGridParams.h"
#include "MphParamsSet/SolverSpecialParams.h"

#include <vector>


using std::vector;


/**
* Класс для считывания заголовка междисциплинарной модели ФРУНД
* В отличие от классической модели ФРУНДа
* здесь типы:
* Temperature = 0
* StressStrain = 1
*
*
* В заголовке лежат все наборы параметров и выбранные значения    
*
* @author Getmanskiy Victor
*/
class FMphHeader : public FElement
{
private:
	// наборы переменных для SolverParams
	vector<SolverIntParams> _intParamsSets;
	vector<SolverGridParams> _gridParamsSets;
	vector<SolverSpecialParams> _specialParamsSets;
	vector<BcParams> _bcParamsSets;
	void Init();
public: 
	FMphHeader();
	FMphHeader(const char* path);
	FMphHeader(ifstream& ifs){Load(ifs);};

	size_t CountIntParamsSets() const {return _intParamsSets.size();}
	size_t CountGridParamsSets() const  {return _gridParamsSets.size();}
	size_t CountSpecialParamsSets() const  {return _specialParamsSets.size();}
	size_t CountBcParamsSets() const  {return _bcParamsSets.size();}

	/** Добавить новый набор параметров интегрирования
	* @return набор параметров интегрирования
	*/
	SolverIntParams& AddNewIntParams();
	
	/** Добавить параметры сетки
	* @return набор параметров сетки
	*/
	SolverGridParams& AddNewGridParams();
	
	/** Добавить новый набор специальных параметров
	* @param type - тип решателя
	* @return набор специальных параметров
	*/
	SolverSpecialParams& AddNewSpecialParams(SolverTypes type);

	/** Добавить граничные условия
	* @param type - тип граничного условия
	* @return набор параметров граничных условий
	*/
	BcParams& AddNewBcParams(BcTypes type);

	/** Получить параметры интегрирования
	* @param id - номер набора параметров
	* @return набор параметров интегрирования
	*/
	SolverIntParams& GetIntParams(size_t id);
	const SolverIntParams& GetIntParams(size_t id) const;

	/** Получить параметры сетки
	* @param id - номер набора параметров
	* @return набор параметров сетки
	*/	
	SolverGridParams& GetGridParams(size_t id);
	const SolverGridParams& GetGridParams(size_t id) const;
	
	/** Получить специальные параметры
	* @param id - номер набора параметров
	* @return набор специальных параметров
	*/
	SolverSpecialParams& GetSpecialParams(size_t id);
	const SolverSpecialParams& GetSpecialParams(size_t id) const;

	/** Получить параметры граничных условий
	* @param id - номер набора параметров
	* @return набор параметров граничных условий
	*/	
	BcParams& GetBcParams(size_t id);
	const BcParams& GetBcParams(size_t id) const;

#pragma region accessors
	vector<BcParams> BcParamsSets() const { return _bcParamsSets; }
	void BcParamsSets(vector<BcParams> val) { _bcParamsSets = val; }

	vector<SolverGridParams> GridParamsSets() const { return _gridParamsSets; }
	void GridParamsSets(vector<SolverGridParams> val) { _gridParamsSets = val; }

	vector<SolverIntParams> IntParamsSets() const { return _intParamsSets; }
	void IntParamsSets(vector<SolverIntParams> val) { _intParamsSets = val; }

	vector<SolverSpecialParams> SpecialParamsSets() const { return _specialParamsSets; }
	void SpecialParamsSets(vector<SolverSpecialParams> val) { _specialParamsSets = val; }
#pragma endregion

// См. комментарии в базовом классе
#pragma region overriden

	virtual
	bool Load
		(
			ifstream& stream
		)	override;

	virtual
	void Save(ofstream& ofs) const override;

#pragma endregion

};


#endif // FMPH_HEADER_H