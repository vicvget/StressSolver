#ifndef FMPH_H

#define FMPH_H


#include "../../../fcore/fcore.h"
#include "FElement.h"
#include "MphParamsSet/SolverParamsBase.h"
#include "MphParamsSet/SolverIntParams.h"
#include "MphParamsSet/SolverGridParams.h"

#include <vector>


using std::vector;


/**
* Класс для междисциплинарных свойств тела (набор дополнительных решателей)
*
* В заголовке лежат все наборы параметров и выбранные значения
*
* @author Getmanskiy Victor
*/
class FMph : public FElement
{
private:
	// наборы переменных для SolverParams
	vector<SolverParamsBase*> _solverParams;
public: 
	FMph(string& fileIndex):FElement(fileIndex){};
	FMph(){}
	~FMph();
	FMph(const char* path){FElement::Load(path);};
	FMph(ifstream& ifs){Load(ifs);};

	/** Добавляет и возвращает новый решатель
	* @param type - тип решателя
	* @return решатель
	*/
	SolverParamsBase* AddSolver(SolverTypes type);
	
	/** Удаляет решатель по номеру
	* @param solverId - неомер решателя
	*/
	void RemoveSolver(size_t solverId);

	/** Возвращает решатель по номеру
	* @param solverId - неомер решателя
	* @return решатель
	*/
	SolverParamsBase* GetSolver(size_t solverId);

	/** Возвращает количество решателей
	* @return количество решателей
	*/
	size_t CountSolvers() const {return _solverParams.size();};

	vector<SolverParamsBase*>& GetSolverParams() { return  _solverParams; }
	SolverParamsBase* GetSolverParams(int solverId) { return  _solverParams[solverId]; }

// См. комментарии в базовом классе
#pragma region overriden

	virtual
	bool Load
		(
			ifstream& stream
		)	override;

	void Save(ofstream& ofs) const;

#pragma endregion
};


#endif // FMPH_H