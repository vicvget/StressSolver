#ifndef SOLVER_INT_PARAMS_H

#define SOLVER_INT_PARAMS_H


#include "SolverParamsBase.h"

/** Набор параметров интегрирования
* (наборы хранятся в заголовке, связка по Id)
*
* @author Getmanskiy Victor
*/
class SolverIntParams: public ParamsSetBase
{
protected:
	// шаг интегрирования
	double _timeStep;
	// количество кадров дял записи в файл результатов
	size_t _numberOfFrames;
public:
	SolverIntParams(string &name);
	SolverIntParams();
	~SolverIntParams() {};

	double GetTimeStep() const {return _timeStep;}
	void SetTimeStep(double timeStep) {_timeStep = timeStep;}

	size_t GetNumberOfFrames() const {return _numberOfFrames;}
	void SetNumberOfFrames(size_t numberOfFrames) {_numberOfFrames = numberOfFrames;}



#pragma region accessors
	double TimeStep() const {return _timeStep;}
	void TimeStep(double val){_timeStep = val;}
	size_t NumberOfFrames() const {return _numberOfFrames;}
	void NumberOfFrames(int val){_numberOfFrames = val;}
#pragma endregion

//TODO:
#pragma region overriden	
	void Init();
	void Save(ofstream& ofs) const;
	void Load(ifstream& ifs);
#pragma endregion
};


#endif // SOLVER_INT_PARAMS_H