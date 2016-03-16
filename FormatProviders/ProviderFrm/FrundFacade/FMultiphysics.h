#ifndef FMULTIPHYSICS_H

#define FMULTIPHYSICS_H


#include "BasicTypes.h"

#include <fstream>
#include <string>
#include <vector>
#include <set>


using std::set;


// Типы междисциплинарных решателей
enum MultiphysicsSolverTypes
{
	MPHHeatTransfer = 0,
	MPHStressDeformation = 1,
	FComplexHeatTransferParams = 2
};

struct FComplexHeatTransferParams
{
	// время расчета
	double time;
	// шаг по времени
	double timeStep;
	// шаг записи
	double recordStep;
};

/** Хранение параметров для теплового расчета
*/
struct FHeatTransferParams
{
	// время расчета
	double time;
	// шаг по времени
	double timeStep;
	// шаг записи
	double recordStep;
	// имя файла с сеткой
	char fileGrid[STRLIMIT];
	
	int _numberOfParams;
	union
	{
		struct
		{
			// теплопроводность
			double TPR;
			// теплоемкость
			double TPL;
			// плотность
			double PLOT;

		};
		double params[5];
	}p;

	FHeatTransferParams()
	{
		_numberOfParams = 5;
	}
	
	/** Загрузка из текстового потока
	* @param ifs - входной текстовый поток
	*/
	void Load(ifstream& ifs);
	
	/** Сохранение в текстовый поток
	* @param ofs - выходной текстовый поток
	*/
	void Save(ofstream& ofs) const;
};

/** Хранение параметров для расчета напряженно-деформированного состояния
*/
struct FStressDeformationParams
{
	// время расчета
	double time;
	// шаг по времени
	double timeStep;
	// шаг записи
	double recordStep;
	// имя файла с сеткой
	char fileGrid[STRLIMIT];
	// коэффициент жесткости
	double K;
	// модуль упругости
	double E;
	// точки сетки, для связи с MBS моделью
	//vector<FGeometryPoint> points;
	// идентификаторы соединительных элементов
	//vector<int> loadIds;

	/** Загрузка из текстового потока
	* @param ifs - входной текстовый поток
	*/
	void Load(ifstream& ifs);
	
	/** Сохранение в текстовый поток
	* @param ofs - выходной текстовый поток
	*/
	void Save(ofstream& ofs) const;
};


/** Класс для работы с mph-Файлами
*
* @author Getmanskiy Victor
*/
class FMultiphysics
{
	// множество типов решателей 0 - Thermal, 1 - StressStrain
	set<int> _solverTypes;
	// параметры теплового расчета
	FHeatTransferParams _heatTransferParams;
	// параметры расчета НДС
	FStressDeformationParams _stressDeformationParams;

public:

	FMultiphysics()
		:
			// TODO: сделать конструктор по умолчанию в классе FStressDeformationParams и убрать эту строчку
			_stressDeformationParams()
	{
	}
	
	FMultiphysics
		(
			ifstream& ifs
		);
	
	/** Определение наличия решателя по его типу
	* @param solverType - тип решателя
	* @return true, если есть
	*/
	bool HasSolver(int solverType) {return _solverTypes.find(solverType) != _solverTypes.end();}

	/** Возвращает указатель на параметры теплового расчета
	* @return nullptr,если нет теплового решателя
	*/
	FHeatTransferParams* GetHeatTransferParams() {return HasSolver(MPHHeatTransfer) ? &_heatTransferParams : nullptr;}
	void SetHeatTransferParams(FHeatTransferParams* val){_heatTransferParams = *val;}
	/** Возвращает указатель на параметры расчета НДС
	* @return nullptr,если нет НДС решателя
	*/
	FStressDeformationParams* GetStressDeformationParams() {return HasSolver(MPHStressDeformation) ? &_stressDeformationParams : nullptr;}
	void SetStressDeformationParams(FStressDeformationParams* val){_stressDeformationParams = *val;}	
	
	set<int> &GetSolverTypes(){return _solverTypes ;}
	void SetSolverTypes(set<int>* val) {_solverTypes = *val;}
	
	/** Загрузка из текстового потока
	* @param ifs - входной текстовый поток
	*/
	void Load(ifstream& ifs);

	/** Сохранение в текстовый поток
	* @param ofs - выходной текстовый поток
	*/
	void Save(ofstream& ofs);
};


#endif // FMULTIPHYSICS_H