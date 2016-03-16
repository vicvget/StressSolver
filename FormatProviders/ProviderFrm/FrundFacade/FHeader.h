#ifndef FHEADER_H

#define FHEADER_H


#include "BasicTypes.h"
#include "FVariablesGroupsSet.h"
#include "FSolveParamsSet.h"
#include "FMeasureParamsSet.h"
#include "FMacro.h"
#include "FElement.h"
#include "../../../fcore/fcore.h" // for filenames

#include "FMphHeader.h"


//TODO: нормализовать комментарии

/**
* Класс для считывания заголовка модели ФРУНД
* В заголовке лежат все наборы параметров и информация о макросах
*
* @author Getmanskiy Victor
*/
class FHeader : public FElement
{
private:
	// имя задачи
	string _taskName;
	// код задачи
	string _taskCode;
	// макросы
	vector<FMacro> _macros;

	FMphHeader* _mphHeader;

	// наборы переменных
	FVariablesGroupsSet variablesGroupsSet;
	// наборы параметров расчета
	FSolveParamsSet solveParamsSet;
	// наборы параметров анализа
	FMeasureParamsSet measureParamsSet;
#pragma region File operation
	/** чтение переменных
	* @param stream - входной поток
	*/
	void ReadVariables(ifstream &stream);
	/** чтение макросов
	* @param stream - входной поток
	*/
	void ReadMacros(ifstream &stream);
	/** чтение параметров расчета
	* @param stream - входной поток
	*/
	void ReadSolverParams(ifstream &stream);
	/** чтение параметров анализа
	* @param stream - входной поток
	*/
	void ReadMeasureParams(ifstream &stream);

	/** запись переменных
	* @param stream - выходной поток
	*/
	void WriteVariables(ofstream &stream) const;
	/** запись макросов
	* @param stream - выходной поток
	*/
	void WriteMacros(ofstream &stream) const;
	/** запись параметров расчета
	* @param stream - выходной поток
	*/
	void WriteSolverParams(ofstream &stream) const;
	/** запись парамтров анализа
	* @param stream - выходной поток
	*/
	void WriteMeasureParams(ofstream &stream) const;
#pragma endregion

public:
	
	FHeader();

	//FHeader(const FHeader& src) {};

	//FHeader(const char* path);
	//FHeader(ifstream& stream);

	FHeader
		(
			const char* path
		)
		:
			_mphHeader(nullptr)
	{
		FElement::Load(path);
	}

	FHeader
		(
			ifstream& ifs
		)
		:
			_mphHeader(nullptr)
	{
		Load(ifs);
	}

	~FHeader();

	bool IsTireExistInMacros() const;

	void AddCurrentIcoParamsSet(const string& name, const IcoParams& icoParams)
	{
		size_t icoParamsId = solveParamsSet.GetIcoParamsCount();
		size_t currentId = solveParamsSet.GetParamSetIdsSetsCount();
		solveParamsSet.AddIcoParamsSet(icoParams);
		
		ParamSetIds newParamSetIds = *solveParamsSet.GetParamSetIdsSet();
		newParamSetIds.name = name;
		newParamSetIds.IcoParamsId = icoParamsId;
		solveParamsSet.AddParamSetIdsSet(newParamSetIds);
		solveParamsSet.CurrentId(currentId);
	}

	string GetDetiledRoadProfile() const
	{
		return _taskCode.substr(0, 5) + "949.gmr";
	}

	string GetRoadProfile() const
	{
		return _taskCode.substr(0, 6) + "49.gmr";
	}

	/** Получить имя задачи
	* @return имя задачи
	*/
	string TaskName() const { return _taskName; }
	
	/** Установить имя задачи
	* @param val - имя задачи
	*/
	void TaskName(string val) { _taskName = val; }

	/** Получить Код задачи
	* @return код задачи
	*/
	string TaskCode() const { return _taskCode; }
	
	/** Установить код задачи
	* @param val - код задачи
	*/
	void TaskCode(string val) { _taskCode = val; }


	FVariablesGroupsSet VariablesGroupsSet() const { return variablesGroupsSet; }
	void VariablesGroupsSet(FVariablesGroupsSet val) { variablesGroupsSet = val; }


	FMeasureParamsSet MeasureParamsSet() const { return measureParamsSet; }
	void MeasureParamsSet(FMeasureParamsSet val) { measureParamsSet = val; }
	
	/** Получить набор групп параметров расчета
	* @return набор групп параметров расчета
	*/
	FSolveParamsSet& GetSolveParamsSet()
	{
		return solveParamsSet;
	}

	const FSolveParamsSet& GetSolveParamsSet() const
	{
		return solveParamsSet;
	}
	
	/** Установить набор групп параметров расчета
	* @param newSolveParamsSet - набор групп параметров расчета
	*/
	void SetSolveParamsSet(FSolveParamsSet& newSolveParamsSet) {solveParamsSet = newSolveParamsSet;}

	/** Возвращает набор управляющих параметров расчета (пробрасывание UPRF во FRUND Solver)
	* @param id - номер набора управляющих параметров
	* @return указатель на структуру управляющих параметров
	*/
	ControlParams* GetSolveControlParams(int id = -1) {return GetSolveParamsSet().GetControlParamsSet(id);}

	/** Возвращает набор специальных параметров расчета
	* @param id - номер набора специальных параметров
	* @return указатель на структуру специальных параметров
	*/
	SpecialParams* GetSolveSpecialParams(int id = -1)
	{
		return GetSolveParamsSet().GetSpecialParamsSet(id);
	}

	const SpecialParams* GetSolveSpecialParams(int id = -1) const
	{
		return GetSolveParamsSet().GetSpecialParamsSet(id);
	}

	/** Возвращает начальные условия
	* @param id - номер варианта начальных условий
	* @return указатель на структуру с начальными условиями
	*/
	IcoParams* GetIcoParams(int id = -1) {return GetSolveParamsSet().GetIcoParamsSet(id);}

	/** Получить группу переменных по индексу группы
	* @param id - индексу группы перменных
	* @return группа переменных
	*/
	FVariablesGroup& GetVariablesGroup(size_t index) {return variablesGroupsSet[index];}

	/** Получить группу переменных по индексу группы
	* @param id - индексу группы перменных
	* @return группа переменных
	*/
	const FVariablesGroup& GetVariablesGroup(size_t index) const { return variablesGroupsSet[index]; }

	/** Получить используемую группу переменных
	* @return группа переменных
	*/
	const FVariablesGroup& GetCurrentVariablesGroup() const {return variablesGroupsSet.CurrentVariablesGroup();}

	/** Получить используемую группу переменных
	* @return группа переменных
	*/
	FVariablesGroup& GetCurrentVariablesGroup() {return variablesGroupsSet.CurrentVariablesGroup();}

	/** Получить количество переменных в группе переменных по индексу группы
	* @param id - индексу группы перменных
	* @return количество переменных в группе
	*/
	size_t GetVariablesCount(size_t index) {return variablesGroupsSet[index].Size();}

	/** Получить набор групп переменных модели
	* @return набор групп переменных
	*/
	FVariablesGroupsSet& GetVariablesGroupsSet() {return variablesGroupsSet;}
	
	/** Получить набор групп переменных модели
	* @return набор групп переменных
	*/
	const FVariablesGroupsSet& GetVariablesGroupsSet() const { return variablesGroupsSet; }

	/** Получить количество групп переменных
	* @return количество групп переменных
	*/
	size_t GetVariablesGroupsCount() const {return variablesGroupsSet.Size();}

	/** Получить количество наборов замеров
	* @return количество наборов замеров
	*/
	size_t GetMeasureParamsSetsCount() const {return measureParamsSet.Size();}

	/** Получить количество макросов
	* @return количество макросов
	*/
	size_t GetMacroCount() const {return _macros.size();}

	/**Получить макрос по индексу
	* @retrun Макрос
	*/
	FMacro& GetFMacro(size_t index) {return _macros[index];}

	/** Получить междисциплинарные параметры
	* @return междисциплинарные параметры или nullptr
	*/
	FMphHeader* GetMph() {return _mphHeader;}

	/** Добавить к модели междисциплинарные параметры
	*/
	void AddMph() {_mphHeader = new FMphHeader();}
	

// См. комментарии в базовом классе
#pragma region overriden

	virtual
	bool Load
		(
			ifstream& stream
		)	override;

	void Load(const char* path);
	void Save(ofstream& ofs) const;
	void SaveByIndex() const;

	virtual
	void Output
		(
			ofstream& stream,
			int stage
		)	const override;

#pragma endregion

	// Вывод в текстовый файловый поток переменных текущего набора
	void OutputVariables(ofstream &stream, int isRType) const;
	// Вывод в текстовый файловый поток макросов
	void OutputMacros(ofstream &stream, int isRType) const;
	// Вывод управляющих параметров для расчета (UPRF)
	void OutputSolveControlParams(ofstream &stream) const;
	// Вывод управляющих параметров для анализа (UPRF)
	void OutputAnazControlParams(ofstream &stream) const;

	/** Формирует управляющие файлы для шага расчет
	* @param isRoadExist - признак наличия дороги в модели
	*/
	void OutputSolveParams(bool isRoadExist) const;

	/** Выбрать набор управляющих параметров для расчета
	* @param paramsSetId - номер набора параметров
	*/
	void SetSolveParamsSet(int paramsSetId);

	/** Выбрать набор управляющих параметров для построения графиков
	* @param paramsSetId - номер набора параметров
	*/
	void SetPostprocessorControlParamsSet(int paramsSetId);

	/** Выбрать набор параметров модели
	* @param paramsSetId - номер набора параметров
	*/
	void SetVariablesParamsSet(int paramsSetId);
};


#endif // FHEADER_H