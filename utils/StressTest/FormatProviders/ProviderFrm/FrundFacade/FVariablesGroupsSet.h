#ifndef FVARIABLES_GROUPS_SET_H

#define FVARIABLES_GROUPS_SET_H


#include "BasicTypes.h"
#include "FParamsSet.h"
#include "../../../fcore/Calculator/Calculator.h"


/**
* Класс для хранения набора параметров модели ФРУНД
*
* @author Getmanskiy Victor
*/
class FVariablesGroup
{

private:
	// имя группы переменных
	string name;
	// параметры
	vector<Parameter> parameters;
public:
	FVariablesGroup();
	~FVariablesGroup();

	const Parameter& operator [] (size_t id) const { return parameters[id]; }
	const vector<Parameter>& Parameters() const { return parameters; }
	void Parameters(const vector<Parameter>& val) { parameters = val; }
	//void CopyParametersTo(vector<Parameter>& parametersVector)
	//	{parametersVector.assign(parameters.begin(), parameters.end());}

	/** Создание группы переменных из фалового потока
	* @param stream - файловый поток с указателем на начало секции параметров
	* @param count - количество параметров
	*/
	FVariablesGroup(ifstream& stream, int count);
	
	/** Копирование группы переменных
	* @param variableGroup - копируемая группа
	*/
	FVariablesGroup(const FVariablesGroup& variableGroup);

	// Указатель на структуры с переменными
	Parameter* GetParamsPointer() {return &parameters[0];}

	// Имя набора переменных
	string Name() const {return name;}

	/** Установка имени набора параметра
	* @param name - имя группы переменных
	*/
	void SetName(string name) {this->name = name;}
	
	size_t Size() const {return parameters.size();}

	/** Копирование группы переменных из файла
	* @param stream - файловый поток с указателем на начало секции параметров
	* @ps - набор параметров, в который копируем
	*/
	void Copy(ifstream& stream, FVariablesGroup& ps);

	/** Вывод в текстовый файловый поток для mdl,ras
	* @param stream - файловый поток
	*/
	void Output(ofstream &stream) const;

	/** Запись группы переменных в выходной поток
	* @param ofs - выходной поток
	*/
	void Write(ofstream& ofs) const;

	/** Запись выражений для группы переменных в выходной поток
	* @param ofs - выходной поток
	*/
	void WriteValues(ofstream& ofs) const;


	/** Добавляет в список переменных новую запись (предусмотреть контроль на уникальность имен в оболочке)
	* не использовать эту функцию из внешних модулей!
	* @param newVariable - переменная
	*/
	void AddNewVariable(Parameter newVariable);

	/** Удаляет переменную по индексу
	* не использовать эту функцию из внешних модулей!
	* @param index - индекс переменной
	*/
	void DeleteVariable(int index);

	/** Модификация имени переменной
	* @param newVariableName - имя переменной
	* @param varId - индекс переменной
	*/
	void ModifyVariableName(string newVariableName, int varId);
	
	/** Модификация описания переменной
	* @param newVariableDescription - описание переменной
	* @param varId - индекс переменной
	*/
	void ModifyVariableDescription(string newVariableDescription, int varId);
	
	/** Модификация значения переменной
	* @param newVariableExpression - значение переменной
	* @param varId - индекс переменной
	*/
	void ModifyVariableExpression(string newVariableExpression, int varId);
};

/**
* Класс для хранения всех наборов параметров модели ФРУНД
*
* @author Getmanskiy Victor
*/
class FVariablesGroupsSet: public FParamsSet
{
	vector<FVariablesGroup> variablesGroups;
	// & ... в третьей строке заголовка
	string _tag;
public:

	FVariablesGroupsSet()
		:
			_tag(" & 1 5 10.000000 10.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000")
	{
		SetCurrentId(0);
		AddNewVariablesGroup("default");
	}

	size_t Size() const {return variablesGroups.size();}
	FVariablesGroup& CurrentVariablesGroup() {return variablesGroups.at(currentId - 1);}
	const FVariablesGroup& CurrentVariablesGroup() const {return variablesGroups.at(currentId - 1);}
	vector<FVariablesGroup> VariablesGroups() const { return variablesGroups; }
	void VariablesGroups(vector<FVariablesGroup> val) { variablesGroups = val; }
	
	size_t CurrentId() const {return currentId - 1;}
	void SetCurrentId(size_t newId) {currentId = newId + 1;}
	
	FVariablesGroup& operator [] (size_t id) { return variablesGroups[id]; }
	const FVariablesGroup& operator [] (size_t id) const { return variablesGroups[id]; }
	//public bool Select(int id) {currentSetId = id;}

	/** Считывание наборов параметров из входного потока
	* @param ifs - входной поток
	*/
	void Read(ifstream& ifs);
	
	/** Запись наборов параметров в выходной поток
	* @param ofs - выходной поток
	*/
	void Write(ofstream& ofs) const;


	/** Добавляет во все списки переменных новую запись 
	*(предусмотреть контроль на уникальность имен в оболочке)
	* @param newVariable - переменная
	*/
	void AddNewVariable(Parameter newVariable);

	/** Удаляет переменную по индексу из всех списков
	* @param index - индекс переменной
	*/
	void DeleteVariable(int index);

	/** Модификация переменной в группе
	* expression модифицируется только для группы, name и description - для всех групп
	* @param newVariable - измененная переменная
	* @param groupId - индекс группы
	* @param varId - индекс переменной
	*/
	void ModifyVariable(Parameter newVariable, int groupId, int varId);

	/** Модификация имени группы переменных
	* @param newName - новое имя группы
	* @param groupId - индекс группы
	*/
	void ModifyVariablesGroupName(string newName, int groupId);

	/** Добавляет ноаую группу переменных
	* @param newName - имя группы
	*/
	void AddNewVariablesGroup(string newName);

	/** Удаляет группу переменных по индексу
	* @param groupId - индекс группы
	*/
	void DeleteVariablesGroup(int groupId);

};


#endif // FVARIABLES_GROUPS_SET_H