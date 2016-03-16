#include "FVariablesGroupsSet.h"

#include "../../../fcore/fcore.h"

#include <fstream>


using std::endl;


FVariablesGroup::FVariablesGroup()
{
	//calculator = new Calc::Calculator();
}

FVariablesGroup::FVariablesGroup(ifstream& stream, int count)
{
	Parameter param;
	//calculator = new Calc::Calculator();
	for(int i = 0; i < count; i++)
	{
		fs::ReadLineString(stream, param.name);
		fs::ReadLineString(stream, param.expression);
		fs::ReadLineString(stream, param.description);
		parameters.push_back(param);
		//calculator->AddParameter(param.name, param.expression);
	}
}

/** Копирование группы переменных
* @param variableGroup - копируемая группа
*/
FVariablesGroup::FVariablesGroup
	(
		const FVariablesGroup& variablesGroup
	)
	:
		name(variablesGroup.name),
		parameters(variablesGroup.parameters)
{
}

FVariablesGroup::~FVariablesGroup()
{
	//delete calculator;
}

// Копировать набор параметров
void FVariablesGroup::Copy(ifstream& stream, FVariablesGroup& ps)
{
	parameters.assign(ps.parameters.begin(), ps.parameters.end());
	for(int i = 0; i < parameters.size(); i++)
	{
		fs::ReadLineString(stream, parameters[i].expression);
		//stream.getline(parameters[i].expression, STRLIMIT);
		//calculator->AddParameter(parameters[i].name, parameters[i].expression);
	}
	
}

// Вывод в текстовый файловый поток
void FVariablesGroup::Output(ofstream& stream) const
{
	for(int i = 0; i < parameters.size(); i++)
	{
		stream << parameters[i].name << '=' 
			<< parameters[i].expression << " ! " 
			<< parameters[i].description << " \n";
	}
}

void FVariablesGroup::Write(ofstream& ofs) const
{
	for(int i = 0; i < parameters.size(); i++)
	{
		ofs << parameters[i].name << endl <<
			parameters[i].expression << endl <<
			parameters[i].description << endl;
	}
}
void FVariablesGroup::WriteValues(ofstream& ofs) const
{
	for(int i = 0; i < parameters.size(); i++)
	{
		ofs << parameters[i].expression << endl;
	}
}

void FVariablesGroupsSet::Write(ofstream& ofs) const
{
	size_t size = 0;

	// TODO: CurrentId имеет тип std::size_t, поэтому не может быть меньше нуля
	if(!(CurrentId()<0))
	{
		size = CurrentVariablesGroup().Size();
		ofs << size << _tag << endl;
		CurrentVariablesGroup().Write(ofs);
		ofs << Size() << ' ' << currentId << endl;
		for(int i = 0; i < Size(); i++)
			ofs << variablesGroups[i].Name() << endl;
		if(Size() > 1) // если меньше 1 значения, то значения не выводятся вообще
			for(int i = 0; i < Size(); i++)
				variablesGroups[i].WriteValues(ofs);
	}
	else
	{
		ofs << size << _tag << endl;
	}
}

/** Считывание параметров из входного потока
* @param ifs - input file stream
*/
void FVariablesGroupsSet::Read(ifstream& stream)
{
	int paramsCount,paramsSetsCount;
	stream >> paramsCount;
	fs::ReadLineString(stream,_tag);
	FVariablesGroup vairableGroup(stream, paramsCount);
	variablesGroups.clear();
	variablesGroups.push_back(vairableGroup);

	char buf[STRLIMIT];

	stream >> paramsSetsCount >> currentId;
	stream.getline(buf, STRLIMIT);
	stream.getline(buf, STRLIMIT);
	variablesGroups[0].SetName(buf);
	if (paramsSetsCount > 1)
	{
		int i = 0;

		while (i < paramsSetsCount - 1)
		{
			FVariablesGroup vairablesGroup1;
			stream.getline(buf, STRLIMIT);
			vairablesGroup1.SetName(buf);
			variablesGroups.push_back(vairablesGroup1);
			i++;
		}
		i = 0;
		while (i < paramsSetsCount)
		{
			variablesGroups[i++].Copy(stream, vairableGroup);
		}
	}
}

/** Добавляет во все списки переменных новую запись 
*(предусмотреть контроль на уникальность имен в оболочке)
* @param newVariable - переменная
*/
void FVariablesGroupsSet::AddNewVariable(Parameter newVariable)
{
	vector<FVariablesGroup>::iterator it = variablesGroups.begin();

	while (it != variablesGroups.end())
	{
		it->AddNewVariable(newVariable);
		++it;
	}
}

/** Удаляет переменную по индексу из всех списков
* @param index - индекс переменной
*/
void FVariablesGroupsSet::DeleteVariable(int id)
{
	vector<FVariablesGroup>::iterator it = variablesGroups.begin();

	while (it != variablesGroups.end())
	{
		it->DeleteVariable(id);
		++it;
	}
}

/** Добавляет в список переменных новую запись (предусмотреть контроль на уникальность имен в оболочке)
* не использовать эту функцию из внешних модулей!
* @param newVariable - переменная
*/
void FVariablesGroup::AddNewVariable(Parameter newVariable)
{
	parameters.push_back(newVariable);
}

/** Удаляет переменную по индексу
* не использовать эту функцию из внешних модулей!
* @param index - индекс переменной
*/
void FVariablesGroup::DeleteVariable(int index)
{
	if(parameters.size() > index)
		parameters.erase(parameters.begin()+index);
}


/** Модификация переменной в группе
* expression модифицируется только для группы, name и description - для всех групп
* @param newVariable - измененная переменная
* @param groupId - индекс группы
* @param varId - индекс переменной
*/
void FVariablesGroupsSet::ModifyVariable(Parameter newVariable, int groupId, int varId)
{
	variablesGroups[groupId].ModifyVariableExpression(newVariable.expression,varId);
	if (strcmp(variablesGroups[groupId][varId].description.c_str(), newVariable.description.c_str()))
	{
		vector<FVariablesGroup>::iterator it = variablesGroups.begin();

		while (it != variablesGroups.end())
		{
			it->ModifyVariableDescription(newVariable.description, varId);
			++it;
		}
	}
	if (strcmp(variablesGroups[groupId][varId].name.c_str(), newVariable.name.c_str()))
	{
		vector<FVariablesGroup>::iterator it = variablesGroups.begin();

		while (it != variablesGroups.end())
		{
			it->ModifyVariableDescription(newVariable.name, varId);
			++it;
		}
	}
}

/** Модификация имени переменной
* @param newVariableName - имя переменной
* @param varId - индекс переменной
*/
void FVariablesGroup::ModifyVariableName(string newVariableName, int varId)
{
	//strcpy(parameters[varId].name.c_str(), newVariableName);
	parameters[varId].name = newVariableName;
}
	
/** Модификация описания переменной
* @param newVariableDescription - описание переменной
* @param varId - индекс переменной
*/
void FVariablesGroup::ModifyVariableDescription(string newVariableDescription, int varId)
{
	//strcpy(parameters[varId].description, newVariableDescription.c_str());
	parameters[varId].description = newVariableDescription;

}
	
/** Модификация значения переменной
* @param newVariableExpression - значение переменной
* @param varId - индекс переменной
*/
void FVariablesGroup::ModifyVariableExpression(string newVariableExpression, int varId)
{
	//strcpy(parameters[varId].expression, newVariableExpression.c_str());
	parameters[varId].expression = newVariableExpression;
}

/** Модификация имени группы переменных
* @param newName - новое имя группы
* @param groupId - индекс группы
*/
void FVariablesGroupsSet::ModifyVariablesGroupName(string newName, int groupId)
{
	variablesGroups[groupId].SetName(newName);
}

/** Добавляет ноаую группу переменных
* @param newName - имя группы
*/
void FVariablesGroupsSet::AddNewVariablesGroup(string newName)
{
	FVariablesGroup newVariablesGroup;
	if(CurrentId() > 0)
	{
		newVariablesGroup = CurrentVariablesGroup();
	}
	newVariablesGroup.SetName(newName);
	variablesGroups.push_back(newVariablesGroup);
}

/** Удаляет группу переменных по индексу
* @param groupId - индекс группы
*/
void FVariablesGroupsSet::DeleteVariablesGroup(int groupId)
{
	variablesGroups.erase(variablesGroups.begin()+groupId);
}
