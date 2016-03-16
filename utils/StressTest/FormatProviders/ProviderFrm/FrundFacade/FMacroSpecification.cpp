#include "FMacroSpecification.h"

#include "../../../fcore/fcore.h"
#include "../../../fcore/exceptions/fcExceptions.h"
#include "../../../fcore/wrappers/FileRoutines.h"


// class MacroSolveParam

/**
* Конструктор
*/
MacroSolveParam::MacroSolveParam()
{

}

/**
* Получить наименование параметра
* @return наименование параметра
*/
const string& MacroSolveParam::GetName() const
{
	return _name;
}

/**
* Получить значение по умолчанию
* @return значение по умолчанию
*/
const string& MacroSolveParam::GetDefaultValue() const
{
	return _defaultValue;
}

/**
* Получить значение по умолчанию в качестве вещественного числа
* @return значение по умолчанию в качестве вещественного числа
*/
double MacroSolveParam::GetDefaultValueAsDouble() const
{
	return atof(_defaultValue.c_str());
}

/**
* Чтение параметра макроса шага "расчет" из xml-документа
* @param macroSolveParamElement - элемент xml-документа, содержащий параметр макроса шага "расчет"
* @return признак успешного (true) или неуспешного (false) чтения
*/
bool MacroSolveParam::FromXmlDocument
	(
		const TiXmlElement* macroSolveParamElement
	)
{
	if (string(macroSolveParamElement->Value()) != "SolveParam")
	{
		return false;
	}

	const char* name;
	
	name = macroSolveParamElement->Attribute("name");
	if (name == nullptr)
	{
		return false;
	}

	const TiXmlElement* defaultValueElement;
	
	defaultValueElement = macroSolveParamElement->FirstChildElement("default");
	if (defaultValueElement == nullptr)
	{
		return false;
	}

	const char* defaultValue;

	defaultValue = defaultValueElement->GetText();
	if (defaultValue == nullptr)
	{
		return false;
	}
	_name = name;
	_defaultValue = defaultValue;

	return true;
}


// vector MacroSolveParams

/**
* Чтение списка параметров макроса шага "расчет" из xml-документа
* @param macroSolveParamsElement - элемент xml-документа, содержащий список параметров макроса
* @return признак успешного (true) или неуспешного (false) чтения
*/
bool FromXmlDocument
	(
		const TiXmlElement* macroSolveParamsElement,
		MacroSolveParams& macroSolveParams
	)
{
	if (string(macroSolveParamsElement->Value()) != "SolveParams")
	{
		return false;
	}

	MacroSolveParams macroSolveParamsTmp;
	const TiXmlElement* macroSolveParamElement;

	macroSolveParamElement = macroSolveParamsElement->FirstChildElement();
	while (macroSolveParamElement != nullptr)
	{
		MacroSolveParam macroSolveParam;

		if (!macroSolveParam.FromXmlDocument(macroSolveParamElement))
		{
			return false;
		}
		macroSolveParamsTmp.push_back(macroSolveParam);
		macroSolveParamElement = macroSolveParamElement->NextSiblingElement();
	}
	macroSolveParams = macroSolveParamsTmp;

	return true;
}


// class MacroModelParam

/**
* Конструктор
*/
MacroModelParam::MacroModelParam()
	:
		_type(Type::Unknown)
{

}

/**
* Получить наименование параметра
* @return наименование параметра
*/
const string& MacroModelParam::GetName() const
{
	return _name;
}

/**
* Получить тип параметра
* @return тип параметра
*/
MacroModelParam::Type MacroModelParam::GetType() const
{
	return _type;
}

/**
* Получить признак, является ли данный параметр номером тела
* @return признак, является ли (true) данный параметр номером тела или нет (false)
*/
bool MacroModelParam::IsBodyNumber() const
{
	return _type == Type::BodyNumber;
}

/**
* Получить признак, является ли данный параметр номером узла тела
* @return признак, является ли (true) данный параметр номером узла тела или нет (false)
*/
bool MacroModelParam::IsNodeNumber() const
{
	return _type == Type::NodeNumber;
}

/**
* Получить наименование переменной
* @return наименование переменной
*/
const string& MacroModelParam::GetVariableName() const
{
	return _variableName;
}

/**
* Получить значение по умолчанию
* @return значение по умолчанию
*/
const string& MacroModelParam::GetDefaultValue() const
{
	return _defaultValue;
}

/**
* Получить значение по умолчанию в качестве целого числа
* @return значение по умолчанию в качестве целого числа
*/
int MacroModelParam::GetDefaultValueAsInteger() const
{
	return atoi(_defaultValue.c_str());
}

/**
* Чтение параметра макроса этапа подготовки модели из xml-документа
* @param macroModelParamElement - элемент xml-документа, содержащий параметр макроса
* @return признак успешного (true) или неуспешного (false) чтения
*/
bool MacroModelParam::FromXmlDocument
	(
		const TiXmlElement* macroModelParamElement
	)
{
	if (string(macroModelParamElement->Value()) != "Param")
	{
		return false;
	}

	const char* name;
	
	name = macroModelParamElement->Attribute("name");
	if (name == nullptr)
	{
		return false;
	}

	const char* typeAsString;
	
	typeAsString = macroModelParamElement->Attribute("type");
	if (typeAsString == nullptr)
	{
		return false;
	}

	int type;

	type = atoi(typeAsString);

	const char* variableName;
	
	variableName = macroModelParamElement->Attribute("var");
	if (variableName == nullptr)
	{
		return false;
	}
	
	const TiXmlElement* defaultValueElement;
	
	defaultValueElement = macroModelParamElement->FirstChildElement("default");
	if (defaultValueElement == nullptr)
	{
		return false;
	}

	const char* defaultValue;

	defaultValue = defaultValueElement->GetText();
	if (defaultValue == nullptr)
	{
		return false;
	}
	_name = name;
	_type = static_cast<Type>(type);
	_variableName = variableName;
	_defaultValue = defaultValue;

	return true;
}


// vector MacroModelParams

/**
* Чтение списка параметров макроса этапа подготовки модели из xml-документа
* @param macroModelParamsElement - элемент xml-документа, содержащий список параметров макроса
* @return признак успешного (true) или неуспешного (false) чтения
*/
bool FromXmlDocument
	(
		const TiXmlElement* macroModelParamsElement,
		MacroModelParams& macroModelParams
	)
{
	if (string(macroModelParamsElement->Value()) != "ModelParams")
	{
		return false;
	}

	MacroModelParams macroModelParamsTmp;
	const TiXmlElement* macroModelParamElement;

	macroModelParamElement = macroModelParamsElement->FirstChildElement();
	while (macroModelParamElement != nullptr)
	{
		MacroModelParam macroModelParam;

		if (!macroModelParam.FromXmlDocument(macroModelParamElement))
		{
			return false;
		}
		macroModelParamsTmp.push_back(macroModelParam);
		macroModelParamElement = macroModelParamElement->NextSiblingElement();
	}
	macroModelParams = macroModelParamsTmp;

	return true;
}


// class FMacroSpecification

/**
* Конструктор
*/
FMacroSpecification::FMacroSpecification()
	:
		_type()
{

}

/**
* Получить тип макроса
* @return тип макроса
*/
int FMacroSpecification::GetType() const
{
	return _type;
}

/**
* Получить имя файла с макросом
* @return имя файла с макросом
*/
const string& FMacroSpecification::GetFileName() const
{
	return _fileName;
}

/**
* Получить описание макроса
* @return описание макроса
*/
const string& FMacroSpecification::GetDescription() const
{
	return _description;
}

/**
* Получить короткое описание макроса
* @return короткое описание макроса
*/
const string& FMacroSpecification::GetShortDescription() const
{
	return _shortDescription;
}

/**
* Получить ??? описание макроса
* @return ??? описание макроса
*/
const string& FMacroSpecification::GetIndexDescription() const
{
	return _indexDescription;
}

/**
* Получить количество параметров макроса шага "расчет"
* @return количество параметров макроса шага "расчет"
*/
int FMacroSpecification::GetSolveParamsCount() const
{
	return _solveParams.size();
}

/**
* Получить параметр макроса шага "расчет" по его индексу
* @param paramIndex - индекс параметра макроса
* @return параметр макроса шага "расчет" с данным индексом
*/
const MacroSolveParam* FMacroSpecification::GetSolveParam
	(
		int paramIndex
	)	const
{
	if ((paramIndex < 0) || (paramIndex >= _solveParams.size()))
	{
		return nullptr;
	}

	return &_solveParams.at(paramIndex);
}

/**
* Получить количество параметров макроса этапа подготовки модели
* @return количество параметров макроса этапа подготовки модели
*/
int FMacroSpecification::GetModelParamsCount() const
{
	return _modelParams.size();
}

/**
* Получить параметр макроса этапа подготовки модели по его индексу
* @param paramIndex - индекс параметра макроса
* @return параметр макроса этапа подготовки модели с данным индексом
*/
const MacroModelParam* FMacroSpecification::GetModelParam
	(
		int paramIndex
	)	const
{
	if ((paramIndex < 0) || (paramIndex >= _modelParams.size()))
	{
		return nullptr;
	}

	return &_modelParams.at(paramIndex);
}

/**
* Чтение описания макроса из xml-документа
* @param macroElement - элемент xml-документа, содержащий описание макроса
* @return признак успешного (true) или неуспешного (false) чтения
*/
bool FMacroSpecification::FromXmlDocument
	(
		const TiXmlElement* macroElement
	)
{
	if (string(macroElement->Value()) != "MacrosResource")
	{
		return false;
	}

	const char* typeAsString;
	
	typeAsString = macroElement->Attribute("type");
	if (typeAsString == nullptr)
	{
		return false;
	}

	int type;

	type = atoi(typeAsString);

	const char* fileName;
	
	fileName = macroElement->Attribute("file");
	if (fileName == nullptr)
	{
		return false;
	}

	const TiXmlElement* descriptionElement;
	
	descriptionElement = macroElement->FirstChildElement("Description");
	if (descriptionElement == nullptr)
	{
		return false;
	}

	const char* description;

	description = descriptionElement->GetText();
	if (description == nullptr)
	{
		return false;
	}

	const TiXmlElement* shortDescriptionElement;
	
	shortDescriptionElement = macroElement->FirstChildElement("ShortDescription");
	if (shortDescriptionElement == nullptr)
	{
		return false;
	}

	const char* shortDescription;

	shortDescription = shortDescriptionElement->GetText();
	if (shortDescription == nullptr)
	{
		return false;
	}

	const TiXmlElement* indexDescriptionElement;
	
	indexDescriptionElement = macroElement->FirstChildElement("IndexDescription");
	if (indexDescriptionElement == nullptr)
	{
		return false;
	}

	const char* indexDescription;

	indexDescription = indexDescriptionElement->GetText();
	if (indexDescription == nullptr)
	{
		return false;
	}

	const TiXmlElement* macroSolveParamsElement;

	macroSolveParamsElement = macroElement->FirstChildElement("SolveParams");
	if (macroSolveParamsElement == nullptr)
	{
		return false;
	}

	MacroSolveParams solveParams;

	if (!::FromXmlDocument(macroSolveParamsElement, solveParams))
	{
		return false;
	}

	const TiXmlElement* macroModelParamsElement;

	macroModelParamsElement = macroElement->FirstChildElement("ModelParams");
	if (macroModelParamsElement == nullptr)
	{
		return false;
	}

	MacroModelParams modelParams;

	if (!::FromXmlDocument(macroModelParamsElement, modelParams))
	{
		return false;
	}
	_type = type;
	_fileName = fileName;
	_description = description;
	_shortDescription = shortDescription;
	_indexDescription = indexDescription;
	_solveParams = solveParams;
	_modelParams = modelParams;

	return true;
}


// map MacroSpecificationMap

/**
* Чтение карты описаний макросов из xml-документа
* @param macroSpecificationMapElement - элемент xml-документа, содержащий описания макросов
* @return признак успешного (true) или неуспешного (false) чтения
*/
bool FromXmlDocument
	(
		const TiXmlElement* macroSpecificationMapElement,
		MacroSpecificationMap& macroSpecificationMap
	)
{
	if (string(macroSpecificationMapElement->Value()) != "MacrosResources")
	{
		return false;
	}

	MacroSpecificationMap macroSpecificationMapTmp;
	const TiXmlElement* macroElement;

	macroElement = macroSpecificationMapElement->FirstChildElement();
	while (macroElement != nullptr)
	{
		FMacroSpecification macroSpecification;

		if (!macroSpecification.FromXmlDocument(macroElement))
		{
			return false;
		}

		const string& fileName = macroSpecification.GetFileName();

		if (macroSpecificationMapTmp.find(fileName) != macroSpecificationMapTmp.end())
		{
			return false;
		}
		macroSpecificationMapTmp.insert
			(
				MacroSpecificationPair
					(
						fileName,
						macroSpecification
					)
			);
		macroElement = macroElement->NextSiblingElement();
	}
	macroSpecificationMap = macroSpecificationMapTmp;

	return true;
}

/**
* Чтение карты описаний макросов из xml-документа
* @param xmlFileName - наименование xml-документа, содержащиго описания макросов
* @return признак успешного (true) или неуспешного (false) чтения
*/
bool FromXmlDocument
	(
		const string& xmlFileName,
		MacroSpecificationMap& macroSpecificationMap
	)
{
	TiXmlDocument xmlFile(xmlFileName.c_str());

	if (!xmlFile.LoadFile())
	{
		exceptions::ThrowFileInvalidFormat(xmlFileName);
	}

	const TiXmlElement* resourcesElement;

	resourcesElement = xmlFile.RootElement();
	if (resourcesElement == nullptr)
	{
		return false;
	}
	if (string(resourcesElement->Value()) != "FrundResources")
	{
		return false;
	}

	const TiXmlElement* macrosElement;

	macrosElement = resourcesElement->FirstChildElement("MacrosResources");
	if (macrosElement == nullptr)
	{
		return false;
	}
	
	return FromXmlDocument(macrosElement, macroSpecificationMap);
}

/**
* Чтение карты описаний макросов из файла ресурсов ФРУНД
* @return признак успешного (true) или неуспешного (false) чтения
*/
bool FromResourcesFile
	(
		MacroSpecificationMap& macroSpecificationMap
	)
{
	macroSpecificationMap.clear();
	try
	{
		string resourcesFileName = fs::CombinePathEnv(DIR_DATA, FILE_RESOURCES);

		if (!fs::FileExist(resourcesFileName))
		{
			exceptions::ThrowFileNotFound(resourcesFileName);
		}
		if
			(
				!FromXmlDocument
					(
						resourcesFileName,
						macroSpecificationMap
					)
			)
		{
			exceptions::ThrowFileNotOpened(resourcesFileName);
		}
	}
	catch (const exceptions::CoreException& exception)
	{
		exceptions::Output(exception);

		return false;
	}

	return true;
}