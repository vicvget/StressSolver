#ifndef F_MACRO_SPECIFICATION_H

#define F_MACRO_SPECIFICATION_H


#include "../../../ExternalModules/TinyXml/tinyxml.h"

#include <string>
#include <vector>
#include <map>


using std::string;
using std::vector;
using std::pair;
using std::map;


/**
* Параметр макроса шага "расчет"
*/
class MacroSolveParam
{
public:

	// Конструкторы и деструктор

	/**
	* Конструктор
	*/
	MacroSolveParam();


	// Селекторы

	/**
	* Получить наименование параметра
	* @return наименование параметра
	*/
	const string& GetName() const;

	/**
	* Получить значение по умолчанию
	* @return значение по умолчанию
	*/
	const string& GetDefaultValue() const;

	/**
	* Получить значение по умолчанию в качестве вещественного числа
	* @return значение по умолчанию в качестве вещественного числа
	*/
	double GetDefaultValueAsDouble() const;


	// Чтение из файла

	/**
	* Чтение параметра макроса шага "расчет" из xml-документа
	* @param macroSolveParamElement - элемент xml-документа, содержащий параметр макроса шага "расчет"
	* @return признак успешного (true) или неуспешного (false) чтения
	*/
	bool FromXmlDocument
		(
			const TiXmlElement* macroSolveParamElement
		);


protected:

	// наименование параметра
	string _name;

	// значение по умолчанию (возможно, следует сделать double)
	string _defaultValue;

};

// список параметров макроса шага "расчет"
typedef vector<MacroSolveParam> MacroSolveParams;


// Чтение из файла

/**
* Чтение списка параметров макроса шага "расчет" из xml-документа
* @param macroSolveParamsElement - элемент xml-документа, содержащий список параметров макроса
* @return признак успешного (true) или неуспешного (false) чтения
*/
bool FromXmlDocument
	(
		const TiXmlElement* macroSolveParamsElement,
		MacroSolveParams& macroSolveParams
	);


/**
* Параметр макроса этапа подготовки модели
*/
class MacroModelParam
{
public:

	// Перечисления

	// типы параметров макроса этапа подготовки модели
	enum class Type
	{
		Unknown, // неизвестный тип
		BodyNumber = 30, // номер тела
		NodeNumber = 32 // номер узла тела
	};

	// Конструкторы и деструктор

	/**
	* Конструктор
	*/
	MacroModelParam();


	// Селекторы

	/**
	* Получить наименование параметра
	* @return наименование параметра
	*/
	const string& GetName() const;

	/**
	* Получить тип параметра
	* @return тип параметра
	*/
	Type GetType() const;

	/**
	* Получить признак, является ли данный параметр номером тела
	* @return признак, является ли (true) данный параметр номером тела или нет (false)
	*/
	bool IsBodyNumber() const;

	/**
	* Получить признак, является ли данный параметр номером узла тела
	* @return признак, является ли (true) данный параметр номером узла тела или нет (false)
	*/
	bool IsNodeNumber() const;

	/**
	* Получить наименование переменной
	* @return наименование переменной
	*/
	const string& GetVariableName() const;

	/**
	* Получить значение по умолчанию
	* @return значение по умолчанию
	*/
	const string& GetDefaultValue() const;

	/**
	* Получить значение по умолчанию в качестве целого числа
	* @return значение по умолчанию в качестве целого числа
	*/
	int GetDefaultValueAsInteger() const;


	// Чтение из файла

	/**
	* Чтение параметра макроса этапа подготовки модели из xml-документа
	* @param macroModelParamElement - элемент xml-документа, содержащий параметр макроса
	* @return признак успешного (true) или неуспешного (false) чтения
	*/
	bool FromXmlDocument
		(
			const TiXmlElement* macroModelParamElement
		);


protected:

	// наименование параметра
	string _name;

	// тип параметра
	Type _type;

	// наименование переменной
	string _variableName;

	// значение по умолчанию (возможно, следует сделать double или int)
	string _defaultValue;

};

// список параметров макроса этапа подготовки модели
typedef vector<MacroModelParam> MacroModelParams;


// Чтение из файла

/**
* Чтение списка параметров макроса этапа подготовки модели из xml-документа
* @param macroModelParamsElement - элемент xml-документа, содержащий список параметров макроса
* @return признак успешного (true) или неуспешного (false) чтения
*/
bool FromXmlDocument
	(
		const TiXmlElement* macroModelParamsElement,
		MacroModelParams& macroModelParams
	);


/**
* Класс описания макроса модели ФРУНД
*/
class FMacroSpecification
{
public:

	// Конструкторы и деструктор

	/**
	* Конструктор
	*/
	FMacroSpecification();


	// Селекторы

	/**
	* Получить тип макроса
	* @return тип макроса
	*/
	int GetType() const;

	/**
	* Получить имя файла с макросом
	* @return имя файла с макросом
	*/
	const string& GetFileName() const;

	/**
	* Получить описание макроса
	* @return описание макроса
	*/
	const string& GetDescription() const;

	/**
	* Получить короткое описание макроса
	* @return короткое описание макроса
	*/
	const string& GetShortDescription() const;

	/**
	* Получить ??? описание макроса
	* @return ??? описание макроса
	*/
	const string& GetIndexDescription() const;

	/**
	* Получить количество параметров макроса шага "расчет"
	* @return количество параметров макроса шага "расчет"
	*/
	int GetSolveParamsCount() const;

	/**
	* Получить параметр макроса шага "расчет" по его индексу
	* @param paramIndex - индекс параметра макроса
	* @return параметр макроса шага "расчет" с данным индексом
	*/
	const MacroSolveParam* GetSolveParam
		(
			int paramIndex
		)	const;

	/**
	* Получить количество параметров макроса этапа подготовки модели
	* @return количество параметров макроса этапа подготовки модели
	*/
	int GetModelParamsCount() const;

	/**
	* Получить параметр макроса этапа подготовки модели по его индексу
	* @param paramIndex - индекс параметра макроса
	* @return параметр макроса этапа подготовки модели с данным индексом
	*/
	const MacroModelParam* GetModelParam
		(
			int paramIndex
		)	const;


	// Чтение из файла

	/**
	* Чтение описания макроса из xml-документа
	* @param macroElement - элемент xml-документа, содержащий описание макроса
	* @return признак успешного (true) или неуспешного (false) чтения
	*/
	bool FromXmlDocument
		(
			const TiXmlElement* macroElement
		);


protected:

	// тип макроса
	int _type;

	// имя файла с макросом
	string _fileName;

	// описание макроса
	string _description;

	// короткое описание макроса
	string _shortDescription;

	// ??? описание макроса
	string _indexDescription;

	// список параметров макроса шага "расчет"
	MacroSolveParams _solveParams;
	
	// список параметров максроса этапа подготовки модели
	MacroModelParams _modelParams;

};

// карта описаний макросов
// ключ - имя файла с макросом
// значение - описание макроса
typedef pair<string, FMacroSpecification> MacroSpecificationPair;
typedef map<string, FMacroSpecification> MacroSpecificationMap;


// Чтение из файла

/**
* Чтение карты описаний макросов из xml-документа
* @param macroSpecificationMapElement - элемент xml-документа, содержащий описания макросов
* @return признак успешного (true) или неуспешного (false) чтения
*/
bool FromXmlDocument
	(
		const TiXmlElement* macroSpecificationMapElement,
		MacroSpecificationMap& macroSpecificationMap
	);

/**
* Чтение карты описаний макросов из xml-документа
* @param xmlFileName - наименование xml-документа, содержащиго описания макросов
* @return признак успешного (true) или неуспешного (false) чтения
*/
bool FromXmlDocument
	(
		const string& xmlFileName,
		MacroSpecificationMap& macroSpecificationMap
	);

/**
* Чтение карты описаний макросов из файла ресурсов ФРУНД
* @return признак успешного (true) или неуспешного (false) чтения
*/
bool FromResourcesFile
	(
		MacroSpecificationMap& macroSpecificationMap
	);


#endif // F_MACRO_SPECIFICATION_H