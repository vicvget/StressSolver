#include "MechanismOccurrenceData.h"

#include "../../../Fcore/Exceptions/fcExceptions.h"
#include "../../../Fcore/wrappers/FileRoutines.h"
#include "../../../ExternalModules/TinyXml/tinyxml.h"


namespace MechanismOccurrence
{

	/**
	* «агрузить данные о теле из элемента xml-документа с именем Body
	* @param bodyElement - элемент xml-документа, из которого загружаютс€ данные
	* @param bodyName - наименование тела
	* @param bodyId - идентификатор тела
	* @return признак успешной (true) или неуспешной (false) загрузки
	*/
	bool LoadFromXmlBodyElement
		(
			const TiXmlElement* bodyElement,
			BodyName& bodyName,
			BodyId& bodyId
		)
	{
		if (string(bodyElement->Value()) != "Body")
		{
			return false;
		}

		const char* bodyNameTmp;
	
		bodyNameTmp = bodyElement->Attribute("name");
		if (bodyNameTmp == nullptr)
		{
			return false;
		}
		bodyName = bodyNameTmp;

		const char* bodyIdAsString;
	
		bodyIdAsString = bodyElement->Attribute("id");
		if (bodyIdAsString == nullptr)
		{
			bodyId = 0;
		}
		else
		{
			bodyId = atoi(bodyIdAsString);
		}

		return true;
	}

	/**
	* «агрузить данные ??? из элемента xml-документа с именем MechanismOccurrenceData
	* @param dataElement - элемент xml-документа, из которого загружаютс€ данные
	* @param data - загруженные данные
	* @return признак успешной (true) или неуспешной (false) загрузки
	*/
	bool LoadFromXmlDataElement
		(
			const TiXmlElement* dataElement,
			MechanismOccurrenceData& data
		)
	{
		if (string(dataElement->Value()) != "MechanismOccurrenceData")
		{
			return false;
		}

		MechanismOccurrenceData dataTmp;
		const TiXmlElement* bodyElement;

		bodyElement = dataElement->FirstChildElement();
		while (bodyElement != nullptr)
		{
			BodyName bodyName;
			BodyId bodyId;

			if (!LoadFromXmlBodyElement(bodyElement, bodyName, bodyId))
			{
				return false;
			}
			if (dataTmp.find(bodyName) != dataTmp.end())
			{
				return false;
			}
			dataTmp.emplace(bodyName, bodyId);
			bodyElement = bodyElement->NextSiblingElement();
		}
		data = dataTmp;

		return true;
	}

	/**
	* «агрузить данные ??? из xml-документа с именем fileName
	* @param fileName - им€ файла, из которого загружаютс€ данные
	* @param data - загруженные данные
	* @return признак успешной (true) или неуспешной (false) загрузки
	*/
	bool LoadFromXmlDocument
		(
			const string& fileName,
			MechanismOccurrenceData& data
		)
	{
		TiXmlDocument xmlFile(fileName.c_str());

		if (!xmlFile.LoadFile())
		{
			exceptions::ThrowFileInvalidFormat(fileName);
		}

		const TiXmlElement* dataElement;

		dataElement = xmlFile.RootElement();
		if (dataElement == nullptr)
		{
			return false;
		}
	
		return LoadFromXmlDataElement(dataElement, data);
	}

	/**
	* «агрузить данные ??? из файла с именем fileName
	* @param fileName - им€ файла, из которого загружаютс€ данные
	* @param data - загруженные данные
	* @return признак успешной (true) или неуспешной (false) загрузки
	*/
	bool Load
		(
			const string& fileName,
			MechanismOccurrenceData& data
		)
	{
		try
		{
			if (!fs::FileExist(fileName))
			{
				exceptions::ThrowFileNotFound(fileName);
			}
			if (!LoadFromXmlDocument(fileName, data))
			{
				exceptions::ThrowFileNotOpened(fileName);
			}
		}
		catch (const exceptions::CoreException& exception)
		{
			exceptions::Output(exception);

			return false;
		}

		return true;
	}

}