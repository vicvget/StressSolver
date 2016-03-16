#include "MechanismOccurrenceData.h"

#include "../../../Fcore/Exceptions/fcExceptions.h"
#include "../../../Fcore/wrappers/FileRoutines.h"
#include "../../../ExternalModules/TinyXml/tinyxml.h"


namespace MechanismOccurrence
{

	/**
	* ��������� ������ � ���� �� �������� xml-��������� � ������ Body
	* @param bodyElement - ������� xml-���������, �� �������� ����������� ������
	* @param bodyName - ������������ ����
	* @param bodyId - ������������� ����
	* @return ������� �������� (true) ��� ���������� (false) ��������
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
	* ��������� ������ ??? �� �������� xml-��������� � ������ MechanismOccurrenceData
	* @param dataElement - ������� xml-���������, �� �������� ����������� ������
	* @param data - ����������� ������
	* @return ������� �������� (true) ��� ���������� (false) ��������
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
	* ��������� ������ ??? �� xml-��������� � ������ fileName
	* @param fileName - ��� �����, �� �������� ����������� ������
	* @param data - ����������� ������
	* @return ������� �������� (true) ��� ���������� (false) ��������
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
	* ��������� ������ ??? �� ����� � ������ fileName
	* @param fileName - ��� �����, �� �������� ����������� ������
	* @param data - ����������� ������
	* @return ������� �������� (true) ��� ���������� (false) ��������
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