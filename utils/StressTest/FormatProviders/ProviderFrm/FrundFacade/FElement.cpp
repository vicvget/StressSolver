#include "FElement.h"

void FElement::DefaultLoad(const char* path)
{
	string strPath(path);
	int id1 = strPath.find_first_of(_ext);
	_fileId = strPath.substr(0,id1);

	ifstream stream(path);
	if(!stream.is_open()) 
		exceptions::ThrowFileNotOpened(path);
	else
	{
		if(!Load(stream)) 
			exceptions::ThrowFileInvalidFormat(path);
		stream.close();
	}
}

void FElement::Save(const char* path) const
{
	ofstream ofs(path);
	if(!ofs.is_open()) 
		exceptions::ThrowFileNotOpened(path);
	else
	{
		Save(ofs);
	}
}

void FElement::Load(const char* path)
{
	DefaultLoad(path);
}

/**
 Функция, сохраняющая элемент ФРУНДа в файл с именем из fileId
*/
void FElement::SaveByIndex() const
{
	string path = _fileId + _ext;
	ofstream stream(path);
	if(!stream.is_open()) 
		exceptions::ThrowFileNotOpened(path);
	else
	{
		Save(stream);
	}
}

bool FElement::LoadByIndex()
{
	string path = _fileId + _ext;
	ifstream ifs(path);
	if(ifs.is_open())
	{
		bool res = (bool)Load(ifs);
		ifs.close();
		return res;
	}
	return false;	
}

/** Запись в файловый поток структуры модели в зависимости от шага моделирования
* не чисто виртуальный из-за MPH
*
* @param stream - поток
* @param stage - шаг моделирования 0 = генерация, 1 = расчет
*/
// virtual
void FElement::Output
	(
		ofstream& /* stream */,
		int /* stage */
	)	const
{
}


FElement::FElement
	(
		const FElement& src
	)
	:
		_name(src._name),
		_fileId(src._fileId),
		_ext(src._ext),
		_id(src._id),
		_number(src._number)
{

}

FElement::FElement()
{
	Init("", EXT_MODEL_ELEMENT, 0, 1);
}

FElement::FElement(int index)
{
	Init("", EXT_MODEL_ELEMENT, index, index + 1);
}

FElement::FElement(const string& fileId)
{
	Init(fileId, EXT_MODEL_ELEMENT, 0, 1);
}

void FElement::Init
	(
		const string& fileId,
		const string& ext,
		int id,
		int number
	)
{
	_fileId = fileId;
	_ext = ext;
	_id = id;
	_number = number;
}

void FElement::CheckFileId() const
{
	if (_fileId.empty())
		exceptions::ThrowMessage("Empty fileId");
}

string FElement::GetPath(const string& ext) const
{
	CheckFileId();
	return _fileId + ext;
}

string FElement::GetPath() const
{
	CheckFileId();
	return _fileId+_ext;
}
