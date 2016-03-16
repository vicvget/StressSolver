#ifndef FELEMENT_H

#define FELEMENT_H


#include "BasicTypes.h"

/**
* Базовый класс для считывания элемента модели ФРУНД
* @author Getmanskiy Victor
*/
class FElement
{
protected:
	// имя
	string _name;
	// индекс файла ( в нотации ФРУНД 100000...)
	string _fileId;
	// расширение
	string _ext;
	// числовой индекс (порядковый номер тела, но для тела почему-то +101 ???)
	int _id;
	// номер
	int _number;
public:

	/** Считывание из файла
	* @param path - путь к файлу
	*/
	virtual void Load(const char* path);

	/** Считывание из файла по умолчанию
	* @param path - путь к файлу
	*/
	void DefaultLoad(const char* path);

	/** Считывание из файлового потока
	* @param stream - поток
	* @return true, Если успешное считывание
	*/
	virtual
	bool Load
		(
			ifstream& stream
		)
		= 0;
	
	/** Запись в файловый поток
	* @param stream - поток
	*/
	virtual
	void Save
		(
			ofstream& ofs
		)	const
		= 0;
		
	/** Запись в файл
	* @param path - файл
	*/
	void Save(const char* path) const;

	/** Запись в файл по файловому индексу
	*/
	void SaveByIndex() const;

	/** Загрузка по файловому индексу
	* @returns true, если загружено
	*/
	bool LoadByIndex();


	/** Запись в файловый поток структуры модели в зависимости от шага моделирования
	* не чисто виртуальный из-за MPH
	*
	* @param stream - поток
	* @param stage - шаг моделирования 0 = генерация, 1 = расчет
	*/
	virtual
	void Output
		(
			ofstream& stream,
			int stage
		)	const;

	FElement();
	FElement(int index);
	FElement(const string& fileIndex);
	FElement(const FElement& src);

	virtual ~FElement() {}

	void Init(const string& fileId, const string& ext, int id, int number);

	void CheckFileId() const;
	string GetPath() const;
	string GetPath(const string& ext) const;

#pragma region GettersSettres // by Sergeev

	string GetFileIndex() const {return _fileId;}

	int Number() const { return _number; }
	void Number(int val) { _number = val; }

	string Name() const { return _name; }
	void Name(string val) { _name = val; }

	string FileId() const { return _fileId; }
	void FileId(string val) { _fileId = val; }

	int Id() const { return _id; }
	void Id(int val) { _id = val; }

	void SetExt(const string& ext) {_ext = ext;}

#pragma  endregion

};


#endif // FELEMENT_H