#ifndef FILE_STREAM_H

#define FILE_STREAM_H


//модуль для определения интерфейсов классов для работы с файлами


#include "Types.h"
#include "Buffer.h"

#include <string>


using std::string;


//размер буфера потока ввода/вывода для файлов UNFORMATTED MS PowerStation
const int PowerStationBufferSize = 0x80;


//перечисление используемых манипуляторов
enum class Manipulator
{
	EndBlock //отправляет сохраненный в буфере фрагмент в файл (UF MS PS only)
};


//класс, позволяющий читать из и записывать в файл
class FileStream
{
public:

	//перечисления

	//режим открытия файла
	enum class OpenMode
	{
		NotOpen, //пока/уже не открыт
		ReadOnly, //открыт только для чтения
		WriteOnly, //открыт только для записи
		ReadWrite //открыт как для чтения, так и для записи
	};

	//режим сдвига маркера
	enum class SeekMode
	{
		FromCurrent, //сдвиг от текущего положения маркера
		FromBegin, //сдвиг от начала файла
		FromEnd //сдвиг от конца файла
	};


	// Конструкторы и деструкторы

	//конструктор
	FileStream
		(
			const string& fileName = string()
		);

	//деструктор
	virtual
	~FileStream();


	// Селекторы

	//получение режима открытия файла
	OpenMode GetOpenMode() const;

	//получение текущей позиции маркера
	Integer GetPosition() const;

	//возвращает истину, если файл открыт
	bool IsOpen() const;

	//возвращает истину, если достигнут конец файла
	bool IsEof();

	//получение имени файла
	const string& GetFileName() const;


	// Модификаторы

	//установка нового имени файла
	void SetFileName
		(
			const string& fileName
		);


	// Функции открытия/закрытия

	//открытие
	virtual bool Open
		(
			OpenMode mode
		);

	//закрытие
	virtual void Close();


	// Дополнительные функции

	//функция перемещения маркера
	bool Seek
		(
			Integer shift,
			SeekMode seekMode = SeekMode::FromBegin
		);


	// Функции ввода-вывода
	
	//чтение данных
	virtual
	Integer Read
		(
			void* data,
			Integer size
		);

	//запись данных
	virtual
	Integer Write
		(
			const void* data,
			Integer size
		);


	// Операторы

	// Операторы ввода-вывода

	// Для символа

	//оператор ввода
	friend FileStream& operator >>
		(
			FileStream& fileStream,
			Char& byte
		);

	//оператор вывода
	friend FileStream& operator <<
		(
			FileStream& fileStream,
			Char byte
		);


	// Для целых значений

	//оператор ввода
	friend FileStream& operator >>
		(
			FileStream& fileStream,
			Integer& number
		);

	//оператор вывода
	friend FileStream& operator <<
		(
			FileStream& fileStream,
			Integer number
		);


	// Для коротких целых значений
	
	//оператор ввода
	friend FileStream& operator >>
		(
			FileStream& fileStream,
			ShortInteger& number
		);

	//оператор вывода
	friend FileStream& operator <<
		(
			FileStream& fileStream,
			ShortInteger number
		);


	// Для вещественных значений

	//оператор ввода
	friend FileStream& operator >>
		(
			FileStream& fileStream,
			float& number
		);

	//оператор вывода
	friend FileStream& operator <<
		(
			FileStream& fileStream,
			float number
		);


	// Для вещественных значений двойной точности

	//оператор ввода
	friend FileStream& operator >>
		(
			FileStream& fileStream,
			double& number
		);

	//оператор вывода
	friend FileStream& operator <<
		(
			FileStream& fileStream,
			double number
		);


	// Для буфера

	//оператор ввода
	friend FileStream& operator >>
		(
			FileStream& fileStream,
			Buffer& buffer
		);

	//оператор вывода
	friend FileStream& operator <<
		(
			FileStream& fileStream,
			const Buffer& buffer
		);


	// Для манипуляторов

	//оператор вывода
	friend FileStream& operator <<
		(
			FileStream& fileStream,
			Manipulator manipulator
		);


protected:

	//поля

	FILE* _file; //указатель на файл

	string _fileName; //наименование файла

	OpenMode _openMode; //текущий режим открытия файла

	Integer _position; //позиция маркера


	//методы

	//вспомогательные методы

	//возвращает режим открытия в виде строки, необходимой для передачи в fopen
	static
	string GetOpenModeAsText
		(
			OpenMode mode
		);

	//функция применения манипулятора при выводе
	virtual
	bool ApplyManipulator
		(
			Manipulator manipulator
		);

};


//класс для работы с файлами типа UNFORMATTED языка Fortran (MS PowerStation)
class PowerStationFileStream
	:
		public FileStream
{
public:
	
	// Конструкторы и деструкторы

	PowerStationFileStream
		(
			const string& fileName = string()
		);

	//деструктор
	virtual
	~PowerStationFileStream() override;


	// Функции открытия/закрытия

	//открытие
	virtual
	bool Open
		(
			OpenMode mode
		)	override;

	//закрытие
	virtual
	void Close() override;


	// Функции ввода-вывода
	
	//чтение данных
	virtual
	Integer Read
		(
			void* data,
			Integer size
		)	override;

	//запись данных
	virtual
	Integer Write
		(
			const void* data,
			Integer size
		)	override;


protected:

	//поля

	Buffer _buffer; //буфер


	//методы

	//функция применения манипулятора при выводе
	virtual
	bool ApplyManipulator
		(
			Manipulator manipulator
		);

	//чтение блока
	void ReadBlock();

	//запись блока
	void WriteBlock
		(
			bool isFullBlock
		);

};


#endif // FILE_STREAM_H