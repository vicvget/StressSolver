#include "FileStream.h"

#include "../../AdditionalModules/AuxiliaryModules/Enumerations.h"

#include <cstdio>


//класс, позволяющий читать из и записывать в файл
//class FileStream

//методы

//конструкторы и деструкторы

//конструктор
FileStream::FileStream
	(
		const string& _fileName
	)
	:
		_file(nullptr),
		_fileName(_fileName),
		_openMode(OpenMode::NotOpen)
{
}

//деструктор
/*virtual*/
FileStream::~FileStream()
{
	Close();
}

//функции-селекторы

//получение режима открытия файла
FileStream::OpenMode FileStream::GetOpenMode() const
{
	return _openMode;
}

//получение текущей позиции маркера
Integer FileStream::GetPosition() const
{
	return _position;
}

//возвращает истину, если файл открыт
bool FileStream::IsOpen() const
{
	return _openMode != OpenMode::NotOpen;
}

//возвращает истину, если достигнут конец файла
bool FileStream::IsEof()
{
	if (_file == nullptr)
	{
		return false;
	}
	if (feof(_file))
	{
		return true;
	}
	else
	{
		return false;
	}
}

//получение имени файла
const string& FileStream::GetFileName() const
{
	return _fileName;
}

//установка нового имени файла
void FileStream::SetFileName
	(
		const string& fileName
	)
{
	if (_openMode == OpenMode::NotOpen)
	{
		_fileName = fileName;
	}
}

//функции открытия/закрытия

//открытие
/*virtual*/
bool FileStream::Open
	(
		OpenMode mode
	)
{
	if (_file == nullptr)
	{
		const string& openModeAsText = GetOpenModeAsText(mode);

		_file = fopen
			(
				_fileName.c_str(),
				openModeAsText.c_str()
			);
	}
	if (_file == nullptr)
	{
		return false;
	}
	else
	{
		_openMode = mode;

		return true;
	}
}

//закрытие
/*virtual*/
void FileStream::Close()
{
	if (_file != nullptr)
	{
		fclose(_file);
		_file = nullptr;
		_openMode = OpenMode::NotOpen;
	}
}

//дополнительные функции

//функция перемещения маркера
bool FileStream::Seek
	(
		Integer shift,
		SeekMode seekMode
	)
{
	Integer origin;

	switch (seekMode)
	{
		case SeekMode::FromCurrent:
			origin = SEEK_CUR;
			break;

		case SeekMode::FromBegin:
			origin = SEEK_SET;
			break;

		case SeekMode::FromEnd:
			origin = SEEK_END;
			break;

		//default:
			// unsupported
			//origin = ToIntegralType(seekMode);
	}
	if (fseek(_file, shift, origin) == 0)
	{
		return true;
	}
	else
	{
		return false;
	}
	
}

//функции ввода-вывода

//чтение данных
/*virtual*/
Integer FileStream::Read
	(
		void* data,
		Integer size
	)
{
	return fread
		(
			data,
			1,
			size,
			_file
		);
}

//запись данных
/*virtual*/
Integer FileStream::Write
	(
		const void* data,
		Integer size
	)
{
	return fwrite
		(
			(void*)data,
			1,
			size,
			_file
		);
}

//вспомогательные методы

//возвращает режим открытия в виде строки, необходимой для передачи в fopen
/*static*/
string FileStream::GetOpenModeAsText
	(
		OpenMode openMode
	)
{
	string openModeAsText;

	if (openMode != OpenMode::NotOpen)
	{
		switch (openMode)
		{
		case OpenMode::ReadOnly:
			openModeAsText = "r";
			break;

		case OpenMode::WriteOnly:
			openModeAsText = "w";
			break;

		case OpenMode::ReadWrite:
			openModeAsText = "r+";
			break;
		}
		openModeAsText += "b";
	}
	
	return openModeAsText;
}

//функция применения манипулятора при выводе
/*virtual*/
bool FileStream::ApplyManipulator
	(
		Manipulator /* manipulator */
	)
{
	return false;
}

//операторы

//операторы ввода-вывода

//для символа

//оператор ввода
/*friend*/
FileStream& operator >>
	(
		FileStream& fileStream,
		Char& byte
	)
{
	fileStream.Read(&byte, sizeof(byte));

	return fileStream;
}

//оператор вывода
/*friend*/
FileStream& operator <<
	(
		FileStream& fileStream,
		Char byte
	)
{
	fileStream.Write(&byte, sizeof(byte));

	return fileStream;
}

//для целых значений

//оператор ввода
/*friend*/
FileStream& operator >>
	(
		FileStream& fileStream,
		Integer& number
	)
{
	fileStream.Read(&number, sizeof(number));

	return fileStream;
}

//оператор вывода
/*friend*/
FileStream& operator <<
	(
		FileStream& fileStream,
		Integer number
	)
{
	fileStream.Write(&number, sizeof(number));

	return fileStream;
}

//для коротких целых значений

//оператор ввода
/*friend*/
FileStream& operator >>
	(
		FileStream& fileStream,
		ShortInteger& number
	)
{
	fileStream.Read(&number, sizeof(number));

	return fileStream;
}

//оператор вывода
/*friend*/
FileStream& operator <<
	(
		FileStream& fileStream,
		ShortInteger number
	)
{
	fileStream.Write(&number, sizeof(number));

	return fileStream;
}

//для вещественных значений

//оператор ввода
/*friend*/
FileStream& operator >>
	(
		FileStream& fileStream,
		float& number
	)
{
	fileStream.Read(&number, sizeof(number));

	return fileStream;
}

//оператор вывода
/*friend*/
FileStream& operator <<
	(
		FileStream& fileStream,
		float number
	)
{
	fileStream.Write(&number, sizeof(number));

	return fileStream;
}

//для вещественных значений двойной точности

//оператор ввода
/*friend*/
FileStream& operator >>
	(
		FileStream& fileStream,
		double& number
	)
{
	fileStream.Read(&number, sizeof(number));

	return fileStream;
}

//оператор вывода
/*friend*/
FileStream& operator <<
	(
		FileStream& fileStream,
		double number
	)
{
	fileStream.Write(&number, sizeof(number));

	return fileStream;
}

//для буфера

//оператор ввода
/*friend*/
FileStream& operator >>
	(
		FileStream& fileStream,
		Buffer& buffer
	)
{
	fileStream.Read(buffer.GetBuffer(), buffer.GetContentsSize());

	return fileStream;
}

//оператор вывода
/*friend*/
FileStream& operator <<
	(
		FileStream& fileStream,
		const Buffer& buffer
	)
{
	fileStream.Write(buffer.GetBuffer(), buffer.GetContentsSize());

	return fileStream;
}

//для манипуляторов

//оператор вывода
/*friend*/
FileStream& operator <<
	(
		FileStream& fileStream,
		Manipulator manipulator
	)
{
	fileStream.ApplyManipulator(manipulator);

	return fileStream;
}


//класс для работы с файлами типа UNFORMATTED языка Fortran (MS PowerStation)
//class PowerStationFileStream

//методы

//конструкторы и деструкторы

PowerStationFileStream::PowerStationFileStream
	(
		const string& fileName
	)
	:
		FileStream(fileName),
		_buffer(PowerStationBufferSize)
{
}

//деструктор
/*virtual*/
PowerStationFileStream::~PowerStationFileStream() // override
{
	Close();
}

//функции открытия/закрытия

//открытие
/*virtual*/
bool PowerStationFileStream::Open
	(
		OpenMode mode
	)	// override
{
	if (FileStream::Open(mode))
	{
		Char ch;

		switch (_openMode)
		{
		case OpenMode::ReadOnly:
			if (FileStream::Read(&ch, sizeof(ch)))
			{
				return true;
			}
			break;
		
		case OpenMode::WriteOnly:
			ch = 0x4B;
			if (FileStream::Write(&ch, sizeof(ch)))
			{
				return true;
			}
			break;
		}
	}
	return false;
}

//закрытие
/*virtual*/
void PowerStationFileStream::Close() // override
{
	Char ch;

	switch (_openMode)
	{
	case OpenMode::ReadOnly:
		break;
	
	case OpenMode::WriteOnly:
		if (!_buffer.IsEmpty())
		{
			WriteBlock(true);
		}
		ch = 0x82;
		FileStream::Write(&ch, sizeof(ch));
		break;
	}
	FileStream::Close();
	_buffer.Flush();
}

//функции ввода-вывода

//чтение данных
/*virtual*/
Integer PowerStationFileStream::Read
	(
		void* data,
		Integer size
	)	// override
{
	if (_buffer.Pop(data, size))
	{
		return size;
	}
	else
	{
		ReadBlock();
		if (_buffer.Pop(data, size))
		{
			return size;
		}
		else
		{
			return 0;
		}
	}
}

//запись данных
/*virtual*/
Integer PowerStationFileStream::Write
	(
		const void* data,
		Integer size
	)	// override
{
	if (_buffer.Push(data, size))
	{
		return size;
	}
	else
	{
		WriteBlock(false);
		if (_buffer.Push(data, size))
		{
			return size;
		}
		else
		{
			return 0;
		}
	}
}

//функция применения манипулятора при выводе
/*virtual*/
bool PowerStationFileStream::ApplyManipulator
	(
		Manipulator manipulator
	)
{
	if (!FileStream::ApplyManipulator(manipulator))
	{
		switch (manipulator )
		{
		case Manipulator::EndBlock:
			WriteBlock(true);
			break;

		default:

			return false;
		}
	}
	return true;
}

//чтение блока
void PowerStationFileStream::ReadBlock()
{
	Char ch;

	_buffer.Flush();
	FileStream::Read(&ch, sizeof(ch));
	if (ch == 0x82)
	{
		FileStream::Read(&ch, sizeof(ch));
		return;
	}
	if (ch == 0x81)
	{
		ch = 0x80;
	}
	_buffer.SetContentsLength(ch);

	FileStream::Read(_buffer.GetBuffer(), _buffer.GetContentsSize());
	FileStream::Read(&ch, sizeof(ch));
}

//запись блока
void PowerStationFileStream::WriteBlock
	(
		bool isFullBlock
	)
{
	Char ch;
	
	if (isFullBlock)
	{
		ch = _buffer.GetContentsSize();
	}
	else
	{
		ch = 0x81;
	}

	FileStream::Write(&ch, sizeof(ch));
	FileStream::Write(_buffer.GetBuffer(), _buffer.GetContentsSize());
	FileStream::Write(&ch, sizeof(ch));
	_buffer.Flush();
}
