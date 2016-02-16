#ifndef fcExceptionsH

#define fcExceptionsH


#ifdef GNUCPP

#include "../wrappers/fcUnistd.h"

#else

#include <string>

using std::string;

#endif


/**
* Исключения фасада системы
*
* @author Victor Getmanskiy
*/
namespace exceptions
{

	/** Исключение об ошибке
	* @param message - сообщение
	*/
	void ThrowMessage
		(
			const string& message
		);

	/** Исключение об отсутствии набора команд компилятора
	* @params shellCommands -  имя набора команд компилятора
	*/
	void ThrowNoShellCommands
		(
			const string& shellCommands
		);

	/** Исключение об отсутствии каталога
	* @param path - путь к каталогу
	*/
	void ThrowFolderNotFound
		(
			const string& path
		);

	/** Исключение об отсутствии файла
	* @param path - путь к файлу
	*/
	void ThrowFileNotFound
		(
			const string& path
		);

	/** Исключение об ошибке открытия файла
	* @param path - путь к файлу
	*/
	void ThrowFileNotOpened
		(
			const string& path
		);

	/** Исключение об ошибке сохранения файла
	* @param path - путь к файлу
	*/
	void ThrowFileNotSaved
		(
			const string& path
		);

	/** Исключение о некорректном расширении файла
	* @param path - путь к файлу
	*/
	void ThrowFileInvalidExtension
		(
			const string& path
		);

	/** Исключение о некорректном формате файла
	* @param path - путь к файлу
	*/
	void ThrowFileInvalidFormat
		(
			const string& path = "undefined"
		);

	/** Исключение о неверной маске файла
	* допускается только маска в формате *.ext
	* @param wildcard - маска
	*/
	void ThrowInvalidWildcard
		(
			const string& wildcard
		);

	/** Исключение о ненайденной строке в файле
	* @param path - путь к файлу
	* @param str - искомая строка
	*/
	void ThrowNoStringInFile
		(
			const string& path,
			const string& str
		);

	/** Исключение о ненайденном символе в DLL файле
	* @param path - путь к файлу
	* @param symbol - символ
	*/
	void ThrowNoSymbolInFile
		(
			const string& path,
			const string& symbol
		);

	/** Исключение о незагруженной Dll
	*/
	void ThrowDllNotLoaded
		(
			const string& filename
		);

	/** Исключение об ошибке исполнения
	* @param module - выполняемый модуль
	* @param errorCode - код возврата
	*/
	void ThrowRunTimeError
		(
			const string& module,
			int errorCode
		);

	/** Исключение об ошибке исполнения
	* @param module - выполняемый модуль
	* @param errorMessage - сообщение об ошибке
	* @param errorCode - код возврата
	*/
	void ThrowRunTimeError
		(
			const string& module,
			const string& errorMessage,
			int errorCode
		);

	/** Исключение о неинициализированном объекте
	* @param module - выполняемый модуль
	* @param object - объект
	*/
	void ThrowObjectNotInitialized
		(
			const string& module,
			const string& object
		);

	/** Исключение об ошибке командной строки
	* @param param - коректный формат строки для вывода сообщения
	*/
	void ThrowInvalidCmd
		(
			const string& param
		);

	/** Исключение об ошибке командного параметра
	* @param param - параметр ключа
	*/
	void ThrowInvalidCmdParam
		(
			const string& param
		);

	/** Исключение об ошибке командного параметра в ключе
	* @param key - ключ
	* @param param - параметр ключа
	*/
	void ThrowInvalidCmdParam
		(
			const string& key,
			const string& param
		);

	/** Исключение об ошибке порядка параметров
	* @param param1 - первый параметтр
	* @param param2 - второй параметтр
	*/
	void ThrowInvalidCmdOrder
		(
			const string& param1,
			const string& param2
		);

	/** Исключение о ненайденном теге
	* @param tag - тэг
	*/
	void ThrowTagNotFound
		(
			const string& tag
		);

	/** Исключение о непрочитанном теге
	* @param tag - тэг
	*/
	void ThrowTagNotResolved
		(
			const string& tag
		);

	/** Исключение об отсутствии переменной окружения
	* @param env - имя переменной окружения
	*/
	void ThrowNoEnviromentVariable
		(
			const string& env
		);

	/** Исключение об отсутствии элемента в коллекции
	* @param element - элемент
	* @param collectionName - коллекция
	*/
	void ThrowNoElementInCollection
		(
			const string& element,
			const string& collectionName
		);

	/** Исключение о выходе за границы коллекции
	* @param collectionName - коллекция
	*/
	void ThrowCollectionOutOfBounds
		(
			const string& collectionName
		);
	
	/** Исключение о неинициализированной модели
	*/
	void ThrowModelNotInitialized();

	/** Исключение о ненайденном теле
	* @param bodyNumber - номер тела
	*/
	void ThrowBodyNotFoundException
		(
			int bodyNumber
		);

#pragma region ExceptionClasses

	/**
	* Базовый класс исключений
	*
	* @author Victor Getmanskiy
	*/
	class CoreException
	{
	public:

		virtual
		~CoreException();

		virtual
		string ToString() const = 0;

	};

	/**
	* Исключение в виде стркового сообщения
	*
	* @author Victor Getmanskiy
	*/
	class MessageException
		:
			public CoreException
	{
	public:

		MessageException
			(
				const string& message
			);

		virtual
		string ToString() const override;

	protected:

		string _message;

	};

	/**
	* Исключение при работе с файлами
	*
	* @author Victor Getmanskiy
	*/
	class FileException
		:
			public MessageException
	{
	public:

		FileException
			(
				const string& path,
				const string& message
			);

		virtual
		string ToString() const override;
	
	protected:

		string _path;

	};

	/**
	* Исключение в коллекции
	*
	* @author Victor Getmanskiy
	*/
	class CollectionException
		:
			public MessageException
	{
	public:

		CollectionException
			(
				const string& collectionName,
				const string& message
			);

		virtual
		string ToString() const override;

	protected:

		string _collectionName;

	};

	/**
	* Исключение при выполнении модуля
	*
	* @author Victor Getmanskiy
	*/
	class RunTimeException
		:
			public MessageException
	{
	public:

		RunTimeException
			(
				int errorCode,
				const string& message
			);

		virtual
		string ToString() const override;

	protected:

		int _errorCode;

	};
	
	/**
	* Исключение при обработке параметров коммандной строки
	*
	* @author Victor Getmanskiy
	*/
	class CmdException
		:
			public MessageException
	{
	public:

		CmdException
			(
				const string& message
			);

	};

#pragma endregion

	/** Вывод сообщения объекта-исключения
	* Вынесено в функцию для быстрого отключения вывода и для замены на cerr
	* @param path - путь к файлу
	*/
	void Output
		(
			const CoreException& exp
		);

}

#endif