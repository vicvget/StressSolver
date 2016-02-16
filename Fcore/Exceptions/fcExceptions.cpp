#include "fcExceptions.h"

#include <sstream>
#include <iostream>

using std::stringstream;
using std::cout;
using std::endl;


char FileExceptionErrors[][256] = 
{
	/* 0 */ "file is not found",
	/* 1 */ "file is not opened",
	/* 2 */ "file is not saved",
	/* 3 */ "invalid file format",
	/* 4 */ "invalid file extension",
	/* 5 */ "folder is not found"
};

namespace exceptions
{

	void ThrowMessage
		(
			const string& message
		)
	{
		throw MessageException(message);
	}

	void ThrowNoShellCommands
		(
			const string& env
		)
	{
		throw MessageException(string("Shell Commands in config file not found: ") + env);	
	}

	void ThrowFolderNotFound
		(
			const string& path
		)
	{
		throw FileException(path, FileExceptionErrors[5]);
	}

	void ThrowFileNotFound
		(
			const string& path
		)
	{
		throw FileException(path, FileExceptionErrors[0]);
	}

	void ThrowFileNotOpened
		(
			const string& path
		)
	{
		throw FileException(path, FileExceptionErrors[1]);
	}

	void ThrowFileNotSaved
		(
			const string& path
		)
	{
		throw FileException(path, FileExceptionErrors[2]);
	}

	void ThrowFileInvalidExtension
		(
			const string& path
		)
	{
		throw FileException(path, FileExceptionErrors[4]);
	}

	void ThrowFileInvalidFormat
		(
			const string& path
		)
	{
		throw FileException(path, FileExceptionErrors[3]);
	}

	void ThrowInvalidWildcard
		(
			const string& wildcard
		)
	{
		throw MessageException(string("Wildcard can be only *.ext: ") + wildcard);
	}

	void ThrowNoStringInFile
		(
			const string& path,
			const string& str
		)
	{
		throw FileException(path, string("string ") + str + " not found");
	}

	void ThrowNoSymbolInFile
		(
			const string& path,
			const string& str
		)
	{
		throw FileException(path, string("symbol ") + str + " not found");
	}

	void ThrowDllNotLoaded
		(
			const string& filename
		)
	{
		throw MessageException(string("Error loading dll: ") + filename);
	}

	void ThrowRunTimeError
		(
			const string& module,
			int errorCode
		)
	{
		throw RunTimeException(errorCode, string("Runtime error in module ") + module);
	}

	void ThrowRunTimeError
		(
			const string& module,
			const string& errorMessage,
			int errorCode
		)
	{
		throw RunTimeException
			(
				errorCode,
				string("Runtime error in module ") + module + string(": ") + errorMessage
			);
	}

	void ThrowObjectNotInitialized
		(
			const string& module,
			const string& object
		)
	{
		throw RunTimeException(-1, string("Object ") + object + " is not initialized in module " + module);
	}

	void ThrowInvalidCmd
		(
			const string& param
		)
	{
		throw CmdException(string("Invalid cmd\n SYNTAX:") + param);
	}

	void ThrowInvalidCmdParam
		(
			const string& param
		)
	{
		throw CmdException(string("Invalid cmd parameter ") + param);
	}

	void ThrowInvalidCmdParam
		(
			const string& key,
			const string& param
		)
	{
		throw CmdException(string("Invalid cmd parameter ") + param + " in key " + key);
	}

	void ThrowInvalidCmdOrder
		(
			const string& param1,
			const string& param2
		)
	{
		throw CmdException(string("Invalid cmd order: ") + param1 + " can't be after " + param2);
	}

	void ThrowTagNotFound
		(
			const string& tag
		)
	{
		throw MessageException(string("Tag not found: ") + tag);
	}

	void ThrowTagNotResolved
		(
			const string& tag
		)
	{
		throw MessageException(string("Tag not resolved: ") + tag);
	}
	
	void ThrowNoEnviromentVariable
		(
			const string& env
		)
	{
		throw MessageException(string("Enviroment variable not found: ") + env);
	}

	void ThrowNoElementInCollection
		(
			const string& element,
			const string& collectionName
		)
	{
		string message = " element " + string(element) + " not found ";

		throw CollectionException(collectionName, message);
	}
	
	void ThrowCollectionOutOfBounds
		(
			const string& collectionName
		)
	{
		throw CollectionException(collectionName, " is out of bounds ");
	}

	void ThrowModelNotInitialized()
	{
		throw MessageException("Model is not initialized");
	}

	/** Исключение о ненайденном теле
	* @param bodyNumber - номер тела
	*/
	void ThrowBodyNotFoundException
		(
			int bodyNumber
		)
	{
		stringstream stream;

		stream << "Body not found: " << bodyNumber;

		throw MessageException(stream.str());
	}

#pragma region CoreExceptionClass

	// virtual
	CoreException::~CoreException()
	{
	}

#pragma endregion

#pragma region MessageExceptionClass

	MessageException::MessageException
		(
			const string& message
		)
		:
			_message(message)
	{
	}

	// virtual
	string MessageException::ToString() const // override
	{
		return _message;
	}

#pragma endregion

#pragma region FileExceptionClass

	FileException::FileException
		(
			const string& path,
			const string& message
		)
		:
			MessageException(message),
			_path(path)
	{
	}

	// virtual
	string FileException::ToString() const // override
	{
		return _path + ':' + _message;
	}

#pragma endregion

#pragma region CollectionExceptionClass

	CollectionException::CollectionException
		(
			const string& collectionName,
			const string& message
		)
		:
			MessageException(message),
			_collectionName(collectionName)
	{
	}

	// virtual
	string CollectionException::ToString() const // override
	{
		return _collectionName + _message;
	}

#pragma endregion

#pragma region RunTimeExceptionClass

	RunTimeException::RunTimeException
		(
			int errorCode,
			const string& message
		)
		:
			MessageException(message),
			_errorCode(errorCode)
	{
	}

	// virtual
	string RunTimeException::ToString() const // override
	{
		stringstream stream;

		stream << "Error code " << _errorCode << ": " << _message << endl;

		return stream.str();
	}

#pragma endregion

#pragma region CmdExceptionClass

	CmdException::CmdException
		(
			const string& message
		)
		:
			MessageException(message)
	{
	}

#pragma endregion

	void Output
		(
			const CoreException& exp
		)
	{
		//TODO: возможно стоит изменить на cerr
		cout << "EXCEPTION: " << exp.ToString() << endl;
	}

}