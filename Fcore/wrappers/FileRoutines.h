#ifndef FILE_ROUTINES_H

#define FILE_ROUTINES_H


#include "fcUnistd.h"

#include <string>
#include <sstream>
#include <fstream>
#include <vector>


using std::string;
using std::stringstream;
using std::ifstream;
using std::ofstream;
using std::ios_base;
using std::vector;


#define NIX_PATH_SEPARATOR '/'
#define EXT_SEPARATOR '.'


namespace fs
{
	/** Создание bak-копии файла
	* @param filename - имя файла
	*/
	void BackupFile(const char* filename);

	/** Добавление линии к файлу
	*/
	void AppendLineToFile(const char* fname, stringstream& ss);

	/** Удаление /\ \/, замена \\ на /
	*/
	string NormalizePath(string path);

	/** Инициализация переменной окружения FRUND
	* @param - входной поток конфигурационного файла
	* @return путь к ФРУНДу
	*/
	string InitializeFromConfigStream(ifstream& ifs);

	/** Объединяет имя файла и имя каталога
	* @param dir - имя каталога
	* @param file - имя файла
	* @param string - полный путь к файлу
	*/
	string CombinePath(const string& dir, const string& file);

	/** Объединяет имя файла и имя каталога из ENV_FRUND
	* @param dir - имя каталога
	* @param file - имя файла
	* @param string - полный путь к файлу
	*/
	string CombinePathEnv(const string& dir, const string& file);

	/** Отделяет файл от каталога в полном пути
	* @param path - полный путь
	* @return имя файла
	*/
	string SplitFileFromPath(const string& inpath);

	/** Отделяет имя файла от расширения
	* @param path - имя файла
	* @return имя файла
	*/
	string SplitFileFromExt(const string& inpath);

	/** Отделяет каталог в полном пути
	* @param path - полный путь
	* @return путь к каталогу без имени файла
	*/
	string SplitDirFromPath(const string& inpath);

	/** Присоединение имени файла к каталогу с рабочей программой
	* @param fileName - имя файла
	* @param argvPath - путь из argv[0], используется только для NIX
	* @return путь из катлога с программой и имени передаваемого файла
	*/
	string GetFullPath
		(
			const string& fileName,
			const string& argvPath = string()
		);

	/** Сроит абсолютный путь по относительному
	* @param path - путь относительно текущего каталога
	* @return абсолютный путь
	*/
	//string GetRealPath(const string& path);

	/**
	* Чтение в строковый поток строки из файлового потока
	* @param ifs - файловый поток
	* @param line - строковый поток
	* @return true, если строка не пустая
	*/
	bool ReadLine(ifstream& ifs, stringstream& line);
	
	/**
	* Чтение в строковый поток строки из файлового потока без исключения
	* @param ifs - файловый поток
	* @param line - строковый поток
	* @return true, если строка не пустая и не конец файла
	*/
	bool ReadLineNoException(ifstream& ifs, stringstream& line);

	/**
	* Чтение строки из файлового потока
	* @param ifs - файловый поток
	* @param line - строка
	*/
	void ReadLineString(ifstream& ifs, string& line);

	/**
	* Чтение строки из файлового потока
	* @param ifs - файловый поток
	* @param line - строка
	*/
	void ReadLineString(ifstream& ifs, char* line);
	

	/**
	* Чтение строки из файлового потока
	* @param ifs - файловый поток
	* @param line - строка
	* @return false if EOF
	*/
	bool ReadLineStringNoException(ifstream& ifs, string& line);
        bool ReadLineStringNoException2(ifstream& ifs, string& line);

	/**
	* Чтение произвольного типа из строки файлового потока
	* @param ifs - файловый поток
	* @param out - считываемое значение
	*/
	template <class T>
	void ReadLineValue(ifstream& ifs, T& out)
	{
		stringstream line;
		if(ReadLine(ifs, line))
			line >> out;
	}

	/**
	* Чтение произвольного типа из строки файлового потока
	* @param ifs - файловый поток
	* @param out - считываемое значение
	*/
	template <class T>
	bool ReadLineValueNoException(ifstream& ifs, T& out)
	{
		stringstream line;
		bool result = ReadLineNoException(ifs, line);
		if(result)
			line >> out;
		return result;
	}

	/**
	* Чтение произвольного типа из строки файлового потока
	* @param ifs - файловый поток
	* @param tag - тэг
	* @param out - считываемое значение
	*/
	template <class T>
	bool ReadLineTagValueNoException(ifstream& ifs, string& tag, T& out)
	{
		stringstream line;
		bool res = ReadLineNoException(ifs, line);
		if(res)
			line >> tag >> out;
		return res;
	}

	/**
	* Чтение произвольного типа из строки файлового потока
	* @param ifs - файловый поток
	* @param out - считываемое значение
	*/
	template <class T>
	void WriteCommentedLine(ofstream& ofs, const T& value, const string& comment)
	{
		ofs << value << " \t// " << comment << std::endl;
	}

	//void WriteCommentedLine(ofstream& ofs, const double& value, const string& comment)
	//{
	//	ofs << fixed << value << " \t! " << comment << std::endl;
	//}

	/**
	* Чтение из строки файлового потока пары ключ-значение
	* @param ifs - файловый поток
	* @param out - считываемое значение
	* @return false, если еонец файла или пустая строка
	*/
	bool ReadLineKeyValue(ifstream& ifs, string& key, string& value);

	/** Проверка существования файла
	* @param path - путь у файлу
	* @return 0, если не существует
	*/
	int FileExist(const char* path);

	/** Проверка файла-флага
	* @param path - путь у файлу
	* @return true, если файл есть и в нем не 0
	*/
	bool CheckFileFlag(const string& path);


	/** Проверка расширения файла
	* @param path - путь у файлу
	* @param ext - проверяемое расширение
	* @return 0, если не существует
	*/
	int CheckExtension(const char* path, const char* ext);

	/** Удаляет файлы по маске типа *.ext
	* @param path - маска или путь к файлу
	*/
	void DeleteFiles(const char* path);

	/** Удаляет файлы по маске типа *.ext с постфиксом
	* @param mask - маска или путь к файлу
	* @param prefix - имя без постфикса
	* @param postfixSeparator - разделитель постфикса
	* @param postfixLength - ограничение на длину постфикса
	*/
	void DeleteFilesWithPostfix(const string& mask, const string& prefix, char postfixSeparator, int postfixLength);

	/** Копирует файл
	* @param srcPath - файл-источник
	* @param destPath - файл, в который копируется
	*/
	void CopySingleFile(const char* srcPath, const char* destPath);

	/** Переименовывает файл
	* @param srcPath - файл-источник
	* @param destPath - файл, в который переименовывается
	*/
	void RenameFile(const char* srcPath, const char* destPath);

	/** Объединение файлов
	* @param path1 - путь к первому файлу
	* @param path2 - путь ко второму файлу
	* @param destPath - файл, в который объединять
	*/
	void MergeFiles(const char* path1, const char* path2, const char destPath);

	/** Объединение файлов
	* @param path1 - путь к первому файлу, к которому добавляем
	* @param path2 - путь ко второму файлу, который дописывается в первый
	*/
	void AppendFile(const char* path1, const char* path2);

	/** Выполнение файла
	* @param process - запускаемый процесс
	* @param cmd - параметры коммандной строки
	* @return код возврата (exit code)
	*/
	int Exec
		(
			const string& process,
			const string& cmd
		);

#ifdef NIX
    /** Получени списка имен файлов
     *  @param ext - расширение файла
     *  @return список файлов
     */
    vector<string> GetListFile(const char* ext);

    /** Замена фильтра файлов на список файлов
     * @param keys список параметров командной строки
     * @param ext расширение файлов
     */
    void FormatCMD(string &keys, const char* ext);
#endif


	/** Сроит абсолютный путь по относительному
	* @param path - путь относительно текущего каталога
	* @return абсолютный путь
	*/
	string GetRealPath(const string& path);

	/** Проверка наличия каталога
	* @param path - путь к каталогу
	* @return true, если каталог существует
	*/
	bool ExistDir(const string& path);
	
	/** Получить стандартный разделитель каталогов
	* @return стандартный разделитель каталогов
	*/
	string DirSeparator();		

	/** Текущий каталог
	* @return текущий каталог
	*/
	string GetCurrentDir();

	/** Перейти в каталог
	* @param path - каталог
	*/
	void SetCurrentDir(const string& path);

	/** Перейти в каталог
	* @param path - каталог
	* @param oldPath - текущий каталог
	*/
	void SetCurrentDir(const string& path, string& oldPath);

	/** Создать каталог
	* @param каталог
	*/
	void MakeDir(const string& path);

	/** Удалить файл
	* @param path - путь к файлу
	*/
	void DeleteFile(const string& path);

	/** Проверить наличие файла
	* @param path - путь к файлу
	* @return - true, если файл существует
	*/
	bool FileExist(const string& path);

	/** Копирует файл в каталог
	* @param file - имя файла
	* @param folder - каталог
	*/
	void CopyFileToFolder(const string& file, const string& folder);

	/** Копирует файлы в каталог (допускаются записи вида *.ext)
	* @param files - список файлов
	* @param folder - каталог
	*/
	void CopyFilesToFolder(const vector<string>& files, const string& folder);


	/** Копирует файлы из каталога в каталог
	* @param inFolder - исходный каталог
	* @param outFolder - каталог, в который копировать
	* если его нет, он создается
	*/
	void CopyFolder(const string& _inFolder, const string& _outFolder);

	/** Удаляет файлы из каталога
	* @param folder - каталог
	*/
	void CleanFolder(const string& _folder);

	string AppendToFileName(const string& fileName, const string& str);
}


#endif // FILE_ROUTINES_H