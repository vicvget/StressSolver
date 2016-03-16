#ifndef DIRECTORY_LISTER_H

#define DIRECTORY_LISTER_H


#include <string>
#include <vector>


using std::string;
using std::vector;


namespace fs
{

	/** Класс для чтения каталога
	*
	* @author Victor Getmanskiy
	*/
	class DirectoryLister
	{
	public:
	
		DirectoryLister();

		/** Конструктор составляет список директории
		* @param path - путь к каталогу
		*/
		DirectoryLister
			(
				const string& path
			);

		DirectoryLister
			(
				const string& path,
				const string& extension
			);

		~DirectoryLister();
	
		/** Проверка наличия файла
		* @param fileName - имя файла
		* @return true, если существует
		*/
		bool HasFile
			(
				const string& fileName
			);

		/** Проверка наличия каталога
		* @param fileName - имя каталога
		* @return true, если существует
		*/
		bool HasFolder
			(
				const string& folderName
			);

		const string& GetFile
			(
				size_t id
			)	const
		{
			return files[id];
		}

		const string& GetFolder
			(
				size_t id
			)	const
		{
			return folders[id];
		}

		size_t GetFilesCount() const
		{
			return files.size();
		}

		size_t GetFoldersCount() const
		{
			return folders.size();
		}

		void SetAllToLowerCase();

		size_t Init();

		int Next
			(
				string& path
			);

		int IsLast();

	protected:

		// subfolders
		vector<string> folders;

		// files
		vector<string> files;

		vector<string>::iterator itPath;

		void ReadCurrentDir();

		/** Чтение текущего каталога
		* @param extension - расширение файлов
		*/
		void ReadCurrentDir
			(
				const string& extension
			);

	};

}


#endif // DIRECTORY_LISTER_H