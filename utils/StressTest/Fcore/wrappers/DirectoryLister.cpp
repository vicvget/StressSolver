#include "DirectoryLister.h"

#include "FileRoutines.h"

#ifdef NIX
#include <glob.h>
#include <sys/stat.h>
#else
#include <io.h>
#endif

#include <algorithm>
#include <cctype>


namespace fs
{

/** Конструктор составляет список директории
* @param path - путь к каталогу
*/
DirectoryLister::DirectoryLister()
{
	ReadCurrentDir();
}

/** Конструктор составляет список директории
* @param path - путь к каталогу
*/
DirectoryLister::DirectoryLister(const string& path)
{
	string currentDir = fs::GetCurrentDir();
	fs::SetCurrentDir(path);
	ReadCurrentDir();
	fs::SetCurrentDir(currentDir);
}

DirectoryLister::DirectoryLister(const string& path, const string& extension)
{
	string currentDir = fs::GetCurrentDir();
	fs::SetCurrentDir(path);
	ReadCurrentDir(extension);
	fs::SetCurrentDir(currentDir);
}

void DirectoryLister::ReadCurrentDir()
{
 #ifdef NIX

    string wildcard ="*";
    glob_t globbuf;

    if(!glob(wildcard.c_str(), 0, NULL, &globbuf)) {
                size_t i=0;
                for(;i<globbuf.gl_pathc;++i) {
                    struct stat st;
                    if((!stat(globbuf.gl_pathv[i],&st))&&( S_ISDIR(st.st_mode)))
                  //      cout << globbuf.gl_pathv[i] << " | dir" <<endl;
                        folders.push_back(string(globbuf.gl_pathv[i]));
                    else
                        //cout << globbuf.gl_pathv[i] << " | file" <<endl;
                        files.push_back(string(globbuf.gl_pathv[i]));
                }
        }
    globfree(&globbuf);
    #else
	struct _finddata_t c_file;
	intptr_t hFile;
	if( (hFile = _findfirst( "*", &c_file )) != -1L )
	{
		do 
		{
			if( c_file.attrib & _A_SUBDIR )
			{
				bool isSysdir = false;
				if(c_file.name[0] == '.')
				{
					if(c_file.name[1] == '.')
						isSysdir = (c_file.name[2] == '\0') ? true : false;
					else
						isSysdir = (c_file.name[1] == '\0') ? true : false;
				}
				if(!isSysdir)
					folders.push_back(c_file.name);
			}
			else
				files.push_back(c_file.name);
		} 
		while( _findnext( hFile, &c_file ) == 0 );
		_findclose( hFile );
	}
#endif
}

void DirectoryLister::ReadCurrentDir
	(
		const string& extension
	)
{
#ifdef NIX
    glob_t globbuf;

    if(!glob(extension.c_str(), 0, NULL, &globbuf)) {
                size_t i=0;
                for(;i<globbuf.gl_pathc;++i) {
                    struct stat st;
                    if((!stat(globbuf.gl_pathv[i],&st))&&( S_ISDIR(st.st_mode)))
                  //      cout << globbuf.gl_pathv[i] << " | dir" <<endl;
                        folders.push_back(string(globbuf.gl_pathv[i]));
                    else
                        //cout << globbuf.gl_pathv[i] << " | file" <<endl;
                        files.push_back(string(globbuf.gl_pathv[i]));
                }
        }
    globfree(&globbuf);
#else
struct _finddata_t c_file;
	intptr_t hFile;
	if( (hFile = _findfirst(extension.c_str(), &c_file)) != -1L )
	{
		do 
		{
			if( c_file.attrib & _A_SUBDIR )
			{
				bool isSysdir = false;
				if(c_file.name[0] == '.')
				{
					if(c_file.name[1] == '.')
						isSysdir = (c_file.name[2] == '\0') ? true : false;
					else
						isSysdir = (c_file.name[1] == '\0') ? true : false;
				}
				if(!isSysdir)
					folders.push_back(c_file.name);
			}
			else
				files.push_back(c_file.name);
		} 
		while( _findnext( hFile, &c_file ) == 0 );
		_findclose( hFile );
	}
#endif
}

/** Проверка наличия файла
* @param fileName - имя файла
* @return true, если существует
*/
bool DirectoryLister::HasFile(const string& fileName)
{
	return find(files.begin(), files.end(), fileName) != files.end();
}

/** Проверка наличия каталога
* @param fileName - имя каталога
* @return true, если существует
*/
bool DirectoryLister::HasFolder(const string& folderName)
{
	return find(folders.begin(), folders.end(), folderName) != folders.end();
}

DirectoryLister::~DirectoryLister()
{
}

void DirectoryLister::SetAllToLowerCase()
{
#ifndef GNUCPP
	for (size_t i = 0; i < files.size(); ++i)
	{
		string t = files[i];
		std::transform(t.begin(), t.end(), t.begin(), std::tolower);
		files[i] = t;
	}
#endif
}

// Установка  пути к первому файлу в каталоге
size_t DirectoryLister::Init()
{
	itPath = files.begin();

	return files.size();
}

// Возврат и переход к следующему
int DirectoryLister::Next(string& path)
{
    if(IsLast())
        return 0;
    path = *itPath;
    ++itPath;
    return 1;
}

// Признак конца
int DirectoryLister::IsLast()
{
    return itPath == files.end() ? 1 : 0;
}
}