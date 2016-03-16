#include "FileRoutines.h"

#include "DirectoryLister.h"
#include "../fcore.h"
#include "../Exceptions/fcExceptions.h"
#include "../wrappers/StringRoutines.h"

#ifdef GNUCPP

#	include <sys/stat.h>

#else

#	include <io.h>
#	include <direct.h>

#endif

#include <algorithm>
#include <stack>


using std::ios;
using std::stack;


#ifndef NIX
// для использования сокетов (Winsock2.h)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <limits.h>
#include <sys/wait.h>
#include <unistd.h>
#include <glob.h>
#endif


#ifdef DeleteFile

#define DELETE_FILE_DEFINE DeleteFile

#undef DeleteFile

#define WINDOWS_FUNCTION_NAME_DEFINED

#endif


namespace fs
{
	bool IsAbsolutePath(const string& path)
	{
#ifndef NIX
		if(path.find_first_of(":") != string::npos)
			return true;
#else
		// TODO: detect \dev ?
#endif
		return false;
	}
	
	string GetAbsolutePath(const string& path)
	{
		string res;
		res = IsAbsolutePath(path) ? path : fs::GetCurrentDir() + NIX_PATH_SEPARATOR + path;
		return res;
	}

	void BackupFile(const char* filename)
	{
		CopySingleFile(filename, (string(filename) + ".bak").c_str());
	}
	
	void AppendLineToFile(const char* fname, stringstream& ss)
	{
		ofstream ofs(fname, std::ios_base::app);
		char line[256];
		ss.getline(&line[0], 256);
		if(ofs.is_open())
		{
			ofs << line << std::endl;
			ofs.close();
		}
		else
		{
			exceptions::ThrowFileNotOpened(fname);
		}			
	}

	string InitializeFromConfigStream(ifstream& ifs)
	{
		// TODO: переделать !
		string sbuf, value;
		fs::ReadLineString(ifs,sbuf);
		if(sbuf.find("force_path=")==string::npos)
			exceptions::ThrowTagNotFound("force_path");
		value = sbuf.substr(sbuf.find("=")+1,sbuf.length());
		if(value != "0")
		{
			fs::ReadLineString(ifs,sbuf);
			string tag = "FRUND=";
			size_t id = sbuf.find(tag);
			if(id == string::npos)
				exceptions::ThrowTagNotFound("FRUND");
			value=sbuf.substr(tag.length()+id,sbuf.length());
		}
		else
		{
			if(getenv(ENV_FRUND) == NULL)
				exceptions::ThrowNoEnviromentVariable(ENV_FRUND);
			value.assign(getenv(ENV_FRUND));
		}	
		return value;
	}

	string NormalizePath(string path)
	{
		replace(path.begin(),path.end(),'\\',NIX_PATH_SEPARATOR);  // for Linux
		while(path.find("//",2) != string::npos)
		{
			path = path.substr(0,path.find("//",2)) + path.substr(path.find("//",2)+1, path.size());
		}

		size_t pos = 0;
		while ( (pos = path.find(";/",pos+2)) != string::npos)
		{
			path.insert(pos+1,"/");
		}
		
		if(path[path.size()-1] == NIX_PATH_SEPARATOR)
			path = path.substr(0,path.size()-1);
		return path;
	}

	/** Отделяет файл от каталога в полном пути
	* @param path - полный путь
	* @return имя файла
	*/
	string SplitFileFromPath(const string& inpath)
	{
		string path = NormalizePath(inpath);
		size_t pos = path.find_last_of(EXT_SEPARATOR);
		if(pos == string::npos)
		{
			// no file
			return string();
		}

		// inpath исправлено на path. Возможно, где-то работать не будет.
		pos = path.find_last_of(NIX_PATH_SEPARATOR);
		if(pos != string::npos)
		{
			// extract file
			return path.substr(pos+1, path.length());
		}

		// no dir, only file
		return inpath;
	}
	
	string SplitFileFromExt(const string& inpath)
	{
 		size_t pos = inpath.find_last_of(".");
		if(pos != string::npos)
		{
			return inpath.substr(0, pos);
		}
		else
			return inpath;	
	}

	string SplitDirFromPath(const string& inpath)
	{
		string path = NormalizePath(inpath);
		size_t pos = path.find_last_of(EXT_SEPARATOR);
		if(pos == string::npos)
		{
			// no file
			return inpath;
		}

		pos = path.find_last_of(NIX_PATH_SEPARATOR);
		if(pos!=string::npos)
		{
			return path.substr(0,pos);
		}
		return string(); // empty string
	}
	
	string CombinePathEnv(const string& dir, const string& file)
	{
		char* env = getenv(ENV_FRUND);
		if (env == NULL)
			exceptions::ThrowNoEnviromentVariable(ENV_FRUND);	
		return CombinePath(CombinePath(env, dir), file);
	}
	
	string CombinePath(const string& dir, const string& file)
	{
		string dirCopy(dir);
		string fileCopy(file);

		replace(dirCopy.begin(), dirCopy.end(), '\\', NIX_PATH_SEPARATOR); // for Linux
		replace(fileCopy.begin(), fileCopy.end(), '\\', NIX_PATH_SEPARATOR); // for Linux
		int newLength = (int)dirCopy.length() - 1;
		
		if (newLength >= 0)
			if (dirCopy[newLength] == NIX_PATH_SEPARATOR)
			{
				while (newLength > 0)
				{
					--newLength;
					if (dirCopy[newLength] != NIX_PATH_SEPARATOR) break;
				}
				if (newLength <= 0) 
					return fileCopy;
			}
		if (newLength == -1)
			return fileCopy;
		return dirCopy.substr(0, newLength + 1) + NIX_PATH_SEPARATOR + fileCopy;
	}


	/** Присоединение имени файла к каталогу с рабочей программой
	* @param fileName - имя файла
	* @param argvPath - путь из argv[0], используется только для NIX
	* @return путь из катлога с программой и имени передаваемого файла
	*/
	string GetFullPath
		(
			const string& fileName,
			const string&
		#ifdef NIX
				argvPath
		#endif
		)
	{
		string fullPath, runDirectory, newRunDirectory;
#ifdef NIX
		newRunDirectory = argvPath;
#else

#	ifdef UNICODE
		//// TODO: исправить использование неинициализированной переменной
		//char szPath[MAX_PATH * sizeof(char)];
		//string s(szPath);
		//int len;
		//int slength = (int)s.length() + 1;
		//len = MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, 0, 0); 
		//wchar_t* buf = new wchar_t[len + 1];
		//MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, buf, len);
		//GetModuleFileName((HMODULE) GetModuleHandle(0), buf, MAX_PATH * sizeof(char));
		//newRunDirectory = szPath;
		wchar_t executablePath[MAX_PATH];
		HMODULE executableHandle = GetModuleHandle(NULL);

		if (executableHandle != NULL)
		{
			DWORD length = GetModuleFileNameW(executableHandle, executablePath, sizeof(executablePath)); 

			if
				(
					!
					(
						(length == sizeof(executablePath)) &&
						(GetLastError() == ERROR_INSUFFICIENT_BUFFER)
					)
				)
			{
				newRunDirectory = executablePath;
			}
		}

#	else
		char executablePath[MAX_PATH];
		HMODULE executableHandle = GetModuleHandle(NULL);

		if (executableHandle != NULL)
		{
			DWORD length = GetModuleFileName(executableHandle, executablePath, sizeof(executablePath)); 

			if
				(
					!
					(
						(length == sizeof(executablePath)) &&
						(GetLastError() == ERROR_INSUFFICIENT_BUFFER)
					)
				)
			{
				newRunDirectory = executablePath;
			}
		}
#	endif
#endif
		replace(newRunDirectory.begin(), newRunDirectory.end(), '\\', NIX_PATH_SEPARATOR);// for Windows

		size_t directorySeparatorPos;

		if ((directorySeparatorPos = newRunDirectory.find_last_of(NIX_PATH_SEPARATOR)) != string::npos)
		{
			newRunDirectory.resize(directorySeparatorPos);
			runDirectory = newRunDirectory;
		}
		fullPath = runDirectory;
		if (!fileName.empty())
		{
			fullPath += NIX_PATH_SEPARATOR + fileName; 
		}
		
		return fullPath;
	}

	bool ReadLine(ifstream& ifs, stringstream& line)
	{
		char buf[STRLIMIT];

		ifs.getline(buf, STRLIMIT);

		size_t count = strlen(buf);

		if (count > 0)
			if (buf[count - 1] == '\r')
				buf[count - 1] = 0;
		if (ifs.eof())
			exceptions::ThrowFileInvalidFormat();
		line.clear();
		line << buf;

		return buf[0] != 0;
	}

	bool ReadLineNoException(ifstream& ifs, stringstream& line)
	{
		char buf[STRLIMIT];

		ifs.getline(buf, STRLIMIT);

		size_t count = strlen(buf);

		if (count > 0)
			if (buf[count - 1] == '\r')
				buf[count - 1] = 0;
		if (ifs.eof())
			return false;
		line.clear();
		line << buf;

		return buf[0] != 0;
	}

	void ReadLineString(ifstream& ifs, char* line)
	{
		ifs.getline(line,STRLIMIT);

		size_t count = strlen(line);

		if (count > 0)
			if (line[count - 1]=='\r')
				line[count - 1]= 0;
		if (ifs.eof())
			exceptions::ThrowFileInvalidFormat();
	}

	void ReadLineString(ifstream& ifs, string& line)
	{
		char buf[STRLIMIT];

		ifs.getline(buf,STRLIMIT);

		size_t count = strlen(buf);

		if (count > 0)
			if (buf[count - 1] == '\r')
				buf[count - 1] = 0;
		if (ifs.eof())
			exceptions::ThrowFileInvalidFormat();
		line = buf;	
	}

	bool ReadLineStringNoException(ifstream& ifs, string& line)
	{
		char buf[STRLIMIT];

		ifs.getline(buf,STRLIMIT);

		size_t count = strlen(buf);

		if (count > 0)
			if (buf[count - 1] == '\r')
				buf[count - 1] = 0;
		if (ifs.eof())
			return false;
		line = buf;	

		return true;
	}

	bool ReadLineStringNoException2(ifstream& ifs, string& line)
	{
		char buf[STRLIMIT];

		ifs.getline(buf,STRLIMIT);

		size_t count = strlen(buf);

		if (count > 0)
			if (buf[count - 1] == '\r')
				buf[count - 1] = 0;
		if ((ifs.eof()) && (count == 0))
			return false;
		line = buf;

		return true;
	}


	/** Проверка существования файла
	* @param path - путь у файлу
	* @return 0, если не существует
	*/
	int FileExist(const char* path)
	{
		int res = access(path, 0) == -1 ? 0 : 1;
		return res;
	}

	/** Проверка расширения файла
	* @param path - путь у файлу
	* @param ext - проверяемое расширение
	* @return 0, если не существует
	*/
	int CheckExtension(const char* path, const char* ext)
	{
		const char* extp;
		if((extp=strstr(path,ext)) != NULL)
		{
			if(extp-path == strlen(path)-strlen(ext))
				return 1;
		}
		return 0;
	}

	void DeleteFiles(const char* path)
	{
		fs::DirectoryLister* fe = new fs::DirectoryLister(".", string(path));
		string filename;
		fe->Init();
		while(!fe->IsLast())
		{
			fe->Next(filename);
			unlink(filename.c_str());
		}
		delete fe;
	}
	
	void DeleteFilesWithPostfix(const string& mask, const string& prefix, char postfixSeparator, int postfixLength)
	{		
		fs::DirectoryLister* fe = new fs::DirectoryLister(".",string(mask));
		string filename;
		fe->Init();
		while(!fe->IsLast())
		{
			fe->Next(filename);
			string fileNameNoExt = fs::SplitFileFromExt(filename);
			int dl = (int)fileNameNoExt.length() - (int)prefix.length();
			//std::cout << "CHECK: " << fileNameNoExt << " to " << prefix << std::endl;

			if(fileNameNoExt.substr(0, prefix.length()) == prefix && dl > 0)
			{
				if(fileNameNoExt[prefix.length()] == postfixSeparator && dl < postfixLength+2)
				{
					//std::cout << "UNLINK: " << filename << std::endl;
					unlink(filename.c_str());
				}
			}
		}
		delete fe;
	}
	void CopySingleFile(const char* srcPath, const char* destPath)
	{
		ifstream ifs(srcPath, ios::binary);
		if(ifs.is_open())
		{
			ofstream ofs(destPath,ios::binary);
			ofs << ifs.rdbuf();			
			ofs.close();
			ifs.close();
		}
		else
		{
			exceptions::ThrowFileNotFound(srcPath);
		}
	}
	
	void RenameFile(const char* srcPath, const char* destPath)
	{
		ifstream ifs(srcPath, ios::binary);
		if(ifs.is_open())
		{
			ifs.close();
			unlink(destPath);
			rename(srcPath, destPath);
		}
		else
		{
			exceptions::ThrowFileNotFound(srcPath);
		}
	}
	
	void MergeFiles(const char* path1, const char* path2, const char* destPath)
	{
		ifstream ifs1(path1);
		if(!ifs1.is_open())
		{
			exceptions::ThrowFileNotFound(path1);
		}
		ifstream ifs2(path2);
		if(!ifs2.is_open())
		{
			ifs1.close();
			exceptions::ThrowFileNotFound(path2);
		}
		
		ofstream ofs(destPath);
		ofs << ifs1.rdbuf();
		ofs << ifs2.rdbuf();
		ofs.close();
		ifs1.close();
		ifs2.close();
	}

	void AppendFile(const char* path1, const char* path2)
	{
		ifstream ifs1(path1);
		if(!ifs1.is_open())
		{
			exceptions::ThrowFileNotFound(path1);
		}
		ifs1.close();
		ifstream ifs2(path2);
		if(!ifs2.is_open())
		{
			exceptions::ThrowFileNotFound(path2);
		}

		ofstream ofs(path1, ios::app);
		ofs.seekp(0, ios::end);
		ofs << ifs2.rdbuf();			
		ofs.close();
		ifs2.close();
	}

	int Exec ( const string& process, const string& cmd)
	{
#ifdef NIX
		char buf[1000];

		cout << "cur_dir: " << getcwd(buf, 1000) << endl;
		cout << "env_path: " << getenv(ENV_PATH) << endl;
		cout << "env_lib: " << getenv(ENV_LIB) << endl;
		cout << "env_include: " << getenv(ENV_INCLUDE) << endl; 
		cout << "process: " << process << endl;
		cout << "cmd: " << cmd << endl;

		pid_t pid = fork();

		if (pid == 0)
		{
			string cmdCopy(cmd);

			cmdCopy = pathStringReplace(cmdCopy, "*.f");
			cmdCopy = pathStringReplace(cmdCopy, "*.o");
			cmdCopy = pathStringReplace(cmdCopy, "*.cpp");

			vector<string> args = tokenize(cmdCopy, " ");
			char** argv = new char*[args.size()+1];
			int i;

			for (i = 0; i < args.size(); i++)
			{
				argv[i] = strdup(args[i].c_str());
				cout << "argv " << argv[i] << endl;
			}
			argv[i] = NULL;

			char *env[] = {"HOME=/usr/home", "LOGNAME=home", "PATH=/usr/bin", NULL };
			int result = execve(process.c_str(), argv, env );

			if (result != 0)
			{
				perror(process.c_str());
			}
			for (i = 0; i < args.size(); i++)
			{
				delete [] argv[i];
			}
			delete [] argv;
			cout << "Execve failed\n";

			return -1;
		}
		else
		{ // parent process
			int res = 0;
			int status;

			do
			{
				res = waitpid(pid, &status, WUNTRACED | WCONTINUED);
				if (res == -1)
				{
					cout << "Waitpid failed\n";
					return -2;
				}
				if (WIFEXITED(status))
				{
					cout << "Exited\n";
					res = WEXITSTATUS(status);
				}
				else if (WIFSIGNALED(status))
				{
					cout << "Killed by signald\n";
					res = WTERMSIG(status);
				}
				else if (WIFSTOPPED(status))
				{
					cout << "Stopped by signal\n";
					//return WSTOPSIG(status);
				} 
				else if (WIFCONTINUED(status))
				{
					cout << "Continued\n";
					//return 0;
				}
			}
			while (!WIFEXITED(status) && !WIFSIGNALED(status));

			return res;
		}
#else
		string cmdCopy(cmd);

#	ifdef GFORTRAN
		cmdCopy = pathStringReplace(cmdCopy, "*.f");
		cmdCopy = pathStringReplace(cmdCopy, "*.o");
		cmdCopy = pathStringReplace(cmdCopy, "*.cpp");
#	endif

		STARTUPINFO si;
		PROCESS_INFORMATION pi;
		ZeroMemory(&si, sizeof(si));
		si.cb = sizeof(si);
		ZeroMemory(&pi, sizeof(pi));
		string sCmd = process;
		sCmd = sCmd + " " + string(cmdCopy);
		char* ccmd = strdup(sCmd.c_str());

#	ifdef UNICODE
		string s(ccmd);
		int len;
		int slength = (int)s.length() + 1;
		len = MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, 0, 0); 
		wchar_t* buf = new wchar_t[len];

		MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, buf, len);
		if
			(
				!CreateProcess
					(
						NULL,	// No module name (use command line)
						buf,	// Command line
						NULL,	// Process handle not inheritable
						NULL,	// Thread handle not inheritable
						FALSE,	// Set handle inheritance to FALSE
						0,		// No creation flags
						NULL,	// Use parent's environment block
						NULL,	// Use parent's starting directory 
						&si,	// Pointer to STARTUPINFO structure
						&pi		// Pointer to PROCESS_INFORMATION structure
					)
			)
#	else
		// Start the child process. 
		if
			(
				!CreateProcess
					(
						NULL,	// No module name (use command line)
						ccmd,	// Command line
						NULL,	// Process handle not inheritable
						NULL,	// Thread handle not inheritable
						FALSE,	// Set handle inheritance to FALSE
						0,		// No creation flags
						NULL,	// Use parent's environment block
						NULL,	// Use parent's starting directory 
						&si,	// Pointer to STARTUPINFO structure
						&pi		// Pointer to PROCESS_INFORMATION structure
					)
		) 
#	endif
		{
			printf( "CreateProcess failed (%lu)\n", GetLastError() );
			free(ccmd);

			return -1;
		}

		DWORD dwExitCode;

		// Wait until child process exits.
		WaitForSingleObject(pi.hProcess, INFINITE);
		if (!GetExitCodeProcess(pi.hProcess, &dwExitCode))
		{
			dwExitCode = -1;
		}
	
		// Close process and thread handles. 
		CloseHandle(pi.hProcess);
		CloseHandle(pi.hThread);
		free(ccmd);

		return dwExitCode;
#endif
	}

#ifdef NIX

    vector<string> GetListFile(const char* ext)
    {
        vector<string> files;
        glob_t globbuf;
        globbuf.gl_offs = 0;
        glob(ext, GLOB_DOOFFS, NULL, &globbuf);
        for (int i = 0; i < (int) (globbuf.gl_pathc); i++)
            files.push_back(globbuf.gl_pathv[i]);
        return files;
    }
    
    void FormatCMD(string &keys, const char* ext)
    {
        list<string> args;
        vector<string> args_copy = tokenize(keys, " ");
        args.insert(args.begin(), args_copy.begin(), args_copy.end());

        list<string>::iterator filesType;
        filesType = find(args.begin(), args.end(), ext);
        if (filesType != args.end())
        {
            vector<string> listFiles =fs::GetListFile(ext);
            args.insert(filesType, listFiles.begin(), listFiles.end());
            
            filesType = find(args.begin(), args.end(), ext);
            args.erase(filesType);
            keys.clear();
            for(list<string>::const_iterator i = args.begin(); i != args.end(); ++i)
                keys += *i + " ";
        }
    }

#endif
// FROM TESTER

	string GetRealPath(const string& inpath)
	{
		string path = inpath;

		if(path.find_first_of(':') != string::npos) // absolute
			return path;

		replace(path.begin(),path.end(),'\\',NIX_PATH_SEPARATOR); // for Linux
		string currentDir = GetCurrentDir();
		while(path.find_first_of("../") == 0)
		{
			path = path.substr(3,path.size()-3);
			size_t pos = currentDir.find_last_of(NIX_PATH_SEPARATOR);
			if(pos != string::npos)
				currentDir = currentDir.substr(0, pos);
		}
		return currentDir + DirSeparator() + path;
	}

	string DirSeparator()
	{
		string separator;
		separator += NIX_PATH_SEPARATOR;
		return separator;
	}

	string GetCurrentDir()
	{
#ifdef NIX
		char buf[PATH_MAX];
		string currentDir = getcwd(buf, PATH_MAX);
		//replace(currentDir.begin(),currentDir.end(),'\\',NIX_PATH_SEPARATOR); // for Linux
		//if(currentDir[currentDir.size()-1] == NIX_PATH_SEPARATOR)
		//	currentDir = currentDir.substr(0,currentDir.size()-1);
		return currentDir;
#else
		char buf[MAX_PATH];           
		string currentDir = _getcwd(buf, MAX_PATH);	
		replace(currentDir.begin(),currentDir.end(),'\\',NIX_PATH_SEPARATOR); // for Linux
		if(currentDir[currentDir.size()-1] == NIX_PATH_SEPARATOR)
			currentDir = currentDir.substr(0,currentDir.size()-1);

		return currentDir;
#endif
	}

	bool ExistDir(const string& path)
	{
		string dir = GetCurrentDir();
#if defined(NIX) || defined (GNUCPP)
		int res = chdir(path.c_str());
#else
		int res = _chdir(path.c_str());
#endif

		SetCurrentDir(dir);
		return res != -1;
	}

	void SetCurrentDir(const string& path)
	{
#ifdef GNUCPP
		if(chdir(path.c_str()) == -1L)
			exceptions::ThrowFolderNotFound(path);
#else
		if(_chdir(path.c_str()) == -1L)
			exceptions::ThrowFolderNotFound(path);
#endif
	}

	void SetCurrentDir(const string& path, string& oldPath)
	{
		oldPath = GetCurrentDir();
		SetCurrentDir(path);
	}

	string GetParentDir(const string& inpath)
	{
		string path = inpath;
		replace(path.begin(),path.end(),'\\',NIX_PATH_SEPARATOR); // for Linux
		size_t pos = path.find_last_of(NIX_PATH_SEPARATOR);
		if(pos != string::npos)
			return path.substr(0, pos);
		else
			return "";
	}

	void MakeDir(const string& path)
	{
		stack<string> pathsToCreate;
		pathsToCreate.push(path);
		string currentParent = GetAbsolutePath(path);

		while(!ExistDir(pathsToCreate.top()))
		{
			currentParent = GetParentDir(currentParent);
			if(currentParent.size() != 0)
				pathsToCreate.push(currentParent);
			else
				break;
		}
		pathsToCreate.pop(); // exclude drive or paret exiting dir
		while(!pathsToCreate.empty())
		{
			currentParent=pathsToCreate.top();
#ifdef NIX
			if(mkdir(currentParent.c_str(), 0777) == -1L)
				exceptions::ThrowFolderNotFound(path);
#elif defined (GNUCPP)
			if(mkdir(currentParent.c_str()) == -1L)
				exceptions::ThrowFolderNotFound(path);
#else
			if(_mkdir(currentParent.c_str()) == -1L)
				exceptions::ThrowFolderNotFound(path);
#endif
			pathsToCreate.pop();
		}
	}

	void DeleteFile(const string& path)
	{
#ifdef GNUCPP
		unlink(path.c_str());
#else
		_unlink(path.c_str());
#endif
	}

	bool FileExist(const string& path)
	{
#ifdef GNUCPP
		return !(access(path.c_str(), 0) == -1L);
#else
		return !(_access(path.c_str(), 0) == -1L);
#endif

	}

	bool CheckFileFlag(const string& path)
	{
		if (FileExist(path))
		{
			ifstream ifs(path);
			if (ifs.is_open())
			{
				int value = 0;
				stringstream line;
				bool res = ReadLineNoException(ifs, line);
				if (res)
					line >> value;
				ifs.close();
				return value == 0 ? false : true;
			}
		}
		return false;
	}

	/** Копирует файл в каталог
	* @param file - имя файла
	* @param folder - каталог
	*/
	void CopyFileToFolder(const string& file,const string& folder)
	{
		string pureFile = SplitFileFromPath(file);
		if (pureFile.empty()) 
			return;
		string destPath = GetAbsolutePath(folder) + NIX_PATH_SEPARATOR + pureFile;
		ifstream ifs(file, ios::binary);

		if (ifs.is_open())
		{
			ofstream ofs(destPath, ios::binary);

			if (!ofs.is_open())
			{
				exceptions::ThrowFileNotOpened(destPath);
			}
			ofs << ifs.rdbuf();
			ofs.close();
			ifs.close();
		}
		else
		{
			exceptions::ThrowFileNotFound(file);
		}
	}

	void CopyFilesToFolder(const vector<string>& files, const string& folder)
	{
		vector<string> filesToCopy;
		vector<string>::const_iterator it = files.begin();
		string filename;
		while(it != files.end())
		{
			fs::DirectoryLister dl(".",*it);
			dl.Init();
			while(!dl.IsLast())
			{
				dl.Next(filename);
				filesToCopy.push_back(filename);
			}
			++it;
		}
		
		it = filesToCopy.begin();
		while(it != filesToCopy.end())
		{
			CopyFileToFolder(*it,folder);
			++it;
		}
	}


	/** Удаляет файлы из каталога
	* @param folder - каталог
	*/
	void CleanFolder(const string& _folder)
	{
		string currentDir = fs::GetCurrentDir();	
		SetCurrentDir(_folder);
		// удаление файлов
		fs::DirectoryLister dlToClean;
		for(size_t i = 0; i < dlToClean.GetFilesCount();i++)
		{
			fs::DeleteFile(dlToClean.GetFile(i));
		}
		SetCurrentDir(currentDir);	
	}


	/** Копирует файлы из каталога в каталог
	* @param inFolder - исходный каталог
	* @param outFolder - каталог, в который копировать
	* если его нет, он создается
	*/
	void CopyFolder(const string& _inFolder, const string& _outFolder)
	{
		fs::DirectoryLister dl(".");
		string testDir = fs::GetCurrentDir();	

		if(!fs::ExistDir(_outFolder))
			fs::MakeDir(_outFolder);
		else
			fs::CleanFolder(_outFolder);

		fs::DirectoryLister dlToCopy(_inFolder);
		SetCurrentDir(_inFolder);
		for(size_t i = 0; i < dlToCopy.GetFilesCount();i++)
		{
			fs::CopyFileToFolder(dlToCopy.GetFile(i), _outFolder);
		}
		fs::SetCurrentDir(testDir);

	}

	string AppendToFileName(const string& fileName, const string& str)
	{
		string tmpStr = fileName;
		size_t pos = fileName.find_last_of(".");
		string res = (pos == string::npos)? tmpStr.append(str): tmpStr.insert(pos,str);
		return res;
	}

	bool ReadLineKeyValue(ifstream& ifs, string& key, string& value)
	{
		stringstream line;
		bool res = ReadLineNoException(ifs, line);
		if(res)
			line >> key >> value;
		return res;
	}
}

#ifdef WINDOWS_FUNCTION_NAME_DEFINED

#define DeleteFile DELETE_FILE_DEFINE

#undef DELETE_FILE_DEFINE

#undef WINDOWS_FUNCTION_NAME_DEFINED

#endif
