#include "StringRoutines.h"

#include "fcUnistd.h"
#include "DirectoryLister.h"
#include "../fcore.h"
#include "Translit.h"
#include "cstdlib"


using std::stringstream;


#pragma region ServiceFunctions


string translit_table[256];
string translit_table_alt[256];

void InitTranslitTable()
{
	for(int i = 0; i < 66; i++)
	{
		size_t id = (unsigned char)rus_chars[i];
		translit_table[id] = trs_chars[i];		
		translit_table_alt[id] = trs_chars_alt[i];
	}	
}

string Translit(const string& str, bool isAlt)
{
	static bool isTableInitialized = false;
	if(!isTableInitialized)
	{
		InitTranslitTable();
		isTableInitialized = true;
	}

	string res;
	stringstream stream;

	for (size_t i = 0; i < str.size(); ++i)
	{
		if(isascii(str[i]))
		{
			stream << str[i];
		}
		else 
		{			
			stream << (isAlt ? translit_table_alt[(unsigned char)str[i]]: translit_table[(unsigned char)str[i]]);
		}
	}
	char buf[STRLIMIT];
	stream.getline(buf,STRLIMIT);
	res = string(buf);
	return res;
}

/** разбивает строку на подстроки по заданным разделителям
* http://www.digitalpeer.com/id/simple
* @param str - строка
* @param delimeters - разделители
* @return массив подстрок
*/
vector<string> tokenize(const string& str,const string& delimiters, bool isEmptyTokens)
{
	vector<string> tokens;

	string::size_type lastPos = 0, pos = 0;  
	size_t count = 0;

	if (str.length() < 1) return tokens;

	// skip delimiters at beginning.  
	lastPos = str.find_first_not_of(delimiters, 0);
	if (isEmptyTokens)
	{
		if ((str.substr(0, lastPos - pos).length()) > 0)
		{
			count = str.substr(0, lastPos - pos).length();

			for (size_t i = 0; i < count; ++i)
				tokens.push_back("");
			if (string::npos == lastPos)
				tokens.push_back("");
		}
	}

	// find first "non-delimiter".
	pos = str.find_first_of(delimiters, lastPos);
	while (string::npos != pos || string::npos != lastPos)
	{
		// found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));

		// skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
		
		if (isEmptyTokens)
		{
			if ((string::npos != pos) && (str.substr(pos, lastPos-pos).length() > 1))
			{
				count = str.substr(pos, lastPos-pos).length();

				for (size_t i = 0; i < count; ++i)
					tokens.push_back("");
			}
		}
		pos = str.find_first_of(delimiters, lastPos);
	}

	return tokens;
}

string pathStringReplace(string source, string ext) 
{
    string filesList;
	fs::DirectoryLister fe(".",ext);
    fe.Init();
    string curPath;
    while (!fe.IsLast()) 
    {
        fe.Next(curPath);
        filesList += curPath + " ";
    }
    trim1(filesList);
	auto pos = source.find(ext);
    if(pos != string::npos)
	{
        source.replace(pos, ext.length(), filesList);
	}
    return source;
}

void trim1(string& str)
{
  string::size_type pos1 = str.find_first_not_of(' ');
  string::size_type pos2 = str.find_last_not_of(' ');
  str = str.substr(
	pos1 == string::npos ? 0 : pos1, 
    pos2 == string::npos ? str.length() - 1 : pos2 - pos1 + 1);
}
#pragma endregion

#ifdef GNUCPP

char* strlwr(char* str)
{
    for(int i = 0; i < strlen(str); i++)
    {
        str[i] = (char)(tolower(str[i]));
    }
    return str;
}

char* strupr(char* str)
{
    for(int i = 0; i < strlen(str); i++)
    {
        str[i] = (char)(toupper(str[i]));
    }
    return str;
}


int strnicmp(const char* str1, const char* str2, unsigned int len)
{
    char* ds1 = strdup(str1);
    char* ds2 = strdup(str2);
    
    strlwr(ds1);
    strlwr(ds2);

    int res = strncmp(ds1,ds2,len);
    free(ds1);
    free(ds2);
    return res;
}


#endif
