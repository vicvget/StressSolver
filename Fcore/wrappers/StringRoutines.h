#ifndef STRING_ROUTINES_H

#define STRING_ROUTINES_H


#include <string>
#include <vector>
#include <cstdio>
#include <cctype>
#include <sstream>


using std::string;
using std::vector;
using std::stringstream;


#pragma region ServiceFunctions

/** разбивает строку на подстроки по заданным разделителям
* http://www.digitalpeer.com/id/simple
* @param str - строка
* @param delimeters - разделители
* @return массив подстрок
*/
vector<string> tokenize(const string& str,const string& delimiters, bool isEmptyTokens=true);

/** считывание числа из строки
* @param str - строка
* @return число
*/
template <class T>
T StringToNumber(string str)
{
	stringstream cmdstream(str);
	T val;
	cmdstream >> val;
	return val;
}

/** считывание числа из строки
* @param str - строка
* @return число
*/
template <class T>
string NumberToString(T val)
{
	stringstream stream;
	stream << val;
	return stream.str();
}

string Translit(const string& str, bool isAlt =  false);

void trim1(string& str);

/** Замена подстроки типа файла нас список файлов в строке 
 */
string pathStringReplace(string source, string ext);

#pragma endregion

#ifdef GNUCPP
#ifdef NIX
/** Перевод строки в нижний регистр
* @param str - исходная строка, модифицируется
* @return переведенная строка
*/
char* strlwr(char* str);

/** Перевод строки в верхний регистр
* @param str - исходная строка, модифицируется
* @return переведенная строка
*/
char* strupr(char* str);

/** Сравнение строк без учета регистра
* @param str1 - строка 1
* @param str2 - строка 2
* @param len - количество сравниваемых символов
* @return 0 = , -1 < , 1 >
*/
int strnicmp(const char* str1, const char* str2, unsigned int len);

#endif
#endif


#endif // STRING_ROUTINES_H