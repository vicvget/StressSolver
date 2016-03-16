#ifndef FMACRO_H

#define FMACRO_H


#include "BasicTypes.h"

//TODO: нормализовать комментарии

/**
* Класс для считывания макроса модели ФРУНД
*
* @author Getmanskiy Victor
*/
class FMacro
{
private:
	//char name[STRLIMIT];
	//char type[STRLIMIT];

	string name;
	string type;

	//TODO: загрузка коментария 
	vector<string> mParams;
	vector<string> rParams;

	// экранные координаты на диаграмме
	int x, y;

public:

	FMacro()
		:
			x(),
			y()
	{
	}

	/** Считывание макроса
	* @param ifs - входной поток
	*/
	FMacro
		(
			ifstream& ifs
		);

	/** Вывод макроса во входной файл описания модели (mdl, ras)
	* @param ofs - выходной поток
	* @param isRType - true для solve step
	*/
	void Output(ofstream& ofs, int isRType) const;

	/** Сохранение макроса
	* @param ofs - выходной поток
	*/
	void Write(ofstream& ofs) const;

	/** Получить количество строк, занимаемых в заголовке (все 3 секции mParams, rParams, Position)
	* @return кличество строк
	*/
	int StrCount() const;

	int GetX() const {return x;}
	void SetX(int val){x = val;}

	int GetY() const {return y;}
	void SetY(int val){y = val;}



	string Type() const { return type; }
	void Type(string val) { type = val; }


	string Name() const { return name; }
	void Name(string val) { name = val; }

	//char* GetName(){return name;}
	//char* GetType(){return type;}

	vector<string>& GetRParams()
	{
		return rParams;
	}

	const vector<string>& GetRParams() const
	{
		return rParams;
	}

	vector<string>& GetMParams()
	{
		return mParams;
	}

	const vector<string>& GetMParams() const
	{
		return mParams;
	}

};


#endif // FMACRO_H