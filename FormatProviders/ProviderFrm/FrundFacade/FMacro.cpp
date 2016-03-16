#include "FMacro.h"

#include "../../../fcore/fcore.h"

#include <cstring>


using std::istringstream;
using std::endl;


FMacro::FMacro
	(
		ifstream& ifs
	)
	:
		x(),
		y()
{	
	fs::ReadLineString(ifs, name);
	fs::ReadLineString(ifs, type);
	string open, param, close;

	fs::ReadLineString(ifs, open);
	fs::ReadLineString(ifs, open);
	fs::ReadLineString(ifs, param);
	fs::ReadLineString(ifs, close);

	while (close.compare("namer"))
	{
		mParams.push_back(param);
		fs::ReadLineString(ifs,open);
		fs::ReadLineString(ifs,param);
		fs::ReadLineString(ifs, close);
	}

	fs::ReadLineString(ifs,open);
	fs::ReadLineString(ifs,param);
	fs::ReadLineString(ifs, close);
	while (close.compare("position"))
	{
		rParams.push_back(param);
		fs::ReadLineString(ifs,open);
		fs::ReadLineString(ifs,param);
		fs::ReadLineString(ifs,close);
	}
	istringstream strx(open);
	istringstream stry(param);
	strx >> x;
	stry >> y;
}

/** Сохранение макроса
* @param ofs - выходной поток
*/
void FMacro::Write(ofstream& ofs) const
{
	ofs << name << endl << 
		type << endl <<
		"namem" << endl;
	for(size_t i = 0; i < mParams.size(); i++)
	{
		ofs << 'm' << endl << 
			mParams[i] << endl <<
			'm' << endl;
	}
	ofs << name << endl << 
		type << endl <<
		"namer" << endl;
	for(size_t i = 0; i < rParams.size(); i++)
	{
		// TODO: воспроизведение бага для соответствия тарому формату
		//ofs << 'r' << endl << 
		ofs << 'r' << endl << 
			rParams[i] << endl <<
			' ' << endl;
	}
	ofs << x << endl << y << endl << "position" << endl;
}

/** Вывод в текстовый файловый поток
* @param ofs - выходной поток
* @param isRType - tree на шаге Solve
*/
void FMacro::Output(ofstream& ofs, int isRType) const
{
	ofs << "@ " << name << endl << 
		'#' << type << ' ' ;
	if(isRType)
	{
		for(size_t i = 0; i < rParams.size(); i++)
			ofs << rParams[i] << ' ';
		ofs << endl;
	}
	else
	{
		for(size_t i = 0; i < mParams.size(); i++)
			ofs << mParams[i] << ' ';
		ofs << endl;
	}
}

/** Получить количество строк, занимаемых в заголовке (все 3 секции mParams, rParams, Position)
* @return количество строк для макросов
*/
int FMacro::StrCount() const
{
	return 3 + mParams.size() + rParams.size();
}

