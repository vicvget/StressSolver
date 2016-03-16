#ifndef PARAMS_SET_BASE_H

#define PARAMS_SET_BASE_H


#include <string>
#include <fstream>


using std::string;
using std::ifstream;
using std::ofstream;


#define DEFAULT_NAME "default"

/** Базовый класс именованного 
* набора параметров (для междисциплинарного расчета)
*
* @author Getmanskiy Victor
*/
class ParamsSetBase
{
protected:
	// имя набора переменных
	string _name;

public:

	ParamsSetBase();
	ParamsSetBase(const string& name);
	ParamsSetBase(const ParamsSetBase& copy):_name(copy._name){}
	virtual ~ParamsSetBase();
	
	virtual void Save(ofstream& ofs) const = 0;
	virtual void Load(ifstream& ifs) = 0;
	virtual void Init() {};

	void DefaultLoad(ifstream& ifs);
	void DefaultSave(ofstream& ofs) const;

#pragma region properties
	string GetName() const {return _name;}
	void SetName(string name) {_name = name;}
#pragma endregion

};


#endif // PARAMS_SET_BASE_H