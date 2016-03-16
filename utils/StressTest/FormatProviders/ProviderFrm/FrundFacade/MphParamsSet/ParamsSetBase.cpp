#include "ParamsSetBase.h"
#include "../../../../fcore/wrappers/FileRoutines.h"

ParamsSetBase::ParamsSetBase(): _name(DEFAULT_NAME)
{
	Init();
}

ParamsSetBase::ParamsSetBase(const string& name): _name(name)
{
	Init();
}


ParamsSetBase::~ParamsSetBase() 
{

}

void ParamsSetBase::DefaultLoad(ifstream& ifs)
{
	fs::ReadLineString(ifs, _name);
}

void ParamsSetBase::DefaultSave(ofstream& ofs) const
{
	ofs << _name << std::endl;
}
