#ifndef FELEMENT_TYPE_H

#define FELEMENT_TYPE_H


#include <string>


using std::string;


namespace FElementType
{
	enum FElementTypeCodes
	{
		FETC_Body=1,
		FETC_JointChar=2,
		FETC_Header=3,
		FETC_CUnit=4,
		FETC_Joint=5,
		FETC_ForceChar=6,
		FETC_Force=7
	};

	string MakeFileId(FElementTypeCodes code, int number);

}


#endif // FELEMENT_TYPE_H