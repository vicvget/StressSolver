#include "../FrundFacade/FElementType.h"
#include "../../../fcore/wrappers/StringRoutines.h"

namespace FElementType
{
	string MakeFileId(FElementTypeCodes code, int number)
	{
		string fileId = "00000000";
		string codeStringPart = NumberToString(code);
		string numberStringPart = NumberToString(number);
		int nZeroes = fileId.length()-codeStringPart.length()-numberStringPart.length();
		if(nZeroes < 0)
		{
			// TODO: exception
		}
		fileId = codeStringPart.append(fileId, codeStringPart.length(), nZeroes).append(numberStringPart);
		return fileId;
	}
}