#ifndef FELEMENT_BUILDER_H

#define FELEMENT_BUILDER_H


#include "../../ProviderFrm/FrundFacade/FElement.h"

#include <string>


using std::string;


class FElementBuilder
{
public:
	static void SetupCommonProperties(FElement& element, const string& name, int id, const string& fileId);
	static void SetupCommonBodyProperties(FElement& element, const string& name, int id, int number, const string& fileId);
};


#endif // FELEMENT_BUILDER_H