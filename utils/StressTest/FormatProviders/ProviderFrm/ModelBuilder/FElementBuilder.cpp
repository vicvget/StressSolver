#include "FElementBuilder.h"

void FElementBuilder::SetupCommonProperties(FElement& element, const string& name, int id, const string& fileId)
{
	element.Name(name);
	element.Id(id);
	element.Number(id+1);
	element.FileId(fileId);
}

void FElementBuilder::SetupCommonBodyProperties(FElement& element, const string& name, int id, int number, const string& fileId)
{
	element.Name(name);
	element.Id(id);
	element.Number(number);
	element.FileId(fileId);
}
