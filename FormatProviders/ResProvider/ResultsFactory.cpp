#include "ResultsFactory.h"

#include "ProviderRezr.h"
#include "ResultsReader.h"

#include <fstream>


/**
* ������� ��������� ������� � ������ ����������� ������� MBS-��������
* @param resultsFileName - ������������ ����� � ������������ ������� MBS-��������
* @return ��������� ��������� ������� � ������ ����������� ������� MBS-��������
*/
AbstractProviderRezr* ResultsFactory::CreateProvider
	(
		const string& resultsFileName
	)
{
	AbstractProviderRezr* resultsProvider{};
	ifstream resultsFile(resultsFileName, ifstream::binary);

	if (resultsFile.is_open())
	{
		char symbol;

		resultsFile.read(&symbol, sizeof(char));
		switch (symbol)
		{
		case 'K':
			resultsProvider = new ProviderRezr(resultsFileName, "fadres.dat");
			break;
			
		case 'V':
			resultsProvider = new ResultsReader(resultsFileName, "fadres.dat");
			break;
		}
		resultsFile.close();
	}

	return resultsProvider;
}

/**
* ������� ����� ��������� �� ��������� ������� � ������ ����������� ������� MBS-��������
* @param resultsFileName - ������������ ����� � ������������ ������� MBS-��������
* @return ����� ��������� �� ��������� ��������� ������� � ������ ����������� ������� MBS-��������
*/
AbstractProviderRezrPointer ResultsFactory::CreateProviderPointer
	(
		const string& resultsFileName
	)
{
	return AbstractProviderRezrPointer{CreateProvider(resultsFileName)};
}