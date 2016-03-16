#include "ResultsFactory.h"

#include "ProviderRezr.h"
#include "ResultsReader.h"

#include <fstream>


/**
* Создать провайдер доступа к данным результатов расчета MBS-решателя
* @param resultsFileName - наименование файла с результатами расчета MBS-решателя
* @return созданный провайдер доступа к данным результатов расчета MBS-решателя
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
* Создать умный указатель на провайдер доступа к данным результатов расчета MBS-решателя
* @param resultsFileName - наименование файла с результатами расчета MBS-решателя
* @return умный указатель на созданный провайдер доступа к данным результатов расчета MBS-решателя
*/
AbstractProviderRezrPointer ResultsFactory::CreateProviderPointer
	(
		const string& resultsFileName
	)
{
	return AbstractProviderRezrPointer{CreateProvider(resultsFileName)};
}