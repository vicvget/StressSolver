#ifndef RESULTS_FACTORY_H

#define RESULTS_FACTORY_H


#include "AbstractProviderRezr.h"


/*
* Пространство имен для функций создания провайдера доступа
* к данным результатов расчета MBS-решателя
*/
namespace ResultsFactory
{

	/**
	* Создать провайдер доступа к данным результатов расчета MBS-решателя
	* @param resultsFileName - наименование файла с результатами расчета MBS-решателя
	* @return созданный провайдер доступа к данным результатов расчета MBS-решателя
	*/
	AbstractProviderRezr* CreateProvider
		(
			const string& resultsFileName
		);

	/**
	* Создать умный указатель на провайдер доступа к данным результатов расчета MBS-решателя
	* @param resultsFileName - наименование файла с результатами расчета MBS-решателя
	* @return умный указатель на созданный провайдер доступа к данным результатов расчета MBS-решателя
	*/
	AbstractProviderRezrPointer CreateProviderPointer
		(
			const string& resultsFileName
		);

}


#endif // RESULTS_FACTORY_H