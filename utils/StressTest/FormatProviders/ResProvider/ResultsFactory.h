#ifndef RESULTS_FACTORY_H

#define RESULTS_FACTORY_H


#include "AbstractProviderRezr.h"


/*
* ������������ ���� ��� ������� �������� ���������� �������
* � ������ ����������� ������� MBS-��������
*/
namespace ResultsFactory
{

	/**
	* ������� ��������� ������� � ������ ����������� ������� MBS-��������
	* @param resultsFileName - ������������ ����� � ������������ ������� MBS-��������
	* @return ��������� ��������� ������� � ������ ����������� ������� MBS-��������
	*/
	AbstractProviderRezr* CreateProvider
		(
			const string& resultsFileName
		);

	/**
	* ������� ����� ��������� �� ��������� ������� � ������ ����������� ������� MBS-��������
	* @param resultsFileName - ������������ ����� � ������������ ������� MBS-��������
	* @return ����� ��������� �� ��������� ��������� ������� � ������ ����������� ������� MBS-��������
	*/
	AbstractProviderRezrPointer CreateProviderPointer
		(
			const string& resultsFileName
		);

}


#endif // RESULTS_FACTORY_H