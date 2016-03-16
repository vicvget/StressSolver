#include "UsingPolicy.h"

#include <algorithm>


// struct DeletionPolicyTypes

// ������ ������������ ������� �������� ������������� ����������� �����
// static
const DeletionPolicyType UsingDeletionPolicies::Policies[] =
	{
		DeletionPolicyTypes::DeleteSmallSubdomains,
		DeletionPolicyTypes::DeleteAllExceptBiggest
	};


// struct DeletionPolicy


// struct DeletionPolicy<DeletionPolicyTypes::DeleteSmallSubdomains>

// ����������� ������ ������������� �����������, ������� ������� ��������� ��� ��������
// static
int DeletionPolicy<DeletionPolicyTypes::DeleteSmallSubdomains>::MinimalSubdomainSize =
	MINIMAL_SUBDOMAIN_SIZE;

/**
* ��������� �������� �������� ��������� ������������� �����������
* @param sizesList - ������ �������� ������������� ����������� �����
* @param belongingList - ������ ������ �������������� ������������� ����������� ������ �����
*/
// static
void DeletionPolicy<DeletionPolicyTypes::DeleteSmallSubdomains>::Apply
	(
		const SubdomainsSizesList& sizesList,
		SubdomainsBelongingList& belongingList
	)
{
	if (sizesList.size() != belongingList.size())
	{
		return;
	}

	SubdomainsSizesList::const_iterator sizesIterator;
	SubdomainsBelongingList::iterator belongingIterator;

	for
		(
			sizesIterator = sizesList.begin(),
			belongingIterator = belongingList.begin();
			sizesIterator != sizesList.end();
			++sizesIterator,
			++belongingIterator
		)
	{
		*belongingIterator = (*sizesIterator >= MinimalSubdomainSize);
	}
}


// struct DeletionPolicy<DeletionPolicyTypes::DeleteAllExceptBiggest>

/**
* ��������� �������� �������� ���� ������������� �����������, ����� ����� �������
* @param sizesList - ������ �������� ������������� ����������� �����
* @param belongingList - ������ ������ �������������� ������������� ����������� ������ �����
*/
// static
void DeletionPolicy<DeletionPolicyTypes::DeleteAllExceptBiggest>::Apply
	(
		const SubdomainsSizesList& sizesList,
		SubdomainsBelongingList& belongingList
	)
{
	if (sizesList.size() != belongingList.size())
	{
		return;
	}

	SubdomainsSizesList::const_iterator maxSizeIterator;
	maxSizeIterator = std::max_element(sizesList.begin(), sizesList.end());
	int maxSize = *maxSizeIterator;
	SubdomainsSizesList::const_iterator sizesIterator;
	SubdomainsBelongingList::iterator belongingIterator;

	std::fill(belongingList.begin(), belongingList.end(), false);
	for
		(
			sizesIterator = sizesList.begin(),
			belongingIterator = belongingList.begin();
			sizesIterator != sizesList.end();
			++sizesIterator,
			++belongingIterator
		)
	{
		if (*sizesIterator == maxSize)
		{
			*belongingIterator = true;
		}
	}
}