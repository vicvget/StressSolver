#ifndef USING_POLICY_H

#define USING_POLICY_H


#include <vector>


using std::vector;


// ������ �������� ����������� �����
typedef vector<int> SubdomainsSizesList;

// ������ ��������� �������������� ����������� �����
typedef vector<bool> SubdomainsBelongingList;


/**
* �����, ���������� ������������ ����� ������� �������� ������������� ����������� �����
*/
struct DeletionPolicyTypes
{

	// ������������ ����� ������� �������� ������������� ����������� �����
	enum DeletionPoliciesEnumeration
	{
		DeleteSmallSubdomains, // ������� ����������, ������� ������������� �������  (���������� �����)
		DeleteBigSubdomains, // ������� ����������, ������� ������������� ������� (���������� �����)
		DeleteAllExceptBiggest // ������� ��� ����������, ����� ���������� �� ���������� �����
	};

	// ��� �������� �������� ������������� ����������� �����
	typedef DeletionPoliciesEnumeration PolicyType;

};

// ��� �������� �������� ������������� ����������� �����
typedef DeletionPolicyTypes::PolicyType DeletionPolicyType;


/**
* �����, ���������� ������ ������������ ������� �������� ������������� ����������� �����
*/
struct UsingDeletionPolicies
{

	// ������ ������������ ������� �������� ������������� ����������� �����
	static
	const DeletionPolicyType Policies[];

};


/**
* �����, ���������� ���������� ������� �������� ������������� ����������� �����
*/
template
	<
		DeletionPolicyType // ��� �������� �������� ������������� ����������� �����
	>
struct DeletionPolicy;


// ����������� ������ ������������� �����������, ������� ������� ��������� ��� �������� �� ���������
#define MINIMAL_SUBDOMAIN_SIZE 10


/**
* �����, ����������� �������� �������� ��������� ������������� �����������
*/
template<>
struct DeletionPolicy<DeletionPolicyTypes::DeleteSmallSubdomains>
{

	// ����������� ������ ������������� �����������, ������� ������� ��������� ��� ��������
	static
	int MinimalSubdomainSize;

	/**
	* ��������� �������� �������� ��������� ������������� �����������
	* @param sizesList - ������ �������� ������������� ����������� �����
	* @param belongingList - ������ ������ �������������� ������������� ����������� ������ �����
	*/
	static
	void Apply
		(
			const SubdomainsSizesList& sizesList,
			SubdomainsBelongingList& belongingList
		);

};


/**
* �����, ����������� �������� �������� ���� ������������� �����������, ����� ����� �������
*/
template<>
struct DeletionPolicy<DeletionPolicyTypes::DeleteAllExceptBiggest>
{

	/**
	* ��������� �������� �������� ���� ������������� �����������, ����� ����� �������
	* @param sizesList - ������ �������� ������������� ����������� �����
	* @param belongingList - ������ ������ �������������� ������������� ����������� ������ �����
	*/
	static
	void Apply
		(
			const SubdomainsSizesList& sizesList,
			SubdomainsBelongingList& belongingList
		);

};


#endif // USING_POLICY_H