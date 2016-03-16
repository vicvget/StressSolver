#ifndef USING_POLICY_H

#define USING_POLICY_H


#include <vector>


using std::vector;


// список размеров подобластей сетки
typedef vector<int> SubdomainsSizesList;

// список признаков принадлежности подобластей сетке
typedef vector<bool> SubdomainsBelongingList;


/**
* Класс, содержащий перечисление типов политик удаления изолированных подобластей сетки
*/
struct DeletionPolicyTypes
{

	// перечисление типов политик удаления изолированных подобластей сетки
	enum DeletionPoliciesEnumeration
	{
		DeleteSmallSubdomains, // удалить подобласти, меньшие определенного размера  (количества узлов)
		DeleteBigSubdomains, // удалить подобласти, большие определенного размера (количества узлов)
		DeleteAllExceptBiggest // удалить все подобласти, кроме наибольшей по количеству узлов
	};

	// тип политики удаления изолированных подобластей сетки
	typedef DeletionPoliciesEnumeration PolicyType;

};

// тип политики удаления изолированных подобластей сетки
typedef DeletionPolicyTypes::PolicyType DeletionPolicyType;


/**
* Класс, содержащий список используемых политик удаления изолированных подобластей сетки
*/
struct UsingDeletionPolicies
{

	// список используемых политик удаления изолированных подобластей сетки
	static
	const DeletionPolicyType Policies[];

};


/**
* Класс, содержащий реализации политик удаления изолированных подобластей сетки
*/
template
	<
		DeletionPolicyType // тип политики удаления изолированных подобластей сетки
	>
struct DeletionPolicy;


// минимальный размер изолированных подобластей, которые следует оставлять при удалении по умолчанию
#define MINIMAL_SUBDOMAIN_SIZE 10


/**
* Класс, реализующий политику удаления маленьких изолированных подобластей
*/
template<>
struct DeletionPolicy<DeletionPolicyTypes::DeleteSmallSubdomains>
{

	// минимальный размер изолированных подобластей, которые следует оставлять при удалении
	static
	int MinimalSubdomainSize;

	/**
	* Применить политику удаления маленьких изолированных подобластей
	* @param sizesList - список размеров изолированных подобластей сетки
	* @param belongingList - список флагов принадлежности изолированных подобластей чистой сетке
	*/
	static
	void Apply
		(
			const SubdomainsSizesList& sizesList,
			SubdomainsBelongingList& belongingList
		);

};


/**
* Класс, реализующий политику удаления всех изолированных подобластей, кроме самой большой
*/
template<>
struct DeletionPolicy<DeletionPolicyTypes::DeleteAllExceptBiggest>
{

	/**
	* Применить политику удаления всех изолированных подобластей, кроме самой большой
	* @param sizesList - список размеров изолированных подобластей сетки
	* @param belongingList - список флагов принадлежности изолированных подобластей чистой сетке
	*/
	static
	void Apply
		(
			const SubdomainsSizesList& sizesList,
			SubdomainsBelongingList& belongingList
		);

};


#endif // USING_POLICY_H