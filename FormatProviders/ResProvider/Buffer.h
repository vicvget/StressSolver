#ifndef BUFFER_H

#define BUFFER_H


//модуль для определения интерфейса класса буфера


#include "Types.h"


//класс буфера
class Buffer
{
public:

	// Конструктор и деструктор

	//конструктор
	Buffer
		(
			Integer size
		);

	//деструктор
	~Buffer();


	// Селекторы

	//получение указателя на буфер
	char* GetBuffer() const;

	//получение размера буфера
	Integer GetSize() const;

	//получение размера содержимого буфера
	Integer GetContentsSize() const;

	//возвращает истину, если буфер заполнен
	bool IsFilled() const;

	//возвращает истину, если буфер пуст
	bool IsEmpty() const;


	// Модификаторы

	//устанавливает новое значение размера содержимого
	void SetContentsLength
		(
			Integer contentsLength
		);

	//функции для заполнения/опустошения буфера

	//извлекает из буфера еще один блок
	bool Pop
		(
			void* block,
			Integer blockSize
		);

	//заносит в буфер еще один блок
	bool Push
		(
			const void* block,
			Integer blockSize
		);

	//опустошает буфер
	void Flush();


protected:

	//конструктор копирования
	Buffer
		(
			const Buffer&
		);


	// Поля

	char* _buffer; //указатель на буфер

	Integer _size; //размер буфера

	Integer _contentsLength; //размер содержимого буфера

	Integer _shift; //сдвиг начала области для чтения

};


#endif // BUFFER_H