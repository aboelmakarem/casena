// Array.h
// Ahmed M. Hussein (amhussein4@gmail.com)
// 10/21/2019

#ifndef ARRAY_H_
#define ARRAY_H_

#include "string.h"
#include "Sort.h"

namespace EZ
{
	template<typename Type> class ElementOperator
	{
	public:
		ElementOperator(){}
		ElementOperator(const ElementOperator<Type>& target_operator){*this = target_operator;}
		virtual ~ElementOperator(){}
		virtual ElementOperator& operator=(const ElementOperator<Type>& target_operator){return *this;}
		virtual void operator()(Type element) const = 0;
	};

	template<typename Type> class Array
	{
	public:
		Array(){Initialize();}
		Array(const Array& array)
		{
			Initialize();
			*this = array;
		}
		Array(const unsigned int& allocation_size)
		{
			Initialize();
			Allocate(allocation_size);
		}
		~Array(){Reset();}
		Array& operator=(const Array& array)
		{
			Reset();
			Allocate(array.size);
			memcpy(elements,array.elements,size*sizeof(Type));
			return *this;
		}
		void Reset()
		{
			if(elements != 0)		delete [] elements;
			Initialize();
		}
		void Allocate(const unsigned int& allocation_size,const bool& set_to_zero = true)
		{
			if(size != allocation_size)
			{
				Reset();
				size = allocation_size;
				elements = new Type[size];
			}
			if(set_to_zero)	memset(elements,0,allocation_size*sizeof(Type));
		}
		unsigned int Size() const{return size;}
		Type& operator[](const unsigned int& index) const{return elements[index];}
		Type Item(const unsigned int& index) const{return elements[index];}
		bool Empty() const{return (size == 0);}
		Type First() const{return elements[0];}
		Type Last() const{return elements[size - 1];}
		void ForEachDelete()
		{
			for(unsigned int i = 0 ; i < size ; i++)
			{
				delete elements[i];
			}
			memset(elements,0,size*sizeof(Type));
		}
		void ForEachDeleteArray()
		{
			for(unsigned int i = 0 ; i < size ; i++)
			{
				delete [] elements[i];
			}
			memset(elements,0,size*sizeof(Type));
		}
		void ForEachApply(const ElementOperator<Type>& target_operator)
		{
			for(unsigned int i = 0 ; i < size ; i++)
			{
				target_operator(elements[i]);
			}
		}
		void Trim(const unsigned int& new_size)
		{
			// This function reduces the size of the array to the new_size if 
			// it is less than the current size and does nothing if it is not. 
			// It removes all elements beyond new_size (without deleting them) 
			// and the array loses all references to them while maintaining all 
			// elements up to new_size. 
			if(elements == 0)			return;
			if(new_size >= size)		return;
			Type* new_elements = new Type[new_size];
			memcpy(new_elements,elements,new_size*sizeof(Type));
			delete [] elements;
			elements = new_elements;
			size = new_size;
		}
		void Sort(){EZ::Algorithms::QuickSort(elements,size);}
		void Sort(unsigned int* sorted_indices){EZ::Algorithms::QuickSort(elements,size,sorted_indices);}
		void Unique()
		{
			EZ::Algorithms::QuickSort(elements,size);
			Type* new_elements = new Type[size];
			new_elements[0] = elements[0];
			unsigned int new_size = 1;
			for(unsigned int i = 1 ; i < size ; i++)
			{
				if(elements[i] != elements[i - 1])	new_elements[new_size++] = elements[i];
			}
			delete [] elements;
			elements = new Type[new_size];
			memcpy(elements,new_elements,new_size*sizeof(Type));
			delete [] new_elements;
			size = new_size;
		}
		void Compress(const Type& zero_element = 0)
		{
			// this function goes over all of the array elements and removes any elements 
			// which are equal to the passed zero_element, then it trims the array to 
			// its new size
			if(elements == 0)			return;
			unsigned int new_size = size;
			for(unsigned int i = 0 ; i < size ; i++)
			{
				if(elements[i] == zero_element)		new_size--;
			}
			if(size == new_size)		return;
			Type* new_elements = new Type[new_size];
			unsigned int index = 0;
			for(unsigned int i = 0 ; i < size ; i++)
			{
				if(elements[i] != zero_element)		new_elements[index++] = elements[i];
			}
			delete [] elements;
			elements = new_elements;
			size = new_size;
		}

	private:
		void Initialize()
		{
			elements = 0;
			size = 0;
		}
		Type* elements;
		unsigned int size;
	};

	template <typename Type> class Array2D
	{
	public:
		Array2D(){Initialize();}
		Array2D(const Array2D& array){*this = array;}
		Array2D(const unsigned int& target_row_count,const unsigned int& target_column_count)
		{
			Initialize();
			Allocate(target_row_count,target_column_count);
		}
		~Array2D(){Reset();}
		Array2D& operator=(const Array2D& array)
		{
			Reset();
			Allocate(array.row_count,array.column_count);
			for(unsigned int j = 0 ; j < row_count ; j++)
			{
				memcpy(elements[j],array.elements[j],column_count*sizeof(Type));
			}
			return *this;
		}
		void Reset()
		{
			ClearElements();
			Initialize();
		}
		void Allocate(const unsigned int& target_row_count,const unsigned int& target_column_count)
		{
			if((row_count != target_row_count) || (column_count != target_column_count))
			{
				Reset();
				row_count = target_row_count;
				column_count = target_column_count;
				elements = new Type*[row_count];
				for(unsigned int j = 0 ; j < row_count ; j++)
				{
					elements[j] = new Type[column_count];
					memset(elements[j],0,column_count*sizeof(Type));
				}
			}
		}
		unsigned int ColumnCount() const{return column_count;}
		unsigned int RowCount() const{return row_count;}
		Type& operator()(const unsigned int& row,const unsigned int& column) const{return elements[row][column];}
		bool Empty() const{return ((row_count == 0) || (column_count == 0));}
		Type RowFirst(const unsigned int& row) const{return elements[row][0];}
		Type RowLast(const unsigned int& row) const{return elements[row][column_count - 1];}
		Type ColumnFirst(const unsigned int& column) const{return elements[0][column];}
		Type ColumnLast(const unsigned int& column) const{return elements[row_count - 1][column];}
		void ForEachDelete()
		{
			for(unsigned int j = 0 ; j < row_count ; j++)
			{
				for(unsigned int i = 0 ; i < column_count ; i++)
				{
					delete elements[j][i];
				}
				memset(elements[j],0,column_count*sizeof(Type));
			}
		}
		void ForEachDeleteArray()
		{
			for(unsigned int j = 0 ; j < row_count ; j++)
			{
				for(unsigned int i = 0 ; i < column_count ; i++)
				{
					delete [] elements[j][i];
				}
				memset(elements[j],0,column_count*sizeof(Type));
			}
		}
		void ForEachApply(const ElementOperator<Type>& target_operator)
		{
			for(unsigned int j = 0 ; j < row_count ; j++)
			{
				for(unsigned int i = 0 ; i < column_count ; i++)
				{
					target_operator(elements[j][i]);
				}
			}
		}

	private:
		void Initialize()
		{
			column_count = 0;
			row_count = 0;
			elements = 0;
		}
		void ClearElements()
		{
			if(elements != 0)
			{
				for(unsigned int j = 0 ; j < row_count ; j++)
				{
					delete [] elements[j];
				}
				delete [] elements;
			}
		}
		Type** elements;
		unsigned int column_count;
		unsigned int row_count;
	};
}

#endif

