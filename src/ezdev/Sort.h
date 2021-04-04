// QuickSort.h
// Ahmed M. Hussein (amhussein4@gmail.com)
// 05/05/2019

#ifndef SORT_H_
#define SORT_H_

namespace EZ
{
	namespace Algorithms
	{
		template<typename Type> unsigned int Partition(Type* array,unsigned const int& start_index,unsigned const int& end_index,const unsigned int& pivot_index)
		{
			// This function runs over the array and places the value at pivot_index in a new location such 
			// that all the array values less than this value are placed before it and all the array values 
			// more than this value are placed after it. This puts the value at its final sorted position. 
			// The function returns the index of that new position. 
			// get the pivot value (comparison value)
			Type pivot_value = array[pivot_index];
			// place it at the end of the array
			Type temp = pivot_value;
			array[pivot_index] = array[end_index];
			array[end_index] = temp;
			unsigned int new_pivot_index = start_index;
			for(unsigned int i = start_index ; i < end_index ; i++)
			{
				if(array[i] < pivot_value)
				{
					if(i != new_pivot_index)
					{
						// swap elements if required, new_pivot_index points at a value that 
						// is always assumed to be greater than the pivot value
						temp = array[new_pivot_index];
						array[new_pivot_index] = array[i];
						array[i] = temp;
					}
					new_pivot_index++;
				}
			}
			// put the pivot in its new location
			temp = array[new_pivot_index];
			array[new_pivot_index] = array[end_index];
			array[end_index] = temp;
			// return the new pivot index
			return new_pivot_index;
		}
		template<typename Type> unsigned int Partition(Type* array,unsigned int* sorted_indices,unsigned const int& start_index,unsigned const int& end_index,const unsigned int& pivot_index)
		{
			// This function runs over the array and places the value at pivot_index in a new location such 
			// that all the array values less than this value are placed before it and all the array values 
			// more than this value are placed after it. This puts the value at its final sorted position. 
			// The function returns the index of that new position. 
			// get the pivot value (comparison value) and its current index
			Type pivot_value = array[pivot_index];
			unsigned int pivot_value_index = sorted_indices[pivot_index];
			// place them at the end of their arrays
			Type temp = pivot_value;
			array[pivot_index] = array[end_index];
			array[end_index] = temp;
			unsigned int temp_index = pivot_value_index;
			sorted_indices[pivot_index] = sorted_indices[end_index];
			sorted_indices[end_index] = temp_index;
			unsigned int new_pivot_index = start_index;
			for(unsigned int i = start_index ; i < end_index ; i++)
			{
				if(array[i] < pivot_value)
				{
					if(i != new_pivot_index)
					{
						// swap elements if required, new_pivot_index points at a value that 
						// is always assumed to be greater than the pivot value
						temp = array[new_pivot_index];
						array[new_pivot_index] = array[i];
						array[i] = temp;
						pivot_value_index = sorted_indices[new_pivot_index];
						sorted_indices[new_pivot_index] = sorted_indices[i];
						sorted_indices[i] = pivot_value_index;
					}
					new_pivot_index++;
				}
			}
			// put the pivot in its new location
			temp = array[new_pivot_index];
			array[new_pivot_index] = array[end_index];
			array[end_index] = temp;
			temp_index = sorted_indices[new_pivot_index];
			sorted_indices[new_pivot_index] = sorted_indices[end_index];
			sorted_indices[end_index] = temp_index;
			// return the new pivot index
			return new_pivot_index;
		}
		template<typename Type> void QuickSort(Type* array,const unsigned int& length,unsigned const int& start_index,unsigned const int& end_index)
		{
			// check the problem validity
			if(length <= 1)		return;
			if(start_index >= end_index)	return;
			if(end_index >= length)			return;
			// pick a pivot, we choose the pivot to be the middle of the array
			unsigned int pivot_index = start_index + (end_index - start_index + 1)/2;
	 		// partition the array into 2 arrays based on the chosen pivot
	 		unsigned int new_pivot_index = Partition(array,start_index,end_index,pivot_index);
	 		// sort each array recursively
	 		QuickSort(array,length,start_index,new_pivot_index - 1);
	 		QuickSort(array,length,new_pivot_index + 1,end_index);
		}
		template<typename Type> void QuickSort(Type* array,const unsigned int& length,unsigned int* sorted_indices,unsigned const int& start_index,unsigned const int& end_index)
		{
			// check the problem validity
			if(length <= 1)		return;
			if(start_index >= end_index)	return;
			if(end_index >= length)			return;
			// pick a pivot, we choose the pivot to be the middle of the array
			unsigned int pivot_index = start_index + (end_index - start_index + 1)/2;
	 		// partition the array into 2 arrays based on the chosen pivot
	 		unsigned int new_pivot_index = Partition(array,sorted_indices,start_index,end_index,pivot_index);
	 		// sort each array recursively
	 		QuickSort(array,length,sorted_indices,start_index,new_pivot_index - 1);
	 		QuickSort(array,length,sorted_indices,new_pivot_index + 1,end_index);
		}
		template<typename Type> void QuickSort(Type* array,const unsigned int& length)
		{
			unsigned int end_index = length - 1;
			QuickSort(array,length,0,end_index);
		}
		template<typename Type> void QuickSort(Type* array,const unsigned int& length,unsigned int* sorted_indices)
		{
			unsigned int end_index = length - 1;
			for(unsigned int i = 0 ; i < length ; i++)
			{
				sorted_indices[i] = i;
			}
			QuickSort(array,length,sorted_indices,0,end_index);
		}
	}
}

#endif

