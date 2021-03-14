// List.h
// Ahmed M. Hussein (amhussein4@gmail.com)
// 07/10/2019

#ifndef LIST_H_
#define LIST_H_

namespace EZ
{
	template<typename Type> class ListItem
	{
	public:
		ListItem(){Initialize();}
		~ListItem(){Reset();}
		void Reset(){Initialize();}
		ListItem<Type>* Previous() const{return this->operator--();}
		void Previous(ListItem<Type>* item){previous = item;}
		ListItem<Type>* Next() const{return this->operator++();}
		void Next(ListItem<Type>* item){next = item;}
		Type operator()() const{return data;}
		void operator()(const Type& target_data){data = target_data;}
		ListItem<Type>* operator++() const{return next;}
		ListItem<Type>* operator--() const{return previous;}
		Type Data() const{return data;}
		Type* DataPointer(){return &data;}
		void Data(Type target_data){data = target_data;}

	private:
		ListItem(const ListItem& item){*this = item;}
		ListItem& operator=(const ListItem& item){return *this;}
		void Initialize()
		{
			previous = 0;
			next = 0;
		}
		Type data;
		ListItem<Type>* previous;
		ListItem<Type>* next;
	};
	
	template<typename Type> class ListItemOperator
	{
	public:
		ListItemOperator(){}
		ListItemOperator(const ListItemOperator<Type>& target_operator){*this = target_operator;}
		virtual ~ListItemOperator(){}
		virtual ListItemOperator& operator=(const ListItemOperator<Type>& target_operator){return *this;}
		virtual void operator()(ListItem<Type>* item) = 0;
	};

	template<typename Type> class List
	{
	public:
		List(){Initialize();}
		List(const List& list){*this = list;}
		~List(){Reset();}
		List& operator=(const List& list)
		{
			Reset();
			ListItem<Type>* list_item = list.start;
			while(list_item != 0)
			{
				PushBack(list_item->Data());
				list_item = list_item->Next();
			}
			return *this;
		}
		void Reset()
		{
			while(size > 0)
			{
				PopFront();
			}
			Initialize();
		}
		unsigned int Size() const{return size;}
		ListItem<Type>* InsertBefore(const Type& data,ListItem<Type>* item)
		{
			// This function inserts an item to the list before the item 
			// pointed at by the item argument. If the list is empty, 
			// it inserts it at the beginning regardless of where item 
			// points. It returns a pointer to the newly inserted item. 
			if((size > 0) && (item == 0))		return 0;
			ListItem<Type>* new_item = new ListItem<Type>;
			new_item->Data(data);
			if(size == 0)
			{
				// empty list, this is the first item to be inserted
				start = new_item;
				end = new_item;
			}
			else
			{
				new_item->Previous(item->Previous());
				new_item->Next(item);
				item->Previous(new_item);
				if(new_item->Previous() != 0)		new_item->Previous()->Next(new_item);
				else 								start = new_item;
			}
			size++;
			return new_item;
		}
		ListItem<Type>* InsertAfter(const Type& data,ListItem<Type>* item)
		{
			// This function inserts an item to the list after the item 
			// pointed at by the item argument. If the list is empty, 
			// it inserts it at the beginning regardless of where item 
			// points. It returns a pointer to the newly inserted item. 
			if((size > 0) && (item == 0))		return 0;
			ListItem<Type>* new_item = new ListItem<Type>;
			new_item->Data(data);
			if(size == 0)
			{
				// empty list, this is the first item to be inserted
				start = new_item;
				end = new_item;
			}
			else
			{
				new_item->Next(item->Next());
				new_item->Previous(item);
				item->Next(new_item);
				if(new_item->Next() != 0)		new_item->Next()->Previous(new_item);
				else 							end = new_item;
			}
			size++;
			return new_item;
		}
		void PushBack(const Type& data){end = InsertAfter(data,end);}
		void PushFront(const Type& data){start = InsertBefore(data,start);}
		ListItem<Type>* DropItem(const ListItem<Type>* item)
		{
			// This function drops the item pointed at by the passed item from the 
			// list and returns a pointer to the item following it if any. 
			if(item == 0)			return 0;
			if(size == 0)			return 0;
			if(size == 1)
			{
				// the list has only one item and the item points to it
				if(start == 0)		return 0;
				if(start != end)	return 0;
				if(start != item)	return 0;
				delete item;
				Initialize();
				return 0;
			}
			ListItem<Type>* previous = item->Previous();
			ListItem<Type>* next = item->Next();
			delete item;
			if(previous == 0)
			{
				next->Previous(0);
				start = next;
			}
			else if(next == 0)
			{
				previous->Next(0);
				end = previous;
			}
			else
			{
				next->Previous(previous);
				previous->Next(next);
			}
			size--;
			return next;
		}
		void PopBack(){DropItem(end);}
		void PopFront(){DropItem(start);}
		ListItem<Type>* Start() const{return start;}
		ListItem<Type>* End() const{return end;}
		bool Empty() const{return (size == 0);}
		Type Front() const{return start->Data();}
		Type Back() const{return end->Data();}
		void ForEachDelete()
		{
			for(ListItem<Type>* item = start ; item != 0 ; item = item->Next())
			{
				delete item->Data();
			}
		}
		void ForEachApply(ListItemOperator<Type>& target_operator)
		{
			for(ListItem<Type>* item = start ; item != 0 ; item = item->Next())
			{
				target_operator(item);
			}
		}
		ListItem<Type>* CircularNext(const ListItem<Type>* item) const
		{
			if(item->Next() != 0)		return item->Next();
			if(item != end)				return 0;
			return start;
		}
		ListItem<Type>* CircularPrevious(const ListItem<Type>* item) const
		{
			if(item->Previous() != 0)		return item->Previous();
			if(item != start)				return 0;
			return end;
		}

	private:
		void Initialize()
		{
			start = 0;
			end = 0;
			size = 0;
		}
		ListItem<Type>* start;
		ListItem<Type>* end;
		unsigned int size;
	};
}

#endif

