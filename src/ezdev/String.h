// String.h
// Ahmed M. Hussein (amhussein4@gmail.com)
// 06/12/2020

#ifndef STRING_H_
#define STRING_H_

#include "List.h"

namespace EZ
{
	class String
	{
	public:
		String();
		String(const char* target_string);
		String(const String& string);
		~String();
		String& operator=(const String& string);
		void Reset();
		void Set(const char* target_string);
		unsigned int Length() const;
		const char* operator()() const;
		String operator+(const char* string) const;
		String operator+(const String& string) const;
		bool operator==(const char* string) const;
		bool operator==(const String& string) const;
		bool operator<(const char* string) const;
		bool operator<(const String& string) const;
		bool operator>(const char* string) const;
		bool operator>(const String& string) const;
		bool operator<=(const char* string) const;
		bool operator<=(const String& string) const;
		bool operator>=(const char* string) const;
		bool operator>=(const String& string) const;
		void Tokenize(EZ::List<String*>& tokens,const char* delimiters) const;

	private:
		void Initialize();
		char* data;
		unsigned int length;
	};
}

#endif


