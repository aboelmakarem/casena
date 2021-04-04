// Dictionary.cpp
// Ahmed M. Hussein (amhussein4@gmail.com)
// 05/23/2019

#include "String.h"
#include "string.h"

namespace EZ
{
	String::String(){Initialize();}
	String::String(const char* target_string)
	{
		Initialize();
		Set(target_string);
	}
	String::String(const String& string){*this = string;}
	String::~String(){Reset();}
	String& String::operator=(const String& string)
	{
		Set(string.data);
		return *this;
	}
	void String::Reset()
	{
		if(data != 0)		delete [] data;
		Initialize();
	}
	void String::Set(const char* target_string)
	{
		if(strlen(target_string) != length)
		{
			if(data != 0)		delete [] data;
			length = strlen(target_string);
			data = new char[length + 1];
		}
		memcpy(data,target_string,length*sizeof(char));
		data[length] = 0;
	}
	unsigned int String::Length() const{return length;}
	const char* String::operator()() const{return data;}
	String String::operator+(const char* string) const
	{
		String result;
		unsigned int other_length = strlen(string);
		char* result_string = new char[length + other_length + 1];
		memcpy(result_string,data,length*sizeof(char));
		memcpy(&result_string[length],string,other_length*sizeof(char));
		result_string[length + other_length] = 0;
		result.Set(result_string);
		delete [] result_string;
		return result;
	}
	String String::operator+(const String& string) const{return (*this + string.data);}
	bool String::operator==(const char* string) const
	{
		// This function returns true if both strings are exactly the same and 
		// have the same lengths
		unsigned int i = 0;
		while(true)
		{
			if((data[i] == 0) && (string[i] == 0))	return true;
			if(data[i] == string[i])				continue;
			if(data[i] < string[i])					return false;
			if(data[i] > string[i])					return false;
			i++;
		}
		return false;
	}
	bool String::operator==(const String& string) const{return (*this == string());}
	bool String::operator<(const char* string) const
	{
		// A string A is less than string B if length of A is less than length of 
		// B or both A and B have the same length but a character in A has a lower 
		// ASCII value that its corresponding character in B
		unsigned int i = 0;
		while(true)
		{
			if((data[i] == 0) && (string[i] == 0))	return false;
			if(data[i] == string[i])				continue;
			if(data[i] < string[i])					return true;
			if(data[i] > string[i])					return false;
			i++;
		}
		return false;
	}
	bool String::operator<(const String& string) const{return (*this < string());}
	bool String::operator>(const char* string) const
	{
		if((*this < string))		return false;
		if((*this == string))		return false;
		return true;
	}
	bool String::operator>(const String& string) const{return (*this > string());}
	bool String::operator<=(const char* string) const{return ((*this < string) || (*this == string));}
	bool String::operator<=(const String& string) const{return (*this <= string());}
	bool String::operator>=(const char* string) const{return ((*this > string) || (*this == string));}
	bool String::operator>=(const String& string) const{return (*this >= string());}
	void String::Tokenize(EZ::List<String*>& tokens,const char* delimiters) const
	{
		unsigned int delimiter_count = strlen(delimiters);
		// 0 mode is scanning non-delimiters and 1 mode is scanning 
		// delimiters
		int mode = 0;
		// detect initial mode
		for(unsigned int j = 0 ; j < delimiter_count ; j++)
		{
			if(data[0] == delimiters[j])
			{
				mode = 1;
				break;
			}
		}
		bool is_delimiter = false;
		bool mode_switch = false;
		unsigned int last_break = 0;
		char* copy = new char[length + 1];
		for(unsigned int i = 1 ; i < length ; i++)
		{
			is_delimiter = false;
			for(unsigned int j = 0 ; j < delimiter_count ; j++)
			{
				if(data[i] == delimiters[j])
				{
					// delimiter found
					is_delimiter = true;
					break;
				}
			}
			mode_switch = false;
			if((mode == 0) && is_delimiter)
			{
				// character is a delimiter, switch mode
				mode = 1;
				mode_switch = true;
			}
			else if(!((mode == 0) || is_delimiter))
			{
				// character is not a delimiter, switch mode
				mode = 0;
				mode_switch = true;
			}
			if(mode_switch)
			{
				if(mode != 0)
				{
					memcpy(copy,&data[last_break],(i - last_break)*sizeof(char));
					copy[i - last_break] = 0;
					tokens.PushBack(new String(copy));
				}
				last_break = i;
			}
		}
		// always treat the end of the string as a delimiter
		if(mode == 0)
		{
			// if the last part of the string is a non-delimiter, then 
			// reaching the end of the string is a mode switch
			memcpy(copy,&data[last_break],(length - last_break)*sizeof(char));
			copy[length - last_break] = 0;
			tokens.PushBack(new String(copy));
		}
		delete [] copy;
	}
	void String::Initialize()
	{
		data = 0;
		length = 0;
	}
}

