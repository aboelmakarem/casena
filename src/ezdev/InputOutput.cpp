// InputOutput.cpp
// Ahmed M. Hussein (amhussein4@gmail.com)
// 04/30/2019

#include "InputOutput.h"
#include "string.h"

namespace EZ
{
	namespace IO
	{
		bool SeekToRealLine(FILE* file)
		{
			// This function will scan the file stream until a non-blank character is found. 
			// If found, the file pointer will be set so that it points at that first non-blank 
			// character and it will return true. Otherwise, it will return false. 
			int character = 0;
			while(true)
			{
				if(feof(file))			return false;
				character = fgetc(file);
				// skip all invalid characters
				if(character > 32)
				{
					// skip all comment characters (asterisks and hash symbols), if a 
					// comment character is found, skip to the end of this line, or to 
					// the end of the file if it is the last line of the file
					if((character == 35) || (character == 42))
					{
						while(true)
						{
							if(feof(file))			return false;
							character = fgetc(file);
							if((character == 10) || (character == 13))		break;
						}
					}
					else
					{
						// found a valid line start, move the file pointer one character 
						// back
						ungetc(character,file);
						return true;
					}
				}
			}
		}
		bool ReadLine(char* line,unsigned int size,FILE* file)
		{
			if(!SeekToRealLine(file))			return false;
			if(fgets(line,size,file) != 0)		return true;
			return false;
		}
		bool SeekToNumber(FILE* file,unsigned int& number)
		{
			// This function keeps seeking through the file until it finds a line that 
			// starts with a number then it returns it. The file pointer at this point points 
			// at the line following that with the number
			char read_line[128];
			number = 0;
			while(true)
			{
				if(!ReadLine(read_line,128,file))		return false;
				// see if the read line starts with a number
				if((read_line[0] >= 48) && (read_line[0] <= 57))
				{
					sscanf(read_line,"%u",&number);
					break;
				}
			}
			return true;
		}
		void PrintLine(const char* line)
		{
			printf("%s\n",line);
			fflush(0);
		}
	}
}

