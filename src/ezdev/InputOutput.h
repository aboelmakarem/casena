// InputOutput.h
// Ahmed M. Hussein (amhussein4@gmail.com)
// 04/30/2019

#ifndef INPUTOUTPUT_H_
#define INPUTOUTPUT_H_

#include "stdio.h"

namespace EZ
{
	namespace IO
	{
		bool SeekToRealLine(FILE* file);
		bool ReadLine(char* line,unsigned int size,FILE* file);
		bool SeekToNumber(FILE* file,unsigned int& number);
		void PrintLine(const char* line);
	}
}

#endif

