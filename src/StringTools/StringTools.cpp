/*
 * StringTools.cpp
 *
 *  Created on: 30 jan 2012
 *      Author: Viktor Kopp
 */

#include "StringTools.h"

std::size_t stripFirstBlanc(std::string& str)
{
	std::size_t nzeros=0;
	while(isspace(str[nzeros++]));
	str.erase(0, nzeros-1);
	return nzeros;
}

std::size_t stripAllBlanc(std::string& str)
{
	std::size_t nzeros=0;
	std::string temp("");
	for(std::size_t i=0; i<str.length(); i++)
	{
		if(isspace(str[i]))
	    {
			nzeros++;
	    }
		else
		{
			temp+=str[i];
		}
	}
	str=temp;
	return nzeros;
}
