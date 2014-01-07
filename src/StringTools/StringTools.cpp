/*
 * StringTools.cpp
 *
 *  Created on: 30 jan 2012
 *      Author: Viktor Kopp
 */

#include "StringTools.h"

std::string stripExtension(std::string& filename)
{
	std::string ext="";//no extension
	// find last '.' in the file
	std::string::size_type pos(filename.find_last_of('.'));
	if (pos != filename.npos)//has an extension
    {
		ext.assign(filename, pos+1, filename.npos);
		// put period into extension
		filename.assign(filename, 0, pos);
    }
	return ext;
}

//TODO Works only for Linux?
std::string stripPath(std::string& filename)
{
	std::string path="";//no path
	// find last '/' in the file
	std::string::size_type pos(filename.find_last_of('/'));
	if (pos != filename.npos)//has an extension
    {
		// put period into extension
		path.assign(filename, 0, pos);
		filename.assign(filename, pos+1, filename.npos);
    }
	return path;
}

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
