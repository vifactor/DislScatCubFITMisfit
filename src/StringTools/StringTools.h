/*
 * StringTools.h
 *
 *  Created on: 30 jan. 2012
 *      Author: Viktor Kopp
 */

#ifndef STRINGTOOLS_H_
#define STRINGTOOLS_H_

#include<string>
#include <cctype>
//splits string into two parts:
//returns part after the last point (extension)
//argument transformed into part before las point(base name)
extern std::string stripExtension(std::string& filename);
extern std::string stripPath(std::string& filename);
extern std::size_t stripFirstBlanc(std::string& str);
extern std::size_t stripAllBlanc(std::string& str);

#endif /* STRINGTOOLS_H_ */
