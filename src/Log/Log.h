/*
 * Log.h
 *
 *  Created on: 26 jan 2011
 *  Changed on: 5 jan 2011
 *      Author: Dreamcatcher
 */

#ifndef LOG_H_
#define LOG_H_

#include <cstdio>
#include <sstream>
#include <ctime>
#include <cstring>

enum TLogLevel
{
	logERROR,
	logWARNING,
	logTIME,
	logINFO,
	logRTIME,
	logDEBUG
};

template <typename OutputPolicy>
class TLog {
public:
	TLog();
	virtual ~TLog();
	std::ostringstream& Get(TLogLevel level=logINFO);
	static TLogLevel & ReportingLevel();
protected:
	std::ostringstream os;
private:
	//TLog(const TLog&);
	//TLog& operator=(const TLog&);
	//returns time of program run
	std::string getRunTime();
	//returns time and date
	std::string getTime();
	std::string getLevel(TLogLevel level);
private:
	static TLogLevel messageLevel;
	static clock_t start;
};

//definitions of template members:
template <typename OutputPolicy>
TLogLevel TLog<OutputPolicy>::messageLevel=logINFO;

template <typename OutputPolicy>
clock_t TLog<OutputPolicy>::start=clock();

template <typename OutputPolicy>
TLog<OutputPolicy>::TLog() {
	// TODO Auto-generated constructor stub

}

template <typename OutputPolicy>
TLog<OutputPolicy>::~TLog()
{
	os<<std::endl;
	OutputPolicy::Output(os.str());
}

template <typename OutputPolicy>
std::string TLog<OutputPolicy>::getTime()
{
    char buf[256];

    time_t tt;
    time(&tt);
    strcpy(buf,ctime(&tt));

    //suppress the line break
    buf[strlen(buf)-1]='\0';
	std::string timeStr=buf;

    return timeStr;
}

template <typename OutputPolicy>
std::string TLog<OutputPolicy>::getRunTime()
{
	std::ostringstream ostmp;
	clock_t now=clock();
	ostmp<<"Time of run:\t"<<(now-start)<<" mcs";//CLOCKS_PER_SEC
    return ostmp.str();
}

template <typename OutputPolicy>
std::string TLog<OutputPolicy>::getLevel(TLogLevel level)
{
	switch(level)
	{
	case logRTIME:
		return getRunTime();
		break;
	case logTIME:
		return getTime();
		break;
	case logDEBUG:
		return "DEBUG";
		break;
	case logINFO:
		return "INFO";
		break;
	case logWARNING:
		return "WARNING";
		break;
	case logERROR:
		return "ERROR";
		break;
	default:
		break;
	}
	return "logUNKNOWN";
}

template <typename OutputPolicy>
std::ostringstream& TLog<OutputPolicy>::Get(TLogLevel level)
{
	//os<<"-"<<getTime();
	//Transform the logLevel to string
	os<<"-"<<getLevel(level);
	//generates level-logDEBUG+1 tab
	os<<" "<<std::string(level>logDEBUG ? 0:logDEBUG-level,' ');
	//messageLevel=level;
	return os;
}

template <typename OutputPolicy>
TLogLevel & TLog<OutputPolicy>::ReportingLevel()
{
	return messageLevel;
}

//implementation of OutputPolicy
class OutputTO
{
public:
	static FILE* & Stream();
	static void Output(const std::string& msg);
};

inline void OutputTO::Output(const std::string& msg)
{
	FILE* pStream=Stream();
	if(!pStream)
		return;
	fprintf(pStream,"%s",msg.c_str());
	fflush(pStream);
}

inline FILE* & OutputTO::Stream()
{
	static FILE* pStream=stderr;
	return pStream;
}

typedef TLog<OutputTO> Log;

#define LOG(level)\
if (level>Log::ReportingLevel()|| !OutputTO::Stream());\
else Log().Get(level)

#endif /* LOG_H_ */
