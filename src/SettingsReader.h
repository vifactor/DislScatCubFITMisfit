/*
 * SettingsReader.h
 *
 *  Created on: 13 dec 2011
 *      Author: kopp
 */

#ifndef SETTINGSREADER_H_
#define SETTINGSREADER_H_

//TODO rename class into settings processor

#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include <cstdlib>

#include <boost/filesystem.hpp>
//#include <dirent.h>

#include <Log.h>
#include <StringTools.h>
#include <ConfigFile.h>

class SettingsReader {
public:
	typedef double Parameter;
	typedef struct
	{
		size_t parameterIndex;//in the allParameters vector
		double lowValue;
		double upValue;
	} FitParameter;
	typedef struct
	{
		double Qxlab, Qzlab;
		double Qxcub, Qycub, Qzcub;
		double nparx, nparz;
		double nperx, nperz;
		std::vector<size_t> sqfsIndices;//in allParameters vector
		std::vector<size_t> dQxIndices;//in allParameters vector
		std::vector<size_t> dQzIndices;//in allParameters vector
	} Reflection;
	typedef struct
	{
		size_t reflIndex;
		double qx;
		double qz;
		double intensity;
	} DataPoint;
	//content of the data file
	typedef struct
	{
		/*which file the data is read from*/
		boost::filesystem::path filename;
		size_t reflIndex;
		size_t nbPoints;
	} DataFileProperty;
public:
	SettingsReader();
	virtual ~SettingsReader();
	bool readSettings(std::string filename="default.inp");
	void initializeFitData(double * f) const;
	void initializeFitParameters(double * x) const;
	void initializeFitParametersBounds(double * lb, double * ub) const;
	void resetFitParameters(const double * x);
	void saveFitData(const double * f, std::string suffix);
	void saveFitParameters(const double * x, const double * c);
	const double * getCalculatorParameters(int irefl) const;
	double getSubstrateLatticeParameter() const {return aSub;}
	double getBackground() const {return background;}
	double getPoissonRatio() const {return nu;}
	double getQxlab(size_t irefl) const {return reflParameters.at(irefl).Qxlab;}
	double getQzlab(size_t irefl) const {return reflParameters.at(irefl).Qzlab;}
	size_t getNbFitParameters() const {return fitParameters.size();}
	size_t getNbDataPoints() const {return dataPoints.size();}
	size_t getNbReflections() const {return reflParameters.size();}
	size_t getNbIterations() const {return nbIterations;}
	size_t getNbLayers() const {return nbLayers;}
	const DataPoint& getDataPoint(size_t ipt) const {return dataPoints.at(ipt);}
private:
	//content of the layer properties file
	enum {idxCOORD,idxRHO0,idxK,idxSIGN,idxG1, idxG2,idxQCCZ,idxNB};
	std::string inpFileName;
	boost::filesystem::path workDir;
	//content of the inp file
	size_t nbIterations;//how may iterations to perform before break
	size_t nbLinesSkip;//how many datapoints to skip in the dataFile
	std::string dataExt;//extension of datafile
	double lambdaX;//X-ray wavelength
	double aSub;//substrate lattice parameter
	double background;//intensity background
	double nu;//Poisson ratio (the same for all layers here)
	size_t nbLayers;
private:
	std::vector<Parameter> allParameters;
	std::vector<FitParameter> fitParameters;
	std::vector<Reflection> reflParameters;
	std::vector<size_t> layerParameterIndices;//in allParameters vector
	std::vector<DataPoint> dataPoints;
	std::vector<DataFileProperty> dataFileProperties;
	bool readStack(const boost::filesystem::path& filename);
	bool readDataFiles(const boost::filesystem::path& dirname, const std::string& ext);
	bool readDataFile(const boost::filesystem::path& filename, DataFileProperty & dfp);
	bool isComment(std::string& str);
	bool parseLine(std::string & str);
	void parseParameter(char * str, std::vector<std::string>& params);
	void registerParameter(std::vector<std::string>& params);
	int registerReflection(std::string);//returns index of the reflection in reflParameters
	bool saveReflections();//saves the reflection array into a file
	bool saveParameters();//saves the parameters array into a file
	bool setupSQFs(const std::string&, Reflection& );
	boost::filesystem::path reflectionsFile;
	boost::filesystem::path parametersFile;
	boost::filesystem::path errorsFile;
private:
	std::map<std::string, size_t> reflParametersIndices;//TODO this should be finally removed as ineffective
private:
	double ** calcParameters;//xs[irefl] - calculator parameters
	void resetCalcParameters();
	void allocateCalcParameters();
	void freeCalcParameters();
};

#endif /* SETTINGSREADER_H_ */
