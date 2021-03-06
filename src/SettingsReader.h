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
#include <sstream>
#include <cmath>
#include <cstdlib>

#include <libconfig.h++>
#include <boost/filesystem.hpp>

class SettingsReader {
public:
	struct Parameter
	{
		Parameter(): m_name(""), m_value(0.0) {}
		Parameter(const std::string& name, double value): m_name(name), m_value(value) {}
		~Parameter() {}

		std::string m_name;
		double m_value;
	};
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
	bool readSettings(std::string filename="default.cfg");
	void initializeFitData(double * f) const;
	void initializeFitParameters(double * x) const;
	void initializeFitParametersBounds(double * lb, double * ub) const;
	void resetFitParameters(const double * x);
	void saveFitData(const double * f, std::string suffix);
	void saveFitParameters(const double * x, const double * cov, const double * info);
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

	void readEngineConfig(const libconfig::Setting& cfg);
	void readSampleConfig(const libconfig::Setting& cfg);
	void readFitParametersConfig(const libconfig::Setting& cfg);
	void readDataConfig(const libconfig::Setting& cfg);

	void registerSampleSetting(const libconfig::Setting& cfg);
	void registerDataSetting(const libconfig::Setting& reflection);
	bool readDataFile(const boost::filesystem::path& filename, DataFileProperty & dfp);

	bool isComment(std::string& str);

	bool saveReflections();//saves the reflection array into a file

	boost::filesystem::path engineCfgFile;
	boost::filesystem::path settingsCfgFile;
	boost::filesystem::path fitCfgFile;
	boost::filesystem::path resultFile;
	std::map<std::string, size_t> reflParametersIndices;
private:
	double ** calcParameters;//xs[irefl] - calculator parameters
	void resetCalcParameters();
	void allocateCalcParameters();
	void freeCalcParameters();
};

#endif /* SETTINGSREADER_H_ */
