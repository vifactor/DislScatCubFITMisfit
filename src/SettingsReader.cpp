/*
 * SettingsReader.cpp
 *
 *  Created on: 13 dec 2011
 *      Author: kopp
 */

#include "SettingsReader.h"
using namespace boost;

SettingsReader::SettingsReader()
{
	lambdaX = 1.54059; //Cu Kalpha wavelength
	aSub = 5.43043; //Si lattice parameter
	//aSub=1.0;//Si lattice parameter
	nu = 1.0 / 3;	//'standard' Poisson ratio
	nbLayers = 0;
	nbLinesSkip = 0;
	calcParameters = NULL;
	background = 0.0;
	nbIterations = 0;
}

SettingsReader::~SettingsReader()
{
	if (!reflParameters.empty())
		saveReflections();
	if(!layerParameterIndices.empty()) saveParameters();
	if(calcParameters) freeCalcParameters();
}

void SettingsReader::allocateCalcParameters()
{
	calcParameters = new double*[getNbReflections()];
	for (size_t irefl = 0; irefl < getNbReflections(); irefl++)
	{
		calcParameters[irefl] = new double[(idxNB + 3) * nbLayers];
	}
}

void SettingsReader::freeCalcParameters()
{
	for (size_t irefl = 0; irefl < getNbReflections(); irefl++)
	{
		delete[] calcParameters[irefl];
	}
	delete[] calcParameters;
}

bool SettingsReader::isComment(std::string& str)
{
	//empty string
	if(str.empty()) return true;
	//comment
	if(str[0]=='#') return true;
	//old reflection indicator
	if(str[0]=='%') return true;

	return false;
}

bool SettingsReader::readSettings(std::string cfgFile)
{
	boost::filesystem::path filename;
	libconfig::Config cfg;

	cfg.setAutoConvert(true);
	engineCfgFile = cfgFile;
	try
	{
		/*engine settings*/
		filename = engineCfgFile;
		cfg.readFile(filename.c_str());
		readEngineConfig(cfg.getRoot());

		/*engine settings*/
		filename = sampleCfgFile;
		cfg.readFile(filename.c_str());
		readSampleConfig(cfg.getRoot());

		/*data settings*/
		filename = dataCfgFile;
		cfg.readFile(filename.c_str());
		readDataConfig(cfg.getRoot());

		/*fit parameters settings*/
		filename = fitCfgFile;
		cfg.readFile(filename.c_str());
		readFitParametersConfig(cfg.getRoot());

	} catch (const libconfig::FileIOException &fioex)
	{
		std::cerr << fioex.what() << " in\t" << filename << std::endl;
		return false;
	} catch (const libconfig::ParseException &pex)
	{
		std::cerr << pex.what() << " in\t" << filename << ":"
						<< pex.getLine() << " - "
						<< pex.getError() << std::endl;
		return false;
	} catch (const libconfig::SettingNotFoundException &nfex)
	{
		std::cerr << nfex.what() << "\t" << nfex.getPath()
						<< " in\t" << filename << std::endl;
		return false;
	} catch (libconfig::SettingTypeException& tex)
	{
		std::cerr << tex.what() << "\t" << tex.getPath() << " in\t" << filename << std::endl;
		return false;
	}

	allocateCalcParameters();
	resetCalcParameters();


	std::cout << "fit parameters:" << std::endl;
	for (size_t i = 0; i < fitParameters.size(); ++i)
	{
		std::cout << allParameters.at(fitParameters.at(i).parameterIndex).m_name
				<< "\t=\t"
				<< allParameters.at(fitParameters.at(i).parameterIndex).m_value
				<< std::endl;
	}

	return true;
}

void SettingsReader::readEngineConfig(const libconfig::Setting& root)
{
	const libconfig::Setting& engineCfg = root["Engine"];

	/*output directory*/
	workDir = engineCfg["workDirectory"].c_str();
	std::cout << "Working directory: \t" << workDir << std::endl;

	/*file with sample properties*/
	sampleCfgFile = engineCfg["sample_cfg"].c_str();
	sampleCfgFile = workDir / sampleCfgFile;
	std::cout << "Sample configuration file: \t" << sampleCfgFile << std::endl;

	/*file with data configuration*/
	dataCfgFile = engineCfg["data_cfg"].c_str();
	dataCfgFile = workDir / dataCfgFile;
	std::cout << "Data configuration file: \t" << dataCfgFile << std::endl;

	/*file with fit parameters listing*/
	fitCfgFile = engineCfg["fit_cfg"].c_str();
	fitCfgFile = workDir / fitCfgFile;
	std::cout << "Fit parameters file: \t" << fitCfgFile << std::endl;

	/*file to write fitted parameters*/
	resultFile = engineCfg["output"].c_str();
	resultFile = workDir / resultFile;
	std::cout << "Output file: \t" << resultFile << std::endl;

	/*how many LM itaretions to perform*/
	nbIterations = engineCfg["nbIterations"];
	std::cout << "Perform  : " << nbIterations << " iterations" << std::endl
			<< std::endl;
}

void SettingsReader::readSampleConfig(const libconfig::Setting& root)
{
	const libconfig::Setting& sampleCfg = root["Sample"];

	nu = sampleCfg["nu"];
	std::cout << "Poisson ratio : "<<nu << std::endl;

	const libconfig::Setting& layersCfg = sampleCfg["layers"];
	std::cout<<"Reading layer stack parameters..." << std::endl;

	nbLayers = layersCfg.getLength();
	for(size_t ilay = 0; ilay < nbLayers; ++ilay)
	{
		registerSampleSetting(layersCfg[ilay]["d"]);
		registerSampleSetting(layersCfg[ilay]["rho0"]);
		registerSampleSetting(layersCfg[ilay]["k"]);
		registerSampleSetting(layersCfg[ilay]["sign"]);
		registerSampleSetting(layersCfg[ilay]["g1"]);
		registerSampleSetting(layersCfg[ilay]["g2"]);
		registerSampleSetting(layersCfg[ilay]["Qccz"]);
	}

	std::cout << "Layer stack consists of " << nbLayers << " layers"
			<< std::endl;
	std::cout << "It is characterized by " << layerParameterIndices.size()
			<< " parameters" << std::endl;
	std::cout << fitParameters.size() << " of them are fit parameters"
			<< std::endl;
	std::cout << std::endl;
}

void SettingsReader::readFitParametersConfig(const libconfig::Setting& cfg)
{

}

void SettingsReader::registerSampleSetting(const libconfig::Setting& stg)
{
	Parameter param;
	FitParameter fitParameter;

	layerParameterIndices.push_back(allParameters.size());

	param.m_name = stg.getPath();
	if(stg.isList())
	{
		param.m_value = stg[0];

		fitParameter.parameterIndex = allParameters.size();
		fitParameter.lowValue = stg[1][0];
		fitParameter.upValue= stg[1][1];
		fitParameters.push_back(fitParameter);
	}
	else
	{
		param.m_value = stg;
	}
	allParameters.push_back(param);
	std::cout << allParameters.back().m_name << "\t" << allParameters.back().m_value;
	if((!fitParameters.empty()) && (fitParameters.back().parameterIndex == allParameters.size() - 1))
		std::cout << "\tfit it!" << std::endl;
	else
		std::cout << std::endl;
}

void SettingsReader::readDataConfig(const libconfig::Setting& root)
{
	const libconfig::Setting& dataCfg = root["Data"];

	background = dataCfg["background"];
	std::cout << "Background level : "<<background<< std::endl;

	nbLinesSkip = dataCfg["nbSkip"];
	std::cout << "Skip every : "<<nbLinesSkip<<" datapoint" << std::endl;

	const libconfig::Setting& files = dataCfg["files"];
	/*loop over reflections*/
	for(int iref = 0; iref < files.getLength(); ++iref)
	{
		const libconfig::Setting& data = files[iref];
		registerDataSetting(data);
	}
}

void SettingsReader::registerDataSetting(const libconfig::Setting& reflection)
{
	const libconfig::Setting& Q = reflection["Q"];
	const libconfig::Setting& files = reflection["files"];
	const libconfig::Setting& layers = reflection["layers"];

	Reflection thisReflection;
	thisReflection.Qxcub = Q[0];
	thisReflection.Qycub = Q[1];
	thisReflection.Qzcub = Q[2];

	//transform the reflection coords from cubic to laboratory coordinate frame
	double Q0xlab = 1.0 / sqrt(2.0) * (thisReflection.Qxcub + thisReflection.Qycub);
	double Q0ylab = 1.0 / sqrt(2.0) * (-thisReflection.Qxcub + thisReflection.Qycub);
	double Q0zlab = thisReflection.Qzcub;

	std::cout<<"Q0lab: "<<Q0xlab<<"\t"<<Q0zlab << std::endl;

	if (Q0ylab != 0)
	{
		std::cerr << "Reflection " << Q.getPath()
				<< " is not usable since Q0y [lab] is equal to " << Q0ylab
				<< std::endl;
	}

	thisReflection.Qxlab=Q0xlab;
	thisReflection.Qzlab=Q0zlab;

	//reciprocal vector's norm
	double normQ = sqrt(Q0xlab * Q0xlab + Q0zlab * Q0zlab);
	//unity vectors parallel and perpendicular to the reciprocal vector
	thisReflection.nparx = Q0xlab / normQ;
	thisReflection.nparz = Q0zlab / normQ;
	thisReflection.nperx = 1.0
			/ sqrt(
					1.0
							+ thisReflection.nparx * thisReflection.nparx
									/ (thisReflection.nparz
											* thisReflection.nparz));
	thisReflection.nperz = -thisReflection.nparx / thisReflection.nparz
			* thisReflection.nperx;

	reflParametersIndices[Q.getPath()] = reflParameters.size();

	/*TODO check if nbLayers == layers.getLength()*/
	thisReflection.sqfsIndices.resize(nbLayers);
	thisReflection.dQxIndices.resize(nbLayers);
	thisReflection.dQzIndices.resize(nbLayers);

	FitParameter fitParameter;
	for (size_t ilayer = 0; ilayer < nbLayers; ilayer++)
	{
		thisReflection.sqfsIndices.at(ilayer) = allParameters.size();

		//sqf default value
		fitParameter.lowValue = 1e-15;
		fitParameter.upValue = 1e15;
		fitParameter.parameterIndex = thisReflection.sqfsIndices.at(ilayer);
		fitParameters.push_back(fitParameter);

		allParameters.push_back(Parameter("", 1.0));
		allParameters.at(thisReflection.sqfsIndices.at(ilayer)).m_name = layers[ilayer]["sqf"].getPath();
		allParameters.at(thisReflection.sqfsIndices.at(ilayer)).m_value = layers[ilayer]["sqf"];

		thisReflection.dQxIndices.at(ilayer) = allParameters.size();
		allParameters.push_back(Parameter("", 1.0));
		allParameters.at(thisReflection.dQxIndices.at(ilayer)).m_name = layers[ilayer]["center"][0].getPath();
		allParameters.at(thisReflection.dQxIndices.at(ilayer)).m_value = layers[ilayer]["center"][0];

		thisReflection.dQzIndices.at(ilayer) = allParameters.size();
		allParameters.push_back(Parameter("", 1.0));
		allParameters.at(thisReflection.dQzIndices.at(ilayer)).m_name = layers[ilayer]["center"][1].getPath();
		allParameters.at(thisReflection.dQzIndices.at(ilayer)).m_value = layers[ilayer]["center"][1];
	}
	reflParameters.push_back(thisReflection);

	for (int ifile = 0; ifile < files.getLength(); ifile++)
	{
		std::cout << Q.getPath() << "\t" << ifile << std::endl;
		DataFileProperty fileContent;
		fileContent.reflIndex = reflParametersIndices[Q.getPath()];
		if(readDataFile(files[ifile].c_str(), fileContent))
		{
			dataFileProperties.push_back(fileContent);
		}
		else
			std::cerr << "File is empty or erroneous: " << files[ifile].c_str() << std::endl;
	}
}

bool SettingsReader::saveReflections()
{
	boost::filesystem::path output_file;
	libconfig::Config cfg;

	// Read the file. If there is an error, report it and exit.
	try
	{
		cfg.readFile(dataCfgFile.c_str());
	} catch (const libconfig::FileIOException &fioex)
	{
		std::cerr << "I/O error while reading file." << std::endl;
		return (EXIT_FAILURE);
	} catch (const libconfig::ParseException &pex)
	{
		std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
				<< " - " << pex.getError() << std::endl;
		return (EXIT_FAILURE);
	}

	// Write out the updated configuration.
	output_file = dataCfgFile;
	output_file.replace_extension("~cfg");
	try
	{
		//reset configuration putting the new fitted values from fitParameters
		for(size_t i = 0; i < fitParameters.size(); ++i)
		{
			if(cfg.exists(allParameters.at(fitParameters.at(i).parameterIndex).m_name))
			{
				cfg.lookup(allParameters.at(fitParameters.at(i).parameterIndex).m_name) = allParameters.at(fitParameters.at(i).parameterIndex).m_value;
			}
		}

		cfg.writeFile(output_file.c_str());
		std::cout << "Updated configuration successfully written to: " << output_file
				<< std::endl;

	} catch (const libconfig::FileIOException &fioex)
	{
		std::cerr << "I/O error while writing file: " << output_file << std::endl;
		return (EXIT_FAILURE);
	}

	return (EXIT_SUCCESS);
}

bool SettingsReader::saveParameters()
{
	std::cout<<"Saving parameters..." << std::endl;
	std::ofstream fout("results.txt");
	if(!fout.is_open()) return false;

	time_t tt;
	time(&tt);
	fout<<"#"<<ctime(&tt);
	fout<<"#d\trho0\tk/60-deg\tsign\tg1\tg2\tQccz" << std::endl;
	for(size_t ilayer = 0; ilayer < nbLayers; ilayer++)
	{
		for (size_t iparam = 0; iparam < idxNB; iparam++)
			fout
					<< allParameters.at(
							layerParameterIndices.at(iparam + idxNB * ilayer)).m_value
					<< "\t";
		fout << std::endl;
	}
	fout.close();
	return true;
}

bool SettingsReader::readDataFile(const boost::filesystem::path& filename, DataFileProperty & dfp)
{
	std::cout<<"File "<< filename << " is being read..." << std::endl;
	std::ifstream fin(filename.c_str());
	if(!fin)
	{
		std::cerr << "File is not found: " << filename << std::endl;
		return false;
	}
	dfp.filename=filename;

	std::string line;
	std::istringstream is;

	size_t nbPointsToRead=0, nbPointsRead=0;
	DataPoint dp;

	while (!fin.eof())
	{
		getline(fin, line);
		//skip the comment
		if (isComment(line)) continue;
		//how many points could be used
		nbPointsToRead++;

		//TODO remove this when interpolation is implemented
		/*skip nbLinesSkip lines*/
		if (nbPointsToRead % (nbLinesSkip + 1) != 0)
			continue;
		//read data point
		is.str(line);
		is >> dp.qx >> dp.qz >> dp.intensity;
		dp.qx = dp.qx * aSub;
		dp.qz = dp.qz * aSub;

		dp.reflIndex = dfp.reflIndex;

		dataPoints.push_back(dp);
		nbPointsRead++;
		is.clear();
	}
	fin.close();
	dfp.nbPoints = nbPointsRead;

	std::cout << "Nb of Points to read:\t" << nbPointsToRead << std::endl;
	std::cout << "Nb of Points read " << nbPointsRead << std::endl;
	std::cout << "Total nbPoints " << dataPoints.size() << std::endl;
	return true;
}

void SettingsReader::initializeFitParameters(double * x) const
{
	for(size_t i=0;i<fitParameters.size();i++)
		x[i] = allParameters.at(fitParameters.at(i).parameterIndex).m_value;
}

void SettingsReader::initializeFitData(double * f) const
{
	for(size_t i=0;i<dataPoints.size();i++)
		f[i]=dataPoints.at(i).intensity;
}

void SettingsReader::initializeFitParametersBounds(double * lb, double * ub) const
{
	for(size_t i=0;i<fitParameters.size();i++)
	{
		lb[i]=fitParameters.at(i).lowValue;
		ub[i]=fitParameters.at(i).upValue;
	}
}

void SettingsReader::resetFitParameters(const double * x)
{
	for (size_t i = 0; i < fitParameters.size(); i++)
		allParameters.at(fitParameters.at(i).parameterIndex).m_value = x[i];
	resetCalcParameters();
}

const double * SettingsReader::getCalculatorParameters(int irefl) const
{
	return calcParameters[irefl];
}

void SettingsReader::resetCalcParameters()
{
	//IntensityCalculator indices
	//enum {idxD, idxRHO0, idxK, idxSIGN, idxG1, idxG2, idxDQCZ, idxSQF, idxDQX, idxDQZ, idxNB};//indices of parameters in parameter array
	//SettingsReader indices
	//enum {idxCOORD,idxRHO0,idxK,idxSIGN,idxG1, idxG2,idxQCCZ,idxNB};
	for (size_t irefl = 0; irefl < getNbReflections(); irefl++)
	{
		for (size_t ilayer = 0; ilayer < nbLayers; ilayer++)
		{
			//TODO change when IntensityCalculator is included
			for (size_t iparam = 0; iparam < idxNB; iparam++)
			{
				calcParameters[irefl][iparam + (idxNB + 3) * ilayer]
						= allParameters.at(layerParameterIndices.at(iparam
								+ idxNB * ilayer)).m_value;
			}
			calcParameters[irefl][idxNB + 0 + (idxNB + 3) * ilayer]
					= allParameters.at(reflParameters.at(irefl).sqfsIndices.at(
							ilayer)).m_value;
			calcParameters[irefl][idxNB + 1 + (idxNB + 3) * ilayer]
					= allParameters.at(reflParameters.at(irefl).dQxIndices.at(
							ilayer)).m_value;
			calcParameters[irefl][idxNB + 2 + (idxNB + 3) * ilayer]
					= allParameters.at(reflParameters.at(irefl).dQzIndices.at(
							ilayer)).m_value;

		}
	}
}

void SettingsReader::saveFitData(const double * f, std::string suffix)
{
	if(f)
	{
		std::cout << "Writing the fitted data..." << std::endl;

		filesystem::path outfile;
		size_t ipoint=0;
		for (size_t ifile = 0;ifile < dataFileProperties.size();ifile++)
		{
			/*trim old extension*/
			outfile = workDir / dataFileProperties.at(ifile).filename.stem();
			/*add new suffix and extention*/
			outfile = filesystem::path(outfile.native() + "_" + suffix + ".out");

			std::cout << "\t...to " << outfile << std::endl;

			double qx, qz;
			double qper, qpar;
			//double Qx=reflParameters.at(dataFileProperties.at(ifile).reflIndex).Qxlab;
			//double Qz=reflParameters.at(dataFileProperties.at(ifile).reflIndex).Qzlab;;
			double nperx=reflParameters.at(dataFileProperties.at(ifile).reflIndex).nperx;
			double nperz=reflParameters.at(dataFileProperties.at(ifile).reflIndex).nperz;
			double nparx=reflParameters.at(dataFileProperties.at(ifile).reflIndex).nparx;
			double nparz=reflParameters.at(dataFileProperties.at(ifile).reflIndex).nparz;

			std::ofstream fout(outfile.c_str());
			fout << "#qper\tqpar\tqx\tqz\tI_init\tI_fin" << std::endl;
			for(size_t ipt=0;ipt<dataFileProperties.at(ifile).nbPoints; ipt++)
			{
				qx=(dataPoints.at(ipoint).qx/(2*M_PI)/*+Qx*/)/aSub;
				qz=(dataPoints.at(ipoint).qz/(2*M_PI)/*+Qz*/)/aSub;
				qpar=nparx*qx+nparz*qz;
				qper=nperx*qx+nperz*qz;

				fout<<qper<<"\t"<<qpar<<"\t"<< qx <<"\t" << qz <<"\t"
						<<dataPoints.at(ipoint).intensity<<"\t"<<f[ipoint]<<"\n";
				ipoint++;
			};
			fout.close();
		}
	}
}

void SettingsReader::saveFitParameters(const double * x, const double * c)
{
	std::ofstream fout(resultFile.c_str());
	resetFitParameters(x);
	fout << "#iparam\tidx\tvalue+/-error" << std::endl;
	for (size_t iparam = 0; iparam < getNbFitParameters(); iparam++)
		fout << iparam << "\t" << fitParameters.at(iparam).parameterIndex
				<< "\t"
				<< allParameters.at((fitParameters.at(iparam).parameterIndex)).m_name
				<< "\t"
				<< allParameters.at((fitParameters.at(iparam).parameterIndex)).m_value
				<< "\t+/-\t" << sqrt(c[iparam + iparam * getNbFitParameters()])
				<< std::endl;
	fout.close();
}
