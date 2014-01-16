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
	if(!reflParameters.empty()) saveReflections();
	if(!layerParameterIndices.empty()) saveParameters();
	if(calcParameters) freeCalcParameters();
}

void SettingsReader::allocateCalcParameters()
{
	LOG(logDEBUG)<<"Allocating calculator parameters";
	calcParameters=new double*[getNbReflections()];
	for(size_t irefl=0;irefl<getNbReflections();irefl++)
	{
		calcParameters[irefl]=new double [(idxNB+3)*nbLayers];
	}
}

void SettingsReader::freeCalcParameters()
{

	LOG(logDEBUG)<<"Freeing calculator parameters";
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
	return false;
}

bool SettingsReader::readSettings(std::string filename)
{
	inpFileName = filename;
	std::cout << "Reading settings..." << inpFileName << std::endl;

	std::string dataFileExt;
	filesystem::path dataDir;
	filesystem::path layerPropsFile;
	try
	{
		ConfigFile config(inpFileName);

		/*directory with data to fit*/
		m_sampleSettings.a0 = sample["a0"];
		config.readInto<filesystem::path>(dataDir, "dataDirectory", "data");
		std::cout<<"Data directory: \t"<< dataDir << std::endl;

		/*output directory*/
		config.readInto<filesystem::path>(workDir, "workDirectory", "work");
		std::cout<<"Working directory: \t"<<workDir << std::endl;

		reflectionsFile = workDir / "reflections.dat";
		parametersFile = workDir / "parameters.dat";
		errorsFile = workDir / "errors.dat";

		config.readInto<filesystem::path>(layerPropsFile, "layerProperties",
				"default.prop");
		std::cout << "Layer stack properties: \t" << layerPropsFile << std::endl;
		layerPropsFile = workDir / layerPropsFile;

		config.readInto<std::string>(dataFileExt, "dataFileExtension", "csv");
		std::cout<<"Data file extension: \t"<<dataFileExt << std::endl;

		config.readInto<double>(background, "background", 0.0);
		std::cout << "Background level : "<<background<< std::endl;

		config.readInto<double>(nu, "nu", 1.0/3);
		std::cout << "Poisson ratio : "<<nu << std::endl;

		config.readInto<size_t>(nbLinesSkip, "nbSkip", 0.0);
		std::cout << "Skip every : "<<nbLinesSkip<<" datapoint" << std::endl;

		config.readInto<size_t>(nbIterations, "nbIterations", 10);
		std::cout << "Perform  : " << nbIterations << " iterations" << std::endl
				<< std::endl;

	}
	catch(ConfigFile::file_not_found& e )
	{
		LOG(logERROR) <<"Error - File '" << e.filename << "' not found.";
		return false;
	}
	catch(ConfigFile::key_not_found& e )
	{
		LOG(logERROR) << "Error - Key '" << e.key << "' not found.";
		return false;
	}

	readStack(layerPropsFile);
	readDataFiles(dataDir, dataFileExt);
	allocateCalcParameters();
	resetCalcParameters();
	return true;
}

bool SettingsReader::readStack(const filesystem::path& filename)
{
	std::cout<<"Reading layer stack parameters..." << std::endl;
	std::ifstream fin(filename.c_str());
	if(!fin)
	{
		LOG(logERROR)<<"File is not found: "<<filename;
		return false;
	}
	std::istringstream is;
	std::string line;
	while(!fin.eof())
	{
			getline(fin,line);
			//skip the comment
			if(isComment(line)) continue;
			//read all data points;
			parseLine(line);
			nbLayers++;
	}
	fin.close();
	std::cout << "Layer stack consists of " << nbLayers << " layers"
			<< std::endl;
	std::cout << "It is characterized by " << layerParameterIndices.size()
			<< " parameters" << std::endl;
	std::cout << fitParameters.size() << " of them are fit parameters"
			<< std::endl;
	return true;
}

bool SettingsReader::parseLine(std::string & str)
{
	char cstr[str.size()+1];
	char delims[]=" \n\t\r";
	char * tok, * saveptr;
	size_t nbToks=0;
	std::vector<std::string> params;

	strcpy(cstr, str.c_str());
	for(tok=strtok_r(cstr,delims, &saveptr);tok;tok=strtok_r(NULL,delims, &saveptr))
	{
		parseParameter(tok, params);
		registerParameter(params);
		nbToks++;
		params.clear();
	}

	if(nbToks!=idxNB)
	{
		LOG(logWARNING)<<"Not enough parameters in the string:\n"<<str;
		return false;
	}

	return true;
}

void SettingsReader::parseParameter(char * str, std::vector<std::string>& params)
{
	char delims[]=":|";
	char * tok, *saveptr;

	for(tok=strtok_r(str,delims, &saveptr);tok;tok=strtok_r(NULL,delims, &saveptr))
	{
		params.push_back(tok);
	}
}

void SettingsReader::registerParameter(std::vector<std::string>& params)
{
	FitParameter fitParameter;
	Parameter parameter;
	if (params.size() == 1)
	{
		layerParameterIndices.push_back(allParameters.size());
		parameter.m_name = "";
		parameter.m_value = atof(params[0].c_str());
		allParameters.push_back(parameter);
	}
	else
	{
		layerParameterIndices.push_back(allParameters.size());

		fitParameter.parameterIndex=allParameters.size();
		fitParameter.lowValue=atof(params[1].c_str());
		fitParameter.upValue=atof(params[2].c_str());
		fitParameters.push_back(fitParameter);

		parameter.m_name = "";
		parameter.m_value = atof(params[0].c_str());
		allParameters.push_back(parameter);
	}
}

int SettingsReader::registerReflection(std::string str)
{
	stripAllBlanc(str);
	str=str.substr(1);//first symbol is '%'

	if (reflParametersIndices.find(str) == reflParametersIndices.end())
	{

		std::vector<std::string> params;
		char cstr[str.size() + 1];
		strcpy(cstr, str.c_str());
		parseParameter(cstr, params);
		if(params.size()<3)
		{
			return -1;
		}

		Reflection thisReflection;
		thisReflection.Qxcub = atof(params[0].c_str());
		thisReflection.Qycub = atof(params[1].c_str());
		thisReflection.Qzcub = atof(params[2].c_str());

		//transform the reflection coords from cubic to laboratory coordinate frame
		double Q0xlab = 1.0 / sqrt(2.0) * (thisReflection.Qxcub + thisReflection.Qycub);
		double Q0ylab = 1.0 / sqrt(2.0) * (-thisReflection.Qxcub + thisReflection.Qycub);
		double Q0zlab = thisReflection.Qzcub;

		std::cout<<"Q0lab: "<<Q0xlab<<"\t"<<Q0zlab << std::endl;

		if (Q0ylab != 0)
		{
			LOG(logWARNING)<<"Reflection "<<str<<" is not usable since Q0y [lab] is equal to "<<Q0ylab;
		}


		thisReflection.Qxlab=Q0xlab;
		thisReflection.Qzlab=Q0zlab;

		//reciprocal vector's norm
		double normQ=sqrt(Q0xlab*Q0xlab+Q0zlab*Q0zlab);
		//unity vectors parallel and perpendicular to the reciprocal vector
		thisReflection.nparx=Q0xlab/normQ;
		thisReflection.nparz=Q0zlab/normQ;
		thisReflection.nperx=1.0/sqrt(1.0+thisReflection.nparx*thisReflection.nparx/(thisReflection.nparz*thisReflection.nparz));
		thisReflection.nperz=-thisReflection.nparx/thisReflection.nparz*thisReflection.nperx;

		setupSQFs(str, thisReflection);

		reflParametersIndices[str]=reflParameters.size();
		reflParameters.push_back(thisReflection);
	}
	return reflParametersIndices.find(str)->second;
}

bool SettingsReader::setupSQFs(const std::string& str, Reflection& reflection)
{
	reflection.sqfsIndices.resize(nbLayers);
	reflection.dQxIndices.resize(nbLayers);
	reflection.dQzIndices.resize(nbLayers);
	std::ifstream fin(reflectionsFile.c_str());

	FitParameter fitParameter;
	for (size_t ilayer = 0; ilayer < nbLayers; ilayer++)
	{
		reflection.sqfsIndices.at(ilayer)=allParameters.size();
		allParameters.push_back(
				Parameter("", 1.0)
				);

		//sqf default value
		fitParameter.lowValue=1e-15;
		fitParameter.upValue=1e15;
		fitParameter.parameterIndex=reflection.sqfsIndices.at(ilayer);
		fitParameters.push_back(fitParameter);

		reflection.dQxIndices.at(ilayer)=allParameters.size();
		allParameters.push_back(
				Parameter("", 0.0)
				);//dQx default value
		reflection.dQzIndices.at(ilayer)=allParameters.size();
		allParameters.push_back(
				Parameter("", 0.0)
				);//dQz default value
	}
	if(!fin.is_open())	return false;
	std::string line;
	std::vector<std::string> params;
	int ilayer;
	double sqf, dQx, dQz;
	while(!fin.eof())
	{
		getline(fin,line);
		//skip the comment
		if(isComment(line)) continue;
		//find reflection
		if(line.find(str)!=std::string::npos)
		{
			char tok[line.size()];
			strcpy(tok, line.c_str());
			parseParameter(tok, params);//"Qx:Qy:Qz:ilayer:sqf:dQx:dQz"
			if(params.size()<7) continue;

			ilayer=atoi(params[3].c_str());
			sqf=atof(params[4].c_str());
			dQx=atof(params[5].c_str());
			dQz=atof(params[6].c_str());
			allParameters.at(reflection.sqfsIndices.at(ilayer)).m_value = sqf;
			allParameters.at(reflection.dQxIndices.at(ilayer)).m_value = dQx;
			allParameters.at(reflection.dQzIndices.at(ilayer)).m_value = dQz;
		}
		params.clear();
	}
	fin.close();
	return true;
}

bool SettingsReader::saveReflections()
{
	std::cout<<"Saving reflections..." << std::endl;
	std::ofstream fout(reflectionsFile.c_str());
	if(!fout.is_open()) return false;
	for(size_t irefl = 0; irefl < reflParameters.size(); irefl++)
	{
		for(size_t ilayer = 0; ilayer < nbLayers; ilayer++)
		{
			fout<<reflParameters.at(irefl).Qxcub<<":"
				<<reflParameters.at(irefl).Qycub<<":"
				<<reflParameters.at(irefl).Qzcub<<":"
				<<ilayer<<":"
				<<allParameters.at(reflParameters.at(irefl).sqfsIndices.at(ilayer)).m_value<<":"
				<<allParameters.at(reflParameters.at(irefl).dQxIndices.at(ilayer)).m_value<<":"
				<<allParameters.at(reflParameters.at(irefl).dQzIndices.at(ilayer)).m_value<<"\n";
		}
	}
	fout.close();
	return true;
}

bool SettingsReader::saveParameters()
{
	std::cout<<"Saving parameters..." << std::endl;
	std::ofstream fout(parametersFile.c_str());
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
		std::cerr << "File is not found: " << filename;
		return false;
	}
	dfp.filename=filename;

	std::string line;
	int reflIndex=-1;
	std::istringstream is;

	size_t nbPointsToRead=0, nbPointsRead=0;
	DataPoint dp;

	while(!fin.eof())
	{
		getline(fin,line);
		if(line.find("%")!=std::string::npos)
		{
			reflIndex=registerReflection(line);
			break;
		}
	}
	if(reflIndex<0)
	{
		LOG(logWARNING)<<"Check reflection in the file:\t"<<filename;
		return false;
	}
	dfp.reflIndex=reflIndex;
	while(!fin.eof())
	{
		getline(fin,line);
		//skip the comment
		if(isComment(line)) continue;
		//how many points could be used
		nbPointsToRead++;
		//skip nbLinesSkip lines
		if(nbPointsToRead%(nbLinesSkip+1)!=0)	continue;
		//read data point
		is.str(line);
		is>>dp.qx>>dp.qz>>dp.intensity;
		dp.qx = dp.qx*aSub;
		dp.qz = dp.qz*aSub;

		dp.reflIndex=reflIndex;

		dataPoints.push_back(dp);
		nbPointsRead++;
		is.clear();
	}
	fin.close();
	dfp.nbPoints=nbPointsRead;

	std::cout << "Nb of Points to read:\t" << nbPointsToRead << std::endl;
	std::cout << "Nb of Points read " << nbPointsRead << std::endl;
	std::cout << "Total nbPoints " << dataPoints.size() << std::endl;
	return true;
}

bool SettingsReader::readDataFiles(const boost::filesystem::path& dirname, const std::string& ext)
{
	std::cout<<"Reading data files...";
	std::vector<std::string> files = std::vector<std::string>();

	filesystem::path path(dirname);
	try
	{
		if ( filesystem::is_directory( path ) )
		{
			filesystem::directory_iterator end_iter;
			for (filesystem::directory_iterator dir_itr(path);
					dir_itr != end_iter; ++dir_itr)
			{
				if (filesystem::is_regular_file(dir_itr->status()))
				{
					filesystem::path filename = dir_itr->path();
					std::string extention = dir_itr->path().extension().c_str();

					if(extention.compare("." + ext) == 0)
					{
						DataFileProperty fileContent;
						if(readDataFile(filename, fileContent))
						{
							dataFileProperties.push_back(fileContent);
						}
						else
							std::cout << "File is empty or erroneous: " << dir_itr->path().filename();
					}

				}
			}
		}
		else
		{
			std::cout << "Error opening directory:\t" << path;
			return false;
		}
	} catch(const std::exception & ex)
	{
		std::cout << path.c_str() << " " << ex.what() << std::endl;
	}

	std::cout << "Nb files read\t" << dataFileProperties.size() << std::endl;
	if (!dataPoints.size())
		return false;
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
	for(size_t i=0;i<fitParameters.size();i++)
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
	std::ofstream fout(errorsFile.c_str());
	resetFitParameters(x);
	fout << "#iparam\tidx\tvalue+/-error" << std::endl;
	for (size_t iparam = 0; iparam < getNbFitParameters(); iparam++)
		fout << iparam << "\t" << fitParameters.at(iparam).parameterIndex
				<< "\t"
				<< allParameters.at((fitParameters.at(iparam).parameterIndex)).m_value
				<< "\t+/-\t" << sqrt(c[iparam + iparam * getNbFitParameters()])
				<< std::endl;
	fout.close();
}
