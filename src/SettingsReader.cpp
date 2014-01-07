/*
 * SettingsReader.cpp
 *
 *  Created on: 13 dec 2011
 *      Author: kopp
 */

#include "SettingsReader.h"

SettingsReader::SettingsReader()
{
	lambdaX=1.54059;//Cu Kalpha wavelength
	aSub=5.43043;//Si lattice parameter
	//aSub=1.0;//Si lattice parameter
	nu=1.0/3;//'standard' Poisson ratio
	nbLayers=0;
	nbLinesSkip=0;
	calcParameters=NULL;
	background = 0.0;
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
	inpFileName=filename;
	LOG(logINFO)<<"Reading settings..."<<inpFileName;

	std::string dataFileExt;
	std::string dataDir;
	std::string layerPropsFile;
	try
	{
		ConfigFile config(inpFileName);

		config.readInto<std::string>(dataDir, "dataDirectory", "data");
		LOG(logINFO)<<"Data directory: \t"<<dataDir;
		dataDir="./"+dataDir+"/";

		config.readInto<std::string>(workDir, "workDirectory", "work");
		LOG(logINFO)<<"Working directory: \t"<<workDir;
		workDir="./"+workDir+"/";

		reflectionsFile=workDir+"reflections.dat";
		parametersFile=workDir+"parameters.dat";;
		errorsFile=workDir+"errors.dat";

		config.readInto<std::string>(layerPropsFile, "layerProperties", "default.prop");
		LOG(logINFO)<<"Layer stack properties: \t"<<layerPropsFile;
		layerPropsFile=workDir+layerPropsFile;

		config.readInto<std::string>(dataFileExt, "dataFileExtension", "csv");
		LOG(logINFO)<<"Data file extension: \t"<<dataFileExt;

		config.readInto<double>(background, "background", 0.0);
		LOG(logINFO) << "Background level : "<<background;

		config.readInto<double>(nu, "nu", 1.0/3);
		LOG(logINFO) << "Poisson ratio : "<<nu;

		config.readInto<size_t>(nbLinesSkip, "nbSkip", 0.0);
		LOG(logINFO) << "Skip every : "<<nbLinesSkip<<" datapoint";

		config.readInto<size_t>(nbIterations, "nbIterations", 10);
		LOG(logINFO) << "Perform  : "<<nbIterations<<" iterations";

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

bool SettingsReader::readStack(const std::string& filename)
{
	LOG(logINFO)<<"Reading layer stack parameters...";
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
	LOG(logINFO)<<"Layer stack consists of "<<nbLayers<<" layers";
	LOG(logINFO)<<"It is characterized by "<<layerParameterIndices.size()<<" parameters";
	LOG(logINFO)<<fitParameters.size()<<" of them are fit parameters";
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
	if (params.size() == 1)
	{
		layerParameterIndices.push_back(allParameters.size());
		allParameters.push_back(atof(params[0].c_str()));
	}
	else
	{
		layerParameterIndices.push_back(allParameters.size());

		fitParameter.parameterIndex=allParameters.size();
		fitParameter.lowValue=atof(params[1].c_str());
		fitParameter.upValue=atof(params[2].c_str());
		fitParameters.push_back(fitParameter);

		allParameters.push_back(atof(params[0].c_str()));
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

		LOG(logINFO)<<"Q0lab: "<<Q0xlab<<"\t"<<Q0zlab;

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
		allParameters.push_back(1.0);

		//sqf default value
		fitParameter.lowValue=1e-15;
		fitParameter.upValue=1e15;
		fitParameter.parameterIndex=reflection.sqfsIndices.at(ilayer);
		fitParameters.push_back(fitParameter);

		reflection.dQxIndices.at(ilayer)=allParameters.size();
		allParameters.push_back(0.0);//dQx default value
		reflection.dQzIndices.at(ilayer)=allParameters.size();
		allParameters.push_back(0.0);//dQz default value
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
			allParameters.at(reflection.sqfsIndices.at(ilayer))=sqf;
			allParameters.at(reflection.dQxIndices.at(ilayer))=dQx;
			allParameters.at(reflection.dQzIndices.at(ilayer))=dQz;
		}
		params.clear();
	}
	fin.close();
	return true;
}

bool SettingsReader::saveReflections()
{
	LOG(logINFO)<<"Saving reflections...";
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
				<<allParameters.at(reflParameters.at(irefl).sqfsIndices.at(ilayer))<<":"
				<<allParameters.at(reflParameters.at(irefl).dQxIndices.at(ilayer))<<":"
				<<allParameters.at(reflParameters.at(irefl).dQzIndices.at(ilayer))<<"\n";
		}
	}
	fout.close();
	return true;
}

bool SettingsReader::saveParameters()
{
	LOG(logINFO)<<"Saving parameters...";
	std::ofstream fout(parametersFile.c_str());
	if(!fout.is_open()) return false;

	time_t tt;
	time(&tt);
	fout<<"#"<<ctime(&tt);
	fout<<"#d\trho0\tk/60-deg\tsign\tg1\tg2\tQccz\n";
	for(size_t ilayer = 0; ilayer < nbLayers; ilayer++)
	{
		for(size_t iparam = 0; iparam < idxNB; iparam++)
			fout<<allParameters.at(layerParameterIndices.at(iparam+idxNB*ilayer))<<"\t";
		fout<<"\n";
	}
	fout.close();
	return true;
}

bool SettingsReader::readDataFile(const std::string& filename, DataFileProperty & dfp)
{
	LOG(logINFO)<<"File"<<filename<< " is being read...";
	std::ifstream fin(filename.c_str());
	if(!fin)
	{
		LOG(logERROR)<<"File is not found: "<<filename;
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
		LOG(logWARNING)<<"Check reflection in the file: "<<filename;
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

	LOG(logINFO)<<"Nb of Points to read:\t"<<nbPointsToRead;
	LOG(logINFO)<<"Nb of Points read "<<nbPointsRead;
	LOG(logINFO)<<"Total nbPoints "<<dataPoints.size();
	return true;
}

bool SettingsReader::readDataFiles(const std::string& dirname, const std::string& ext)
{
	LOG(logINFO)<<"Reading data files...";
	std::vector<std::string> files = std::vector<std::string>();
	std::string dir = dirname;// current directory

	DIR *dp;
	unsigned char isFile=0x8;
	struct dirent *dirp;
	if((dp = opendir(dir.c_str())) == NULL)
	{
		LOG(logERROR) << "Error opening directory: " << dir;
		return false;
	}

	while ((dirp = readdir(dp)) != NULL)
	{
		if(dirp->d_type != isFile)
			continue;
		std::string filename=dirp->d_name;
		std::string filenameext=stripExtension(filename);

		if(filenameext.compare(ext)!=0) continue;

		DataFileProperty fileContent;
		if(!readDataFile(dirname+filename+"."+ext, fileContent))
		{
			LOG(logWARNING)<<"File is empty or erroneous: "<<filename;
		}
		else
		{
			dataFileProperties.push_back(fileContent);
		}

	}
	closedir(dp);

	LOG(logINFO)<<"Nb of Files read "<<dataFileProperties.size();
	if(!dataPoints.size()) return false;
	return true;
}

void SettingsReader::initializeFitParameters(double * x) const
{
	for(size_t i=0;i<fitParameters.size();i++)
		x[i]=allParameters.at(fitParameters.at(i).parameterIndex);
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
		allParameters.at(fitParameters.at(i).parameterIndex)=x[i];
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
								+ idxNB * ilayer));
			}
			calcParameters[irefl][idxNB + 0 + (idxNB + 3) * ilayer]
					= allParameters.at(reflParameters.at(irefl).sqfsIndices.at(
							ilayer));
			calcParameters[irefl][idxNB + 1 + (idxNB + 3) * ilayer]
					= allParameters.at(reflParameters.at(irefl).dQxIndices.at(
							ilayer));
			calcParameters[irefl][idxNB + 2 + (idxNB + 3) * ilayer]
					= allParameters.at(reflParameters.at(irefl).dQzIndices.at(
							ilayer));

		}
	}
}

void SettingsReader::saveFitData(const double * f, std::string suffix)
{
	if(f)
	{
		LOG(logINFO) << "Writing the fitted data...";

		std::string filename;
		size_t ipoint=0;
		for (size_t ifile = 0;ifile < dataFileProperties.size();ifile++)
		{
			filename=dataFileProperties.at(ifile).filename;
			//stripPath(filename);
			stripExtension(filename);
			filename=workDir+filename+"_"+suffix+".out";

			LOG(logINFO) << "\t...to "<<filename;

			double qx, qz;
			double qper, qpar;
			double Qx=reflParameters.at(dataFileProperties.at(ifile).reflIndex).Qxlab;
			double Qz=reflParameters.at(dataFileProperties.at(ifile).reflIndex).Qzlab;;
			double nperx=reflParameters.at(dataFileProperties.at(ifile).reflIndex).nperx;
			double nperz=reflParameters.at(dataFileProperties.at(ifile).reflIndex).nperz;
			double nparx=reflParameters.at(dataFileProperties.at(ifile).reflIndex).nparx;
			double nparz=reflParameters.at(dataFileProperties.at(ifile).reflIndex).nparz;

			std::ofstream fout(filename.c_str());
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
	fout<<"#iparam\tidx\tvalue+/-error"<<std::endl;
	for(size_t iparam=0;iparam<getNbFitParameters();iparam++)
		fout<<iparam<<"\t"<<fitParameters.at(iparam).parameterIndex
						<<"\t"
						<<allParameters.at((fitParameters.at(iparam).parameterIndex))
						<<"\t+/-\t"
						<<sqrt(c[iparam+iparam*getNbFitParameters()])
						<<std::endl;
	fout.close();
}
