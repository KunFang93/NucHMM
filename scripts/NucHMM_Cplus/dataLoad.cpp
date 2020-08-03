#include "dataLoad.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <list>
#include <sstream>
#include <string>
#include <vector>

#include "constants.hpp"
#include "macros.hpp"
#include "options.hpp"

using std::atoi;
using std::ceil;
using std::cerr;
using std::endl;
using std::cout;
using std::exit;
using std::exp;
using std::fabs;
using std::ifstream;
using std::list;
using std::ostream;
using std::pow;
using std::size_t;
using std::strcpy;
using std::string;
using std::stringstream;
using std::strtoul;
using std::swap;
using std::vector;

namespace HMM
{
	int** dataLoad::readFile(const string& fname, int interval, ostream& outs, int readShift, bool useFivePrime)
	{
		unsigned long* chrLengths = options::getOptions()->chrLengths;
		if (readShift == -1) readShift = interval / 2;
		int** m = new int*[options::getOptions()->chromosomes];
		for (int i = 0; i < options::getOptions()->chromosomes; i++)
		{
			unsigned long numBins = (chrLengths[i] / interval) + 1;
			int* arr = new int[numBins];
			for (unsigned long j=0; j<numBins; j++)
			{
				arr[j] = 0;
			}
			m[i] = arr;
		}
		string chr;
		string junk;
		unsigned long bpLoc;
		unsigned long bin;
		ifstream file(fname.c_str());
		unsigned long lineNum = 0;
		if (file.is_open())
		{
			string line;
			while(!! getline(file,line))
			{
			    lineNum++;
				if (!line.empty())
				{
					//string buf;
					//stringstream ss(line);
					vector<string> tokens;
					// inspired by http://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
					string delim = "\t";
					size_t start = 0;
					size_t end = line.find(delim);
					while (end != string::npos)
					{
						tokens.push_back(line.substr(start, end - start));
						start = end + delim.length();
						end = line.find(delim, start);
					}
					tokens.push_back(line.substr(start, end));
					//while (ss >> buf) tokens.push_back(buf);
					if (tokens.size() >= 6)
					{
						chr = tokens[0];
						int chrNum = -1; //so it'll be loud and clear with a segfault if there's something wrong here
						if (options::getOptions()->chrNames.count(chr) == 0)
						{
#ifdef WARN_INVALID_CHROMOSOME
                            outs << "File " << fname << ", line " << lineNum << ": Invalid chromosome " << chr << ", read excluded."endl;
#endif
                            continue;
						}
						else
						{
							chrNum = options::getOptions()->chrNames[chr];
						}
						unsigned long lower = strtoul(tokens[1].c_str(),NULL,0);
						unsigned long higher = strtoul(tokens[2].c_str(),NULL,0);
						if (lower > higher)
						{
							unsigned long temp = lower;
							lower = higher;
							higher = temp;
						}
						if (tokens[1].c_str()[0] == '-') //lower is unsigned, can't check < 0
						{
							outs << "File " << fname << ", line " << lineNum << ": Read before beginning of chromosome, read excluded." << endl;
							continue;
						}
						if (higher > chrLengths[chrNum])
						{
							outs << "File " << fname << ", line " << lineNum << ": Read exceeds end of chromosome, read excluded." << endl;
							continue;
						}
						char strand = tokens[5].c_str()[0];
						if (!useFivePrime)
						{
							bpLoc = (lower + higher) / 2;
						}
						else
						{
							if (strand == '-')
							{
								bpLoc = higher;
							}
							else
							{
								bpLoc = lower;
							}
						}
						if (strand == '-')
						{
							if (readShift > 0)
							{
								if (bpLoc < (unsigned long) readShift) bpLoc = 0;
								else bpLoc -= readShift;
							}
						}
						else if (strand == '+')
						{
							bpLoc += readShift;
						}
						else
						{
#ifdef EXCLUDE_UNKNOWN_STRAND
							outs << "Line " << lineNum << ": Strand " << strand << " unknown, read excluded." << endl;
							continue;
#else
							outs << "Line " << lineNum << ": Strand " << strand << " unknown, assuming +." << endl;
							bpLoc += readShift;
#endif
						}
						if (bpLoc > chrLengths[chrNum]) bpLoc = chrLengths[chrNum];
						bin = bpLoc / interval;
						m[chrNum][bin]++;
					}
				}
			}
			file.close();
		}
		else
		{
			cerr << "Error opening file " << fname << endl;
			exit(EXIT_FAILURE);
		}
		return m;
	}

	unsigned long dataLoad::totalReads(int** reads, int interval)
	{
		unsigned long count = 0;
		for (int i = 0; i < options::getOptions()->chromosomes; i++)
		{
			unsigned long numBins = (options::getOptions()->chrLengths[i] / interval) + 1;
			for (unsigned long j=0; j<numBins; j++)
			{
				count += reads[i][j];
			}
		}
		return count;
	}

	long double dataLoad::getMeanMarks(unsigned long totalReads, int interval)
	{
		unsigned long totalGenomeLength = 0;
		for (int i=0; i < options::getOptions()->chromosomes; i++)
		{
			totalGenomeLength += options::getOptions()->chrLengths[i];
		}
		return (long double) totalReads / ((long double) totalGenomeLength / (long double) interval);
	}

	int dataLoad::getThreshold(long double meanMarks, double thresholdToCallMark)
	{
		bool flag = false;
		int i = (int) ceil(meanMarks);
		double p;
		while (!flag)
		{
			p = poisson(i, meanMarks);
			if (p < thresholdToCallMark) flag = true;
			else i++;
		}
		return i;
	}

	long double dataLoad::getContinuousThreshold(long double meanMarks, double thresholdToCallMark)
	{
		bool flag = false;
		int i = (int) ceil(meanMarks);
		double p;
		while (!flag)
		{
			p = poisson(i, meanMarks);
			if (p < thresholdToCallMark) flag = true;
			else i++;
		}
		long double actualThreshold;
		return continuousThresholdRec(meanMarks, i - 1, i, numThresholdSearchIterations, actualThreshold, thresholdToCallMark);
	}

	long double dataLoad::continuousThresholdRec(long double& meanMarks, long double lower, long double upper, int iterationNum, long double& curThreshold, double& thresholdToCallMark)
	{
		long double curLoc = (upper + lower) / 2.0;
		curThreshold = continuousPoisson(curLoc, meanMarks);
		iterationNum--;
		if (iterationNum > 0)
		{
			long double nextThreshold = thresholdToCallMark * -1.0; //to ensure that on the oft chance that curThreshold == thresholdToCallMark, curLoc is not overwritten
			long double nextLoc = 0;
			if (curThreshold < thresholdToCallMark)
			{
				nextLoc = continuousThresholdRec(meanMarks, lower, curLoc, iterationNum, nextThreshold, thresholdToCallMark);
			}
			else if (curThreshold > thresholdToCallMark)
			{
				nextLoc = continuousThresholdRec(meanMarks, curLoc, upper, iterationNum, nextThreshold, thresholdToCallMark);
			}
			if (fabs(nextThreshold - thresholdToCallMark) < fabs(curThreshold - thresholdToCallMark)) curLoc = nextLoc;
		}
		return curLoc;
	}

	bool** dataLoad::getMarks(int** reads, int threshold, int interval, int* numBinsMarked)
	{
	    if(numBinsMarked != NULL) *numBinsMarked = 0;
		bool** m = new bool*[options::getOptions()->chromosomes];
		for (int i=0; i < options::getOptions()->chromosomes; i++)
		{
			unsigned long numBins = (options::getOptions()->chrLengths[i] / interval) + 1;
			bool* arr = new bool[numBins];
			int* numReads = reads[i];
			for (unsigned long j = 0; j < numBins; j++)
			{
				if (numReads[j] >= threshold)
				{
				    if (numBinsMarked != NULL) (*numBinsMarked)++;
				    arr[j] = true;
				}
				else arr[j] = false;
			}
			m[i] = arr;
		}
		return m;
	}

	bool** dataLoad::getMarks(double** reads, double threshold, int interval, int* numBinsMarked)
	{
	    if (numBinsMarked != NULL) *numBinsMarked = 0;
		bool** m = new bool*[options::getOptions()->chromosomes];
		for (int i=0; i < options::getOptions()->chromosomes; i++)
		{
			unsigned long numBins = (options::getOptions()->chrLengths[i] / interval) + 1;
			bool* arr = new bool[numBins];
//			int* numReads = reads[i];
			for (unsigned long j = 0; j < numBins; j++)
			{
				if (reads[i][j] >= threshold)
				{
				    if (numBinsMarked != NULL) (*numBinsMarked)++;
				    arr[j] = true;
				}
				else arr[j] = false;
			}
			m[i] = arr;
		}
		return m;
	}

	void dataLoad::putMarksInData(bool** marks,  int** data, int bit, int interval)
	{
//#pragma omp parallel for
		for (int i=0; i < options::getOptions()->chromosomes; i++)
		{
			unsigned long numBins = (options::getOptions()->chrLengths[i] / interval) + 1;
			for (unsigned long j = 0; j < numBins; j++)
			{
				if (marks[i][j]) data[i][j] += (1 << bit);
			}
		}
	}

	void dataLoad::loadPrecompData(ifstream& is, int** data)
	{
		string line;
		int chr, bin;
		while (!! getline(is, line))
		{
			if (!line.empty())
			{
				stringstream ss(line);
				ss >> chr;
				ss >> bin;
				ss >> data[chr][bin-1];
			}
		}
	}

	void dataLoad::loadPrecompData_out(ifstream& is, int** data, int** start, int** end)
	{
		string line;
		int chr, bin;
		while (!!getline(is, line))
		{
			if (!line.empty())
			{
				stringstream ss(line);
				ss >> chr;
				ss >> bin;
				ss >> data[chr][bin-1];
				
				ss >> start[chr][bin-1];
				ss >> end[chr][bin-1];
				
			}
		}
	}

	double** dataLoad::readFileDistributeReads(const string& fname, unsigned long& totalReads, ostream& outs, int interval)
	{
		unsigned long* chrLengths = options::getOptions()->chrLengths;
		totalReads = 0;
		double** m = new double*[options::getOptions()->chromosomes];
		for (int i = 0; i < options::getOptions()->chromosomes; i++)
		{
			unsigned long numBins = (chrLengths[i] / interval) + 1;
			double* arr = new double[numBins];
			for (unsigned long j=0; j<numBins; j++)
			{
				arr[j] = 0.0;
			}
			m[i] = arr;
		}
		string chr;
		string junk;
		unsigned long fivePrime;
		unsigned long threePrime;
//		unsigned long bin;
		ifstream file(fname.c_str());
		unsigned long lineNum = 0;
		if (file.is_open())
		{
			string line;
			while(!! getline(file,line))
			{
			    lineNum++;
				if (!line.empty())
				{
					//string buf;
					//stringstream ss(line);
					vector<string> tokens;
					string delim = "\t";
					size_t start = 0;
					size_t end = line.find(delim);
					while (end != string::npos)
					{
						tokens.push_back(line.substr(start, end - start));
						start = end + delim.length();
						end = line.find(delim, start);
					}
					tokens.push_back(line.substr(start, end));
					//while (ss >> buf) tokens.push_back(buf);
					if (tokens.size() >= 6)
					{
						chr = tokens[0];
//						string chrStr = chr.substr(3,chr.length());
						int chrNum = -1; //so it'll be loud and clear with a segfault if there's something wrong here
						if (options::getOptions()->chrNames.count(chr) == 0)
						{
#ifdef WARN_INVALID_CHROMOSOME
							outs << "File " << fname << ", line " << lineNum << ": Invalid chromosome " << chr << endl;
#endif
							continue;
						}
						else
						{
							chrNum = options::getOptions()->chrNames[chr];
						}
						fivePrime = strtoul(tokens[1].c_str(),NULL,0);
						threePrime = strtoul(tokens[2].c_str(),NULL,0);
						if (fivePrime > threePrime) swap(fivePrime, threePrime); //we don't care about the strand
						int lowestBin = fivePrime / interval;
						if (fivePrime < chrLengths[chrNum] && chrNum != 25)
						{
							totalReads++;
							int highestBin = threePrime / interval;
							double lowestBinAddition = ((double)(interval - (fivePrime - (lowestBin * interval)))) / ((double)interval);
							double highestBinAddition = ((double)(threePrime - (highestBin * interval))) / ((double) interval);
							if (lowestBin == highestBin)
							{
								m[chrNum][lowestBin] += (highestBin - lowestBin) / interval;
							}
							else
							{
								m[chrNum][lowestBin] += lowestBinAddition;
								if (threePrime < chrLengths[chrNum])
								{
									for (int i = lowestBin + 1; i < highestBin; i++)
									{
										m[chrNum][i]++;
									}
									m[chrNum][highestBin] += highestBinAddition;
								}
								else
								{
									outs << "Read exceeds length of chromosome." << endl;
									for (unsigned int i = lowestBin + 1; i < chrLengths[chrNum] / interval; i++)
									{
										m[chrNum][i]++;
									}
								}
							}
						}
						else if (chrNum != 25) outs << fivePrime << " " << chrLengths[chrNum] << " Read does not fit chromosome." << endl;
					}
				}
			}
			file.close();
		}
		else
		{
			cerr << "Error opening file " << fname << endl;
			exit(EXIT_FAILURE);
		}
		return m;
	}

	//edited from my patch to BELT
	void dataLoad::getChrlenFromListFile(ifstream& input)
	{
		int numChrom;
		list<string> chrNameList;
		list<unsigned long> chrLengthList;
		string line = "";
		while (getline(input, line))
		{
			stringstream ss(line);
			string chrName;
			unsigned long chrLength;
			ss >> chrName;
			ss >> chrLength;
			chrNameList.push_back(chrName);
			chrLengthList.push_back(chrLength);
		}
		numChrom = chrNameList.size();
		options::getOptions()->chromosomes = numChrom;
		unsigned long* chrLengths = new unsigned long[numChrom];
		for (int i = 0; i < numChrom; i++)
		{
			options::getOptions()->chrNames[chrNameList.front()] = i;
			options::getOptions()->chrIDs[i] = chrNameList.front();
			chrNameList.pop_front();
			chrLengths[i] = chrLengthList.front();
			chrLengthList.pop_front();
		}
		if (options::getOptions()->chrLengths != NULL) delete[] options::getOptions()->chrLengths;
		options::getOptions()->chrLengths = chrLengths;
	}

	//revised
	void dataLoad::getChrBinFromListFile(ifstream& input)
	{
		int numChrom;
		list<string> chrNameList;
		list<unsigned long> chrBinNumList;
		string line = "";
		while (getline(input, line))
		{
			stringstream ss(line);
			string chrName;
			unsigned long chrbin;
			ss >> chrName;
			ss >> chrbin;
			chrNameList.push_back(chrName);
			chrBinNumList.push_back(chrbin);
		}
		numChrom = chrNameList.size();
		options::getOptions()->chromosomes = numChrom;
		unsigned long* chrBins = new unsigned long[numChrom];
		for (int i = 0; i < numChrom; i++)
		{
			chrNameList.pop_front();
			chrBins[i] = chrBinNumList.front();
			chrBinNumList.pop_front();
		}
		if (options::getOptions()->chrBins != NULL) delete[] options::getOptions()->chrBins;
		options::getOptions()->chrBins = chrBins;
	}
} //namespace HMM
