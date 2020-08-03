#include "writeResults.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "constants.hpp"
#include "macros.hpp"
#include "options.hpp"

#include "viterbi.hpp"

using std::cerr;
using std::endl;
using std::flush;
using std::ofstream;
using std::string;
using std::stringstream;

namespace HMM
{

	void resultsWriter::writeResults(int** stateMap, int** outputMap, string stateFname, string outputFname, int interval)
	{
        unsigned long* chrLengths = options::getOptions()->chrLengths;
#ifdef HMM_WRITE_OUTPUT_MAP
		ofstream outputFile;
		outputFile.open(outputFname.c_str());
#endif
#ifdef HMM_WRITE_STATE_MAP
		ofstream stateFile;
		stateFile.open(stateFname.c_str());
#endif
		if (
#ifdef HMM_WRITE_OUTPUT_MAP
			outputFile.is_open() &&
#endif
#ifdef HMM_WRITE_STATE_MAP
			stateFile.is_open() &&
#endif
			true)
		{
			for (int i = 0; i < options::getOptions()->chromosomes; i++)
			{
				string chrName = options::getOptions()->chrIDs[i];
				unsigned long intervals = (chrLengths[i] / interval) + 1;
				for (unsigned long j = 0; j < intervals; j++)
				{
						unsigned long startChr = j * interval; //start at 0, BED is half-open
						unsigned long endChr = (j + 1) * interval;
						if (endChr > chrLengths[i]) endChr = chrLengths[i]; //to keep UCSC happy
#ifdef HMM_WRITE_STATE_MAP
						stateFile << chrName << "\t" << startChr << "\t" << endChr << "\tHMM_MPI\t" << stateMap[i][j] << "\t." << endl;
#endif
#ifdef HMM_WRITE_OUTPUT_MAP
						outputFile << chrName << "\t" << startChr << "\t" << endChr << "\tHMM_MPI\t" << outputMap[i][j] << "\t." << endl;
#endif
				}
			}
#ifdef HMM_WRITE_OUTPUT_MAP
			outputFile.close();
#endif
#ifdef HMM_WRITE_STATE_MAP
			stateFile.close();
#endif
		}
		else
		{
			cerr << "Error opening an output file." << endl;
		}
	}

	void resultsWriter::writeResults_n(int** stateMap, int** outputMap, string stateFname, string outputFname, unsigned long* numbins, int** start, int** end)
	{
		unsigned long* chrLengths = options::getOptions()->chrLengths;
#ifdef HMM_WRITE_OUTPUT_MAP
		ofstream outputFile;
		outputFile.open(outputFname.c_str());
#endif
#ifdef HMM_WRITE_STATE_MAP
		ofstream stateFile;
		stateFile.open(stateFname.c_str());
#endif
		if (
#ifdef HMM_WRITE_OUTPUT_MAP
			outputFile.is_open() &&
#endif
#ifdef HMM_WRITE_STATE_MAP
			stateFile.is_open() &&
#endif
			true)
		{
			for (int i = 0; i < options::getOptions()->chromosomes; i++)
			{
				string chrName = options::getOptions()->chrIDs[i];
				unsigned long intervals = numbins[i];
				for (unsigned long j = 0; j < intervals; j++)
				{
					if(end[i][j] == 0)
					{
						continue;
					}
					else
					{
#ifdef HMM_WRITE_STATE_MAP
						stateFile << chrName << "\t" << start[i][j] << "\t" << end[i][j] << "\tHMM_MPI\t" << stateMap[i][j]+1 << "\t." << endl;
#endif
#ifdef HMM_WRITE_OUTPUT_MAP
						outputFile << chrName << "\t" << start[i][j] << "\t" << end[i][j] << "\tHMM_MPI\t" << outputMap[i][j] << "\t." << endl;
#endif			
					}
				}
			}
#ifdef HMM_WRITE_OUTPUT_MAP
			outputFile.close();
#endif
#ifdef HMM_WRITE_STATE_MAP
			stateFile.close();
#endif
		}
		else
		{
			cerr << "Error opening an output file." << endl;
		}
	}

    void resultsWriter::writePrecompData(ofstream& os, int** data, int interval)
    {
        for (int i = 0; i < options::getOptions()->chromosomes; i++)
        {
            unsigned long intervals = (options::getOptions()->chrLengths[i] / interval) + 1;
            for (unsigned long j = 0; j < intervals; j++)
            {
                os << i << "\t" << j << "\t" << data[i][j] << endl;
            }
        }
    }

} //namespace HMM
