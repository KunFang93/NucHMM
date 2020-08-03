#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <omp.h>

#include "constants.hpp"
#include "dataLoad.hpp"
#include "HMMIO.hpp"
#include "macros.hpp"
#include "options.hpp"
#include "stateLD.hpp"
#include "viterbi.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::ifstream;
//using std::log2;
using std::string;

using namespace HMM;

#define PER_MARK

int main(int argc, char* argv[])
{
#ifdef PER_MARK
    cout << "Output per mark." << endl << flush;
#endif

#ifdef PER_OBSERVATION
    cout << "Output per observation." << endl << flush;
#endif
    if (argc == 1)
    {
        cerr << "Usage: HMM_MPI_MARK_STATS [custom genome file] [binWidth] [precomp data file] [HMM file]" << endl;
	cerr << "Enter \"NONE\" for custom genome file for default (hg18)" << endl;
        return EXIT_FAILURE;
    }
    int interval = atoi(argv[2]);
	string customGenomeFname(argv[1]);
	if (customGenomeFname == "NONE") options::setHg18ChrLengthArray();
	else
	{
		cout << "Loading custom genome." << endl << flush;
		ifstream genomeFile(customGenomeFname.c_str());
		if (genomeFile.is_open())
		{
			dataLoad::getChrlenFromListFile(genomeFile);
			genomeFile.close();
		}
		else
		{
			cerr << "Error: cannot open " << customGenomeFname << endl;
			return EXIT_FAILURE;
		}
	}

	cout << "Loading raw data." << endl << flush;
	string inputFilenames[2];
	inputFilenames[0] = argv[3];
	inputFilenames[1] = argv[4];
	//a << x = a * 2^x, but be careful of overflow
	int** data = dataLoad::initDataArray<int>(interval);
	int** calledStates = new int*[options::getOptions()->chromosomes];
	ifstream dataFile(inputFilenames[0].c_str());
	dataLoad::loadPrecompData(dataFile, data);
	dataFile.close();


	cout << "Loading HMM." << endl << flush;
	State** states;
//	long double** transProbs;
//	long double** emitProbs;
	long double* startProbs;
	int numStates = 0;
	int numOutputs = 0;
	ifstream file(inputFilenames[1].c_str());
	HMMIO::inputHMM(file, states, numStates, numOutputs, startProbs);
	int files = log2(numOutputs);
	file.close();
//	std::ofstream outFile("testHMM.rawhmm");
//	outputHMM(outFile, states, numStates, numOutputs, startProbs);
//	exit(1);

	cout << "Calling states." << endl << flush;
	for (int i = 0; i < options::getOptions()->chromosomes; i++)
	{
		unsigned long intervals = (options::getOptions()->chrLengths[i] / interval) + 1;
		calledStates[i] = viterbiAlg::viterbi<long double>(numStates, numOutputs, states, startProbs, data[i], intervals);
	}

	long double** expOutputFreq;
	expOutputFreq = new long double*[numStates];
	for (int i = 0; i < numStates; i++)
	{
#ifdef PER_MARK
		expOutputFreq[i] = new long double[files];
		for (int k = 0; k < files; k++)
#endif
#ifdef PER_OBSERVATION
        expOutputFreq[i] = new long double[numOutputs];
        for (int k = 0; k < numOutputs; k++)
#endif
		{
#ifdef PER_MARK
			expOutputFreq[i][k] = 0.0l;
			for (int j = (1 << k); j < numOutputs; j++)
			{
				if (((1 << k) & j) > 0)
				{
					expOutputFreq[i][k] += states[i]->getEmitProb(j);
				}
			}
#endif
#ifdef PER_OBSERVATION
				expOutputFreq[i][k] = states[i]->getEmitProb(k);
#endif
		}
	}

	//observed output frequency, per state per output
	//they're long doubles since they will be divided later on
	long double** obsOutputFreq;
	obsOutputFreq = new long double*[numStates];
	for (int i = 0; i < numStates; i++)
	{
#ifdef PER_MARK
		obsOutputFreq[i] = new long double[files];
		for (int j = 0; j < files; j++)
#endif
#ifdef PER_OBSERVATION
		obsOutputFreq[i] = new long double[numOutputs];
		for (int j = 0; j < numOutputs; j++)
#endif
		{
			obsOutputFreq[i][j] = 0.0l;
		}
	}

#ifdef PER_MARK
	unsigned long* totalOutputs = new unsigned long[files];
#endif
#ifdef PER_OBSERVATION

	unsigned long* totalOutputs = new unsigned long[numOutputs];
#endif
	unsigned long* stateFrequency = new unsigned long[numStates];
#ifdef PER_MARK
	for (int j = 0; j < files; j++)
#endif
#ifdef PER_OBSERVATION
	for (int j = 0; j < numOutputs; j++)
#endif
	{
		totalOutputs[j] = 0ul;
	}
	for (int j = 0; j < numStates; j++)
	{
		stateFrequency[j] = 0ul;
	}
#pragma omp parallel for
	for (int chrID = 0; chrID < options::getOptions()->chromosomes; chrID++)
	{
		//these avoid atomic/critical blocks in the chromosome location loop
#ifdef PER_MARK
		unsigned long* totalOutputInChr = new unsigned long[files];
#endif
#ifdef PER_OBSERVATION
		unsigned long* totalOutputInChr = new unsigned long[numOutputs];
#endif
		unsigned long* stateFrequencyInChr = new unsigned long[numStates];
#ifdef PER_MARK
		for (int j = 0; j < files; j++)
#endif
#ifdef PER_OBSERVATION
		for (int j = 0; j < numOutputs; j++)
#endif
		{
			totalOutputInChr[j] = 0ul;
		}
		for (int j = 0; j < numStates; j++)
		{
			stateFrequencyInChr[j] = 0ul;
		}
		for (unsigned long loc = 0; loc < options::getOptions()->chrLengths[chrID] / interval; loc++)
		{
#ifdef PER_MARK
			for (int k = 0; k < files; k++)
			{
				if (((1 << k) & data[chrID][loc]) > 0)
				{
					totalOutputInChr[k]++;
#pragma omp critical
					obsOutputFreq[calledStates[chrID][loc]][k]++;
				}
			}
#endif
#ifdef PER_OBSERVATION
			totalOutputInChr[data[chrID][loc]]++;
#pragma omp critical
            obsOutputFreq[calledStates[chrID][loc]][data[chrID][loc]]++;
#endif
			stateFrequencyInChr[calledStates[chrID][loc]]++;
		}
#ifdef PER_MARK
#pragma omp critical
		for (int j = 0; j < files; j++) //./HMM_MARK_STATS MCF7_ERalpha_E2.bed MCF7_ERalpha.bed MCF7_PolII.bed MCF7_h3k4me2.bed MCF7_h3k27me3.bed MCF7_h3k9me2.bed MCF7_methylation.bed HMM3-26-11.rawhmm | tee -a 4-18-11.out
#endif
#ifdef PER_OBSERVATION
#pragma omp critical
		for (int j = 0; j < numOutputs; j++)
#endif
		{
			totalOutputs[j] += totalOutputInChr[j];
		}
#pragma omp critical
		for (int j = 0; j < numStates; j++)
		{
			stateFrequency[j] += stateFrequencyInChr[j];
		}
		delete[] totalOutputInChr;
		delete[] stateFrequencyInChr;
	}

//cout << "State frequencies:" << endl;
	for (int i = 0; i < numStates; i++)
	{
#ifdef PER_MARK
		for (int j = 0; j < files; j++)
#endif

#ifdef PER_OBSERVATION
		cout << stateFrequency[i] << endl;
		for (int j = 0; j < numOutputs; j++)
#endif
		{
			obsOutputFreq[i][j] /= ((long double) stateFrequency[i]);
		}
	}

	for (int i = 0; i < numStates; i++)
	{
//		cout << "State " << i << ":" << endl;
#ifdef PER_MARK
		for (int j = 0; j < files; j++)
		{
//			cout << "Exp:\t" << expOutputFreq[i][j] << "\tObs:\t" << obsOutputFreq[i][j] << endl;
            cout << expOutputFreq[i][j] << "\t" << obsOutputFreq[i][j] << endl;
		}
#endif
#ifdef PER_OBSERVATION
        for (int j = 0; j < numOutputs; j++)
        {
            cout << expOutputFreq[i][j] << "\t" << obsOutputFreq[i][j] << endl;
        }
#endif
	}

//	for (int i = 0; i < numOutputs; i++)
//	{
//		cout << "Total amount of mark " << i << ": " << totalOutputs[i] << endl;
//	}

	for (int chrID = 0; chrID < options::getOptions()->chromosomes; chrID++)
	{
		delete[] data[chrID];
	}
	delete[] data;
	for (int i = 0; i < numStates; i++)
	{
		delete[] obsOutputFreq[i];
		delete states[i];
	}
	delete[] obsOutputFreq;
	delete[] states;
	delete[] startProbs;
	delete[] totalOutputs;
	delete[] stateFrequency;
	options::cleanupOptions();
	return EXIT_SUCCESS;
}
