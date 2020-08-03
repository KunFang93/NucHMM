#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <omp.h>

#include "constants.hpp"
#include "dataLoad.hpp"
#include "HMMIO.hpp"
#include "macros.hpp"
#include "masterUtils.hpp"
#include "options.hpp"
#include "stateLD.hpp"
#include "viterbi.hpp"
#include "writeResults.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::ifstream;
using std::string;

using namespace HMM;

int main(int argc, char** argv)
{
    if (argc == 1)
    {
        cout << "Usage NucHMM-output_result customGenome binNumFile states.bed outputs.bed data.precomp HMM.rawhmm" << endl;
        cout << "Enter \"NONE\" for customGenome if using the default (hg19)" << endl;
        return EXIT_FAILURE;
    }
    options::setHg19ChrLengthArray();
    string customGenomeFname(argv[1]);
	string binnumFname(argv[2]);
    string stateFname(argv[3]);
    string outputFname(argv[4]);
    string dataFname(argv[5]);
    string HMMFname(argv[6]);
    if (customGenomeFname == "NONE") options::setHg19ChrLengthArray();
    else
    {
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
	if (binnumFname != "NONE")
	{
		ifstream binFile(binnumFname.c_str());
		if (binFile.is_open())
		{
			dataLoad::getChrBinFromListFile(binFile);
			binFile.close();
		}
		else
		{
			cerr << "Error: cannot open " << binnumFname << endl;
			return EXIT_FAILURE;
		}
	}
    int chromosomes = options::getOptions()->chromosomes;
    unsigned long* chrLengths = options::getOptions()->chrLengths;
	unsigned long* chrBins = options::getOptions()->chrBins;
	int** data = dataLoad::initDataArray_n<int>(chrBins);
	int** start = dataLoad::initDataArray_n<int>(chrBins);
	int** end = dataLoad::initDataArray_n<int>(chrBins);

	cout << "Loading data.....";
	ifstream dataFile(dataFname.c_str());
	dataLoad::loadPrecompData_out(dataFile, data, start, end);
	dataFile.close();
	cout << "Done." << endl;

	State** states;
	int numStates;
	int numOutputs;
	long double* startProbs;

	cout << "Loading HMM.....";
	ifstream HMMfile(HMMFname.c_str());
	HMMIO::inputHMM(HMMfile, states, numStates, numOutputs, startProbs);
	HMMfile.close();
	cout << "Done." << endl;

	int** calledStates = new int*[chromosomes];
	cout << "Calling states....." << flush;
#pragma omp parallel for
	for (int i = 0; i < chromosomes; i++)
	{
		unsigned long intervals = chrBins[i];
		calledStates[i] = viterbiAlg::viterbi<long double>(numStates, numOutputs, states, startProbs, data[i], intervals);
	}
	cout << "Done." << endl;

	cout << "Writing results....." << flush;
	resultsWriter::writeResults_n(calledStates, data, stateFname, outputFname, chrBins, start, end);
	cout << "Done." << endl;

    options::cleanupOptions();
	return EXIT_SUCCESS;
}
