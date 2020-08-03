#include "master.hpp"

#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <sys/time.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "constants.hpp"
#include "dataLoad.hpp"
#include "getPrH.hpp"
#include "HMMIO.hpp"
#include "macros.hpp"
#include "masterUtils.hpp"
#include "options.hpp"
#include "processChr.hpp"
#include "reestimator.hpp"
#include "stateLD.hpp"

using std::cerr;
using std::endl;
using std::flush;
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
using boost::mt19937;
using boost::uniform_real;
using boost::variate_generator;

namespace HMM
{
    master::master()
    {
        numStates = options::getOptions()->numStates;
        numOutputs = options::getOptions()->numOutputs;
    }

    int master::run(ostream& outs)
    {
        options* opts = options::getOptions();
        if (opts->customGenome)
        {
            ifstream genomeFile(opts->customGenomeFile.c_str());
            if (genomeFile.is_open())
            {
                dataLoad::getChrlenFromListFile(genomeFile);
                genomeFile.close();
            }
            else
            {
                cerr << "Error: cannot open " << opts->customGenomeFile << endl;
                return EXIT_FAILURE;
            }
        }
        else options::setHg19ChrLengthArray();
		if (opts->numbin)
		{
			ifstream binFile(opts->numbinFile.c_str());
			if (binFile.is_open())
			{
				dataLoad::getChrBinFromListFile(binFile);
				binFile.close();
			}
			else
			{
				cerr << "Error: cannot open " << opts->numbinFile << endl;
				return EXIT_FAILURE;
			}
		}
		else
		{
			cerr << "Error: Num of Bins File Needed" << endl;
			return EXIT_FAILURE;
		};
        int chromosomes = options::getOptions()->chromosomes;
        unsigned long* chrLengths = opts->chrLengths;
		unsigned long* chrBins = opts->chrBins;
		cout << chrBins[0] << endl;
        int numDataSets = opts->numDataSets;
        const int reps = opts->iterations;
        const int numSequences = chromosomes * numDataSets;
        long double minProbToHandle = opts->minProb;

        string* initialHMMFilename;
        if (opts->usingCustomInitialHMM) initialHMMFilename = &opts->initialHMMFile;
        else initialHMMFilename = NULL;
        string* finalHMMFilename;
        if (opts->outputHMMFile) finalHMMFilename = &opts->finalHMMFile;
        else finalHMMFilename = NULL;

        int** data = dataLoad::initMultipleDataSetDataArray_n<int>(numDataSets, chrBins);
        {
			bool error = false;
//#ifdef _OPENMP
//			omp_set_num_threads(30);
//#endif
#pragma omp parallel for shared(error)
            for (int i = 0; i < numDataSets; i++)
            {
                int** curLoc = data + (i * chromosomes);
                ifstream input(opts->fnames[i]); //three arguments before these
                if (input.is_open())
                {
                    dataLoad::loadPrecompData(input, curLoc);
                    input.close();
                }
                else
                {
#pragma omp critical
                    {
                        cerr << "Error: cannot open " << opts->fnames[i] << endl;
                        error = true;
                    }
                }
            }
            if (error) return EXIT_FAILURE;
        }

        unsigned long totalIntervals = 0;
        unsigned long* intervals = new unsigned long[numSequences];
        for (int i = 0; i < numSequences; i++)
        {
            totalIntervals += chrLengths[i % chromosomes];
			intervals[i] = chrBins[i % chromosomes];            
        }
        timeval curTime;
        gettimeofday(&curTime, NULL);
        mt19937 gen( static_cast<unsigned int>(curTime.tv_usec));
        variate_generator<mt19937, uniform_real<> > rand(gen, uniform_real<double>(0.0, 1.0));

        long double* startProb;
        State** HMM;
        if (initialHMMFilename != NULL)
        {
            ifstream HMMInput(initialHMMFilename->c_str());
            if (HMMInput.is_open())
            {
                HMMIO::inputHMM(HMMInput, HMM, numStates, numOutputs, startProb);
                HMMInput.close();
            }
            else
            {
                cerr << "Error: cannot open " << initialHMMFilename << endl;
                return EXIT_FAILURE;
            }
        }
        else masterUtils::generateRandomHMM(numStates, numOutputs, rand, HMM, startProb);
        long double*** gammaRowSums = new long double**[numSequences];
        long double*** deltas = new long double**[numSequences];
        // ifNdef
        for (int seq = 0; seq < numSequences; seq++)
        {
            gammaRowSums[seq] = new long double*[numStates];
            for (int i = 0; i < numStates; i++)
            {
                gammaRowSums[seq][i] = new long double[numStates];
            }
        }
        for (int i = 0; i < numSequences; i++)
        {
            deltas[i] = new long double*[numStates];
            for (int j = 0; j < numStates; j++)
            {
                deltas[i][j] = new long double[intervals[i]];
            }
        }
        long double lastBIC = 0l;
        {
            long double logLikelihood = 0.0l;
//#ifdef _OPENMP
//			omp_set_num_threads(30);
//#endif
#pragma omp parallel for reduction (+ : logLikelihood)
            for (int i = 0; i < numSequences; i++)
            {
                logLikelihood += prHCalculator::getPrH<long double>(numStates, HMM, startProb, data[i], intervals[i]).getLog();
            }
            lastBIC = masterUtils::calcBIC(numStates, numOutputs, logLikelihood, totalIntervals);
            outs << "Log likelihood: " << logLikelihood << endl;
            outs << "Initial BIC: " << lastBIC << endl;
            outs << endl;
        }
        for (int iteration = 0; iteration < reps; iteration += (reps != -1 ? 1 : 0))
        {
#pragma omp parallel for
            for (int i = 0; i < numSequences; i++)
            {
                processChr::baumWelch_canonical_scaled<long double>(numOutputs, numStates, HMM, startProb, data[i], intervals[i], gammaRowSums[i], deltas[i]);
            }

            //re-estimation
            masterUtils::deleteHMMStates(numStates, HMM);
            HMM = reestimator::getNextParams_threaded(numStates, numOutputs, numSequences, deltas, gammaRowSums, data, intervals, startProb);
            //delete all the partially trained HMMs
            if (opts->zeroProbs)
            {
                //check if any probs are below thresholdToZeroProb, if so, set to 0 and redistribute
                //the 2 makes sure that no states can be expelled from the model (2 to prevent two states from only transitioning to each other)
                masterUtils::zeroProbs(HMM, numStates, numStates, &State::getTransProb, &State::setTransProb, minProbToHandle, 2);
                //the 1 makes sure that no outputs are considered to be impossible to occur
                masterUtils::zeroProbs(HMM, numStates, numOutputs, &State::getEmitProb, &State::setEmitProb, minProbToHandle, 1);
            }
            else if (opts->enforceMinProb)
            {
                masterUtils::forceMinimumProbs(HMM, numStates, numOutputs, startProb, minProbToHandle);
            }
            if (opts->printAllProbs)
            {
                outs << "Iteration " << iteration + 1 << ":" << endl;
                HMMIO::printHMM(outs, HMM, numStates, numOutputs, startProb);
            }
	    else
	    {
		outs << "Iteration " << iteration + 1 << endl;
	    }
            if (opts->bicEveryIteration || iteration + 1 == reps)
            {
                if (!opts->printAllProbs)
                {
                    outs << "Iteration " << iteration + 1 << ":" << endl;
                }
                long double logLikelihood = 0.0l;
#pragma omp parallel for reduction (+ : logLikelihood)
                for (int i = 0; i < numSequences; i++)
                {
                    logLikelihood += prHCalculator::getPrH<long double>(numStates, HMM, startProb, data[i], intervals[i]).getLog();
                }
                long double curBIC = masterUtils::calcBIC(numStates, numOutputs, logLikelihood, totalIntervals);
                outs << "Log likelihood: " << logLikelihood << endl;
                outs << "BIC: " << curBIC << endl;
                outs << endl;
                if (lastBIC - curBIC < 0)
                {
                    //this would indicate a serious bug (or abuse of the -z option or a nonsensical input HMM),
                    //but that'd be better than invalid results and a retraction years later
                    cerr << "ERROR: BIC increase detected." << endl;
                }
                else if (lastBIC - curBIC < options::getOptions()->minBICChange)
                {
                    //change in BIC below minimum, stop training
                    iteration = reps;
                }
                lastBIC = curBIC;
            }
        }
        if (finalHMMFilename != NULL)
        {
            ofstream HMMOutput(finalHMMFilename->c_str());
            HMMIO::outputHMM(HMMOutput, HMM, numStates, numOutputs, startProb);
            HMMOutput.close();
        }

        //free memory
        masterUtils::deleteHMMStates(numStates, HMM);
        delete[] startProb;
        delete[] intervals;
        for (int i = 0; i < numSequences; i++)
        {
            for (int j = 0; j < numStates; j++)
            {
                delete[] deltas[i][j];
            }
            delete[] deltas[i];
        }
        delete[] deltas;
        for (int i = 0; i < numSequences; i++)
        {
            for (int j = 0; j < numStates; j++)
            {
                delete[] gammaRowSums[i][j];
            }
            delete[] gammaRowSums[i];
        }
        delete[] gammaRowSums;
        for (int i = 0; i < numSequences; i++)
        {
            delete[] data[i];
        }
        delete[] data;
        return EXIT_SUCCESS;
    }
}
