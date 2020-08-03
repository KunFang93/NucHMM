#include "test.hpp"

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
//#include <set>
#include <vector>

#include <sys/time.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "dataLoad.hpp"
#include "getPrH.hpp"
#include "HMMIO.hpp"
#include "macros.hpp"
#include "masterUtils.hpp"
#include "processChr.hpp"
#include "reestimator.hpp"
#include "stateLD.hpp"

using std::cerr;
using std::endl;
using boost::mt19937;
using boost::uniform_real;
using boost::variate_generator;

using boost::mt19937;
using boost::uniform_real;
using boost::variate_generator;
using std::cerr;
using std::cout;
using std::endl;
using std::floor;
using std::pow;
//using std::set;
using std::sqrt;
using std::string;
using std::time;
using std::vector;

//zero probs if they drop below thresholdToZeroProb
//#define ZERO_PROBS
//enforce a minimum prob of minProbToForce, if you use both this and ZERO_PROBS you're just silly :p
#define FORCE_MINIMUM_PROB

namespace HMM
{
    testLearn::testLearn()
    {}

    //a test driver for the HMM learning code that does not use MPI
    int testLearn::test(int iterations)
    {
        cout << "Poisson function tests:" << endl;
        cout << dataLoad::poisson(1, 3) << endl;
        cout << dataLoad::poisson(2, 3) << endl;
        cout << dataLoad::continuousPoisson(1, 3) << endl;
        cout << dataLoad::continuousPoisson(1.5, 3) << endl;
        cout << dataLoad::continuousPoisson(2, 3) << endl;

        int numSequences = 20;
        int numStates = 5;
        int numOutputs = 6;
        int testSeqLength = 5000;
        long double minProbToHandle = 1.0E-6l;
        int** data = new int*[numSequences];
        timeval curTime;
        gettimeofday(&curTime, NULL);
        mt19937 gen( static_cast<unsigned int>(curTime.tv_usec));
        variate_generator<mt19937, uniform_real<> > rand(gen, uniform_real<double>(0.0, 1.0));
        std::srand(time(NULL));
        unsigned long* intervals = new unsigned long[numSequences];
        unsigned long totalIntervals = 0;
        for (int i = 0; i < numSequences; i++)
        {
            intervals[i] = testSeqLength + 1000 * i;
//            cout << intervals[i] << endl;
            totalIntervals += intervals[i];
        }
        for (int seq = 0; seq < numSequences; seq++)
        {
            data[seq] = new int[intervals[seq]];
            for (unsigned int i = 0; i < intervals[seq]; i++)
            {
                //add an artificial bias so there's actually something to learn besides a flat distribution
                if (std::rand() % 3 == 0) data[seq][i] = 2;
                else data[seq][i] = (std::rand() % numOutputs);
            }
        }

        long double* startProb;
        State** HMM;
        masterUtils::generateRandomHMM(numStates, numOutputs, rand, HMM, startProb);
        long double*** gammaRowSums = new long double**[numSequences];
        long double*** deltas = new long double**[numSequences];
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
        long double logLikelihood = 0;
#pragma omp parallel for
        for (int i = 0; i < numSequences; i++)
        {
            //likelihood under the PRIOR HMM
            logLikelihood += prHCalculator::getPrH<long double>(numStates, HMM, startProb, data[i], intervals[i]).getLog();
        }
        cout << "Initial BIC: " << masterUtils::calcBIC(numStates, numOutputs, logLikelihood, totalIntervals) << endl;
        cout << endl;
        for (int iteration = 0; iteration < iterations; iteration++)
        {
            int seqsRemaining = numSequences;
#pragma omp parallel shared(seqsRemaining)
            {
#pragma omp for
                for (int i = 0; i < numSequences; i++)
                {
                    processChr::baumWelch_canonical_scaled<long double>(numOutputs, numStates, HMM, startProb, data[i], intervals[i], gammaRowSums[i], deltas[i]);
#pragma omp atomic
                    seqsRemaining--;
                }
            }

            //re-estimation
            masterUtils::deleteHMMStates(numStates, HMM);
            HMM = reestimator::getNextParams_threaded(numStates, numOutputs, numSequences, deltas, gammaRowSums, data, intervals, startProb);
#if defined ZERO_PROBS
            //check if any probs are below thresholdToZeroProb, if so, set to 0 and redistribute
            //the 2 makes sure that no states can be expelled from the model (2 to prevent two states from only transitioning to each other)
            masterUtils::zeroProbs(HMM, numStates, numStates, &State::getTransProb, &State::setTransProb, minProbToHandle, 2);
            //the 1 makes sure that no outputs are considered to be impossible to occur
            masterUtils::zeroProbs(HMM, numStates, numOutputs, &State::getEmitProb, &State::setEmitProb, minProbToHandle, 1);
#elif defined FORCE_MINIMUM_PROB
            masterUtils::forceMinimumProbs(HMM, numStates, numOutputs, startProb, minProbToHandle);
#endif
            cout << "Iteration " << iteration + 1 << ":" << endl;
            HMMIO::printHMM(cout, HMM, numStates, numOutputs, startProb, 1);
            long double logLikelihood = 0;
#pragma omp parallel for
            for (int i = 0; i < numSequences; i++)
            {
                //likelihood under the PRIOR HMM
                logLikelihood += prHCalculator::getPrH<long double>(numStates, HMM, startProb, data[i], intervals[i]).getLog();
            }
            cout << "BIC: " << masterUtils::calcBIC(numStates, numOutputs, logLikelihood, totalIntervals) << endl;
            cout << endl;
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
