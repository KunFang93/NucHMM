#include "masterUtils.hpp"

#include <cmath>
#include <iostream>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include "constants.hpp"
#include "lnsnum.hpp"
#include "macros.hpp"
#include "options.hpp"
#ifdef FILTER_RUNS
#include "repetitivePattern.hpp"
#endif
#include "stateLD.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

using boost::mt19937;
using boost::uniform_real;
using boost::variate_generator;

using std::fabs;
using std::ostream;
using std::pow;
using std::sqrt;

namespace HMM
{
	long double _pure masterUtils::calcBIC(int numStates, int numOutputs, Lnsnum* prH, unsigned long totalIntervals)
	{
		Lnsnum corpusLikelihood = 1;
		for (int i = 0; i < options::getOptions()->chromosomes; i++)
		{
			corpusLikelihood *= prH[i];
		}
		return (long double) -2.0 * corpusLikelihood.getLog() + ((((long double) numStates) * (((long double) numStates - 1.0) + ((long double) numOutputs - 1.0)) + ((long double) numStates)) * log((long double) totalIntervals));
	}

	long double _pure masterUtils::calcBIC(int numStates, int numOutputs, long double logLikelihood, unsigned long totalIntervals)
	{
		return (long double) -2.0 * logLikelihood + ((((long double) numStates) * (((long double) numStates - 1.0) + ((long double) numOutputs - 1.0)) + ((long double) numStates)) * log((long double) totalIntervals));
	}

	void masterUtils::zeroProbs(State** states, int numStates, int numDataEntries, long double (State::*getDatum) (int) const, void (State::*setDatum) (int, long double), long double minProbToHandle, int numConnectionsToPreserve)
	{
#pragma omp parallel for
		for (int state = 0; state < numStates; state++)
		{
			State* curState = states[state];
			int numZeros = 0;
			int numLowerThanThreshold = 0;
			for (int j = 0; j < numDataEntries; j++)
			{
				//this may be a bug due to float precision, I'll fix it later if it is
				if ((curState->*getDatum)(j) == 0.0) numZeros++;
				else if ((curState->*getDatum)(j) < minProbToHandle) numLowerThanThreshold++;
			}
			if (numLowerThanThreshold > 0)
			{
				for (int j = 0; j < numDataEntries; j++)
				{
					if ((curState->*getDatum)(j) > 0.0 &&
						(curState->*getDatum)(j) < minProbToHandle)
					{
						//make sure we don't excise states or observed outputs
						int connections = 0;
						for (int refState = 0; (connections <= numConnectionsToPreserve) && (refState < numStates); refState++)
						{
							if (((states[refState])->*getDatum)(j) >= minProbToHandle) connections++;
						}
						if (connections > numConnectionsToPreserve)
						{
							long double probToTransfer = (curState->*getDatum)(j) / ((long double) (numZeros + numLowerThanThreshold));
							(curState->*setDatum)(j, 0.0);
							for (int k = 0; k < numDataEntries; k++)
							{
								long double datum = (curState->*getDatum)(k);
								if (datum >= minProbToHandle)
								{
									(curState->*setDatum)(k, datum + probToTransfer);
								}
							}
						}
					}
				}
			}
		}
	}

	void masterUtils::forceMinimumProbs(State** states, int numStates, int numOutputs, long double* startProbs, long double minProbToHandle)
	{
		#pragma omp parallel for
		for (int state = 0; state < numStates; state++)
		{
			State* curState = states[state];
			for (int i = 0; i < numStates; i++)
			{
				long double prob = curState->getTransProb(i);
				curState->setTransProb(i, minProbToHandle + (1 - (numStates * minProbToHandle)) * prob);
			}
			for (int i = 0; i < numOutputs; i++)
			{
				long double prob = curState->getEmitProb(i);
				curState->setEmitProb(i, minProbToHandle + (1 - (numOutputs * minProbToHandle)) * prob);
			}
		}
		for (int state = 0; state < numStates; state++)
		{
			long double prob = startProbs[state];
			startProbs[state] = minProbToHandle + (1 - (numStates * minProbToHandle)) * prob;
		}
	}

	void masterUtils::generateRandomHMM(int numStates, int numOutputs, variate_generator<mt19937, uniform_real<> >& rand, State**& states, long double*& startProb)
	{
		states = new State*[numStates];
		for (int i = 0; i < numStates; i++)
		{
			states[i] = new State(numStates, numOutputs);
			states[i]->initRandom();
		}
		startProb = new long double[numStates];
		long double count = 0;
		//avoids a potential bug: what if the variate_generator passed to this was initialized at the same second as the static variate_generator held by HMM::State?
		for (int j = 0; j < (numStates * numOutputs) + (numStates * numStates); j++)
		{
			rand();
		}
		for (int j = 0; j < numStates; j++)
		{
			startProb[j] = rand();
			count += startProb[j];
		}
		for (int j = 0; j < numStates; j++)
		{
			startProb[j] /= count;
		}
	}

	void masterUtils::deleteHMMStates(int numStates, State**& HMM)
	{
		for (int j = 0; j < numStates; j++)
		{
			delete HMM[j];
		}
		delete[] HMM;
	}

} //namespace HMM
