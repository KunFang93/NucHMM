#include <iostream>

#include "constants.hpp"
#include "macros.hpp"
#include "stateLD.hpp"

#include "HMMIO.hpp"

using std::endl;
using std::flush;
using std::istream;
using std::ostream;

namespace HMM
{
	void HMMIO::outputHMM(ostream& os, State** states, int numStates, int numOutputs, long double* startProbs)
	{
		os << numStates << " " << numOutputs << endl;
		for (int i = 0; i < numStates; i++)
		{
			for (int j = 0; j < numStates; j++)
			{
				os << states[i]->getTransProb(j) << " ";
			}
			os << endl;
		}
		for (int i = 0; i < numStates; i++)
		{
			for (int j = 0; j < numOutputs; j++)
			{
				os << states[i]->getEmitProb(j) << " ";
			}
			os << endl;
		}
		for (int i = 0; i < numStates; i++)
		{
			os << startProbs[i] << " ";
		}
		os << endl;
	}

	void HMMIO::inputHMM(istream& is, State**& states, int& numStates, int& numOutputs, long double*& startProbs)
	{
		numStates = 0;
		numOutputs = 0;
		is >> numStates;
		is >> numOutputs;
		states = new State*[numStates];
		for (int i = 0; i < numStates; i++)
		{
			states[i] = new State(numStates, numOutputs);
		}
		for (int i = 0; i < numStates; i++)
		{
			for (int j = 0; j < numStates; j++)
			{
				long double field;
				is >> field;
				states[i]->setTransProb(j, field);
			}
		}
		for (int i = 0; i < numStates; i++)
		{
			for (int j = 0; j < numOutputs; j++)
			{
				long double field;
				is >> field;
				states[i]->setEmitProb(j, field);
			}
		}
		startProbs = new long double[numStates];
		for (int i = 0; i < numStates; i++)
		{
			is >> startProbs[i];
		}
	}

#define HMMIO_PRINT_HMM_BODY \
		for (int i = 0; i < numStates; i++) \
		{ \
			for (int j = 0; j < numStates; j++) \
			{ \
				os << "Transition probability of state " << i \
				<< " to state " << j << ": " \
				<< states[i]->getTransProb(j) << endl; \
			} \
		} \
		for (int i = 0; i < numStates; i++) \
		{ \
			for (int j = 0; j < numOutputs; j++) \
			{ \
				os << "Emission probability of state " << i << \
				 ", output " << j << ": " \
				<< states[i]->getEmitProb(j) << endl; \
			} \
		} \
		for (int i = 0; i < numStates; i++) \
		{ \
			os << "Initial probability of state " << i << ": " << \
			startProb[i] << endl; \
		} \
		os << flush;

    void HMMIO::printHMM(ostream& os, State** states, int numStates, int numOutputs, long double* startProb)
	{
		os << "HMM learned for " << numStates << " states:" << endl;
		HMMIO_PRINT_HMM_BODY
	}

    void HMMIO::printHMM(ostream& os, State** states, int numStates, int numOutputs, long double* startProb, int numHMM)
	{
		os << "HMM " << numHMM << " learned for " << numStates << " states:" << endl;
		HMMIO_PRINT_HMM_BODY
	}
}

#undef HMMIO_PRINT_HMM_BODY
