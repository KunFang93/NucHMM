#pragma once

// #define HMM_GETPRH_CPP_OUTPUT_DEBUG

#include "constants.hpp"
#include "lnsnum.hpp"
#include "macros.hpp"
#include "stateLD.hpp"

using std::cout;
using std::endl;

namespace HMM
{
	class prHCalculator
	{
	public:
        /** Calculates the PrH value of an HMM over an observation sequence.
        * @param numStates the number of states in the HMM
        * @param states the states of the HMM
        * @param startProb the start probs of the HMM
        * @param data one of the observation sequences used to train the HMM
        * @param intervals the number of bins in the sequence
        * @return the PrH value
        */
		template<typename T> static Lnsnum getPrH(int numStates, State** states, T* startProb, int* data, unsigned long intervals);

	private:
        /** Gets the last alpha values for the last observation in an observation sequence.
        * @param states the states of the HMM
        * @param chr the observation sequence
        * @param numStates the number of states in the HMM
        * @param startProb the start probs of the HMM
        * @param intervals the number of bins in the sequence
        * @return the alpha values
        */
		template<typename T> static Lnsnum* getLastAlphas(State** states, int* chr, int numStates, T* startProb, unsigned long intervals);
	};

	template<typename T> Lnsnum prHCalculator::getPrH(int numStates, State** states, T* startProb, int* data, unsigned long intervals)
	{
		Lnsnum* alpha;
		Lnsnum prH;

		//get alpha and prH values
#ifdef HMM_GETPRH_CPP_OUTPUT_DEBUG
		cout << "Processing alpha values..." << endl;
#endif
		alpha = getLastAlphas(states, data, numStates, startProb, intervals);
		prH = 0;
		for (int j = 0; j < numStates; j++)
		{
			prH += alpha[j];
		}
#ifdef HMM_GETPRH_CPP_OUTPUT_DEBUG
		cout << "log(prH): " << prH.getLog() << endl;
#endif
		delete[] alpha;
		return prH;
	}

	template<typename T> Lnsnum* prHCalculator::getLastAlphas(State** states, int* chr, int numStates, T* startProb, unsigned long intervals)
	{
		Lnsnum* tempAlpha = new Lnsnum[numStates];
		Lnsnum* alpha = new Lnsnum[numStates];
		Lnsnum* tempPtr;
		for (int i = 0; i < numStates; i++)
		{
			alpha[i] = startProb[i] * states[i]->getEmitProb(chr[0]);
		}
		for (unsigned long j = 1; j < intervals; j++)
		{
//#pragma omp parallel for
			for (int k = 0; k < numStates; k++)
			{
				tempAlpha[k] = 0;
				for (int l = 0; l < numStates; l++)
				{
					tempAlpha[k] += alpha[l] * states[l]->getTransProb(k) * states[k]->getEmitProb(chr[j]);
				}
			}
			tempPtr = tempAlpha;
			tempAlpha = alpha;
			alpha = tempPtr;
		}
#ifdef HMM_GETPRH_CPP_OUTPUT_DEBUG
for (int k = 0; k < numStates; k++)
{
	for (unsigned long j = 0; j < intervals; j++)
	{
		cout << "alpha(" << k << ", " << j << "): " << log(alpha[k][j]) << " " << alpha[k][j] << endl;
	}
	cout    << "state " << k << ": " << states[k]->getTransProb(0) << " " << states[k]->getTransProb(1) << " "
		<< states[k]->getEmitProb(0) << " " << states[k]->getEmitProb(1) << endl;
}
#endif
		delete[] tempAlpha;
		return alpha;
	}


} //namespace HMM
