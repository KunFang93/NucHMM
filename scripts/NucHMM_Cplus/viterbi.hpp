#pragma once

#include <iostream>
#include <fstream>
#include <omp.h>

#include "constants.hpp"
#include "lnsnum.hpp"
#include "macros.hpp"
#include "stateLD.hpp"

using std::cout;
using std::endl;
using std::ofstream;

namespace HMM
{
    class viterbiAlg
    {
    public:
        /**
        * Implements the Viterbi decoder algorithm for calling states in an observation sequence from a Markov model.
        * @param T the data type of the model probabilities
        * @param numStates the number of states in the model
        * @param numOutputs the number of outputs of the model
        * @param states the states of the model
        * @param startProb the start probs of the model
        * @param data the observation sequence
        * @param intervals the length of the observation sequence
        * @param numThreads 1 + the number of threads this method is allowed to spawn
        * @return the most likely state sequence
        */
        template<typename T> static int* viterbi(int numStates, int numOutputs, State** states, T* startProb, int* data, unsigned long intervals);
    };

	template<typename T> int* viterbiAlg::viterbi(int numStates, int numOutputs, State** states, T* startProb, int* data, unsigned long intervals)
	{
		int* stateMap = new int[intervals];
		Lnsnum** delta = new Lnsnum*[numStates];
		int** psi = new int*[numStates];
		for (int i = 0; i < numStates; i++)
		{
			delta[i] = new Lnsnum[intervals];
			delta[i][0] = startProb[i] * states[i]->getEmitProb(data[0]);
			psi[i] = new int[intervals];
			psi[i][0] = 0;
		}
		//forward sweep
		for (unsigned long t = 1; t < intervals; t++)
		{
			for (int j = 0; j < numStates; j++)
			{
				Lnsnum tempMax = delta[0][t-1] * states[0]->getTransProb(j);
				int maxIndex = 0;
				for (int i = 1; i < numStates; i++)
				{
					Lnsnum temp = delta[i][t-1] * states[i]->getTransProb(j);
					if (temp > tempMax)
					{
						tempMax = temp;
						maxIndex = i;
					}
				}
				delta[j][t] = tempMax * states[j]->getEmitProb(data[t]);
				psi[j][t] = maxIndex;
			}
		}
		//get q(tMax)
		Lnsnum tempMax = delta[0][intervals-1];
		int maxIndex = 0;
		for (int i = 1; i < numStates; i++)
		{
			if (delta[i][intervals-1] > tempMax)
			{
				tempMax = delta[i][intervals-1];
				maxIndex = i;
			}
		}
		stateMap[intervals - 1] = maxIndex;
		for (unsigned long t = intervals - 1; t > 0; t--) //must use --t since it is unsigned
		{
			stateMap[t-1] = psi[stateMap[t]][t];
		}
		//free memory
		for (int i = 0; i < numStates; i++)
		{
			delete[] psi[i];
			delete[] delta[i];
		}
		delete[] psi;
		delete[] delta;
		return stateMap;
	}
} //namespace HMM
