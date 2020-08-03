#pragma once
/*
#ifdef _OPENMP
#include <omp.h>
#endif
*/
#include "constants.hpp"
#include "macros.hpp"
#include "stateLD.hpp"

namespace HMM
{
    class processChr
    {
    public:
        template <typename T> static void baumWelch_canonical_scaled(int numOutputs, int numStates, State* __restrict__ * __restrict__ HMM, T* __restrict__ startProb, int* __restrict__ data, unsigned long intervals, T**& gammaRowSums, T* __restrict__ * __restrict__ alphasOrDeltas);
    };

	template <typename T> void processChr::baumWelch_canonical_scaled(int numOutputs, int numStates, State* __restrict__ * __restrict__ HMM, T* __restrict__ startProb, int* __restrict__ data, unsigned long intervals, T**& gammaRowSums, T* __restrict__ * __restrict__ alphasOrDeltas)
	{
		State* __restrict__ * __restrict__ states = HMM;
		//perhaps cache this?
		T* __restrict__ * __restrict__ *  __restrict__ stateTransEmit = new T* __restrict__ * [numOutputs];
//#ifdef _OPENMP
//		omp_set_num_threads(numThreads);
//#endif
//#pragma omp parallel for
		for (int k = 0; k < numOutputs; k++)
		{
			stateTransEmit[k] = new T* __restrict__ [numStates];
			for (int i = 0; i < numStates; i++)
			{
				stateTransEmit[k][i] = new T[numStates];
				for (int j = 0; j < numStates; j++)
				{
					stateTransEmit[k][i][j] = states[i]->getTransProb(j) * states[j]->getEmitProb(k);
				}
			}
		}

		//I don't think this is truly the PrH (as it would be without scaling), but it gets set to the proper value for the HMM calculations
		T prH = 0;
		T* scaleFactors = new T[intervals];
		T** betas = new T*[numStates];
		T*** gammas = new T**[numStates];
		for (int i = 0; i < numStates; i++)
		{
			betas[i] = new T[intervals];
			gammas[i] = new T*[numStates];
			for (int j = 0; j < numStates; j++)
			{
				gammas[i][j] = new T[intervals - 1];
				gammaRowSums[i][j] = 0;
			}
		}
		//forward (alphas)
		scaleFactors[0] = 0;
		for (int i = 0; i < numStates; i++)
		{
			alphasOrDeltas[i][0] = startProb[i] * states[i]->getEmitProb(data[0]);
			scaleFactors[0] += alphasOrDeltas[i][0];
//			cout << alphasOrDeltas[i][0] << "\n";
		}
		for (int i = 0; i < numStates; i++)
		{
			alphasOrDeltas[i][0] /= scaleFactors[0];
		}
		for (unsigned long j = 1; j < intervals; j++)
		{
			scaleFactors[j] = 0;
//#pragma omp parallel
			{
//#pragma omp for
				for (int s = 0; s < numStates; s++)
				{
					T tempSum = 0;
					for (int t = 0; t < numStates; t++)
					{
						tempSum += alphasOrDeltas[t][j-1] * stateTransEmit[data[j]][t][s]; //alphasOrDeltas[t][j-1] * states[t]->getTransProb(s) * states[s]->getEmitProb(data[j]);
					}
					alphasOrDeltas[s][j] = tempSum;
					scaleFactors[j] += tempSum;
//					cout << tempSum << "\n";
				}
//#pragma omp for
				for (int s = 0; s < numStates; s++)
				{
					alphasOrDeltas[s][j] /= scaleFactors[j];
				}
			}
		}
//		std::cerr << "Forward finished.\n";
		//PrH
		for (int i = 0; i < numStates; i++)
		{
			prH += alphasOrDeltas[i][intervals - 1];
		}
//		cout << prH << std::endl;
//		std::cerr << "PrH finished.\n";
		//backward (betas)
		for (int i = 0; i < numStates; i++)
		{
			betas[i][intervals - 1] = 1 / scaleFactors[intervals - 1];
		}
		//j + 1 used because it is unsigned, and checking >= 0 would always be true!
		for (unsigned long jPlusOne = intervals - 1; jPlusOne > 0; jPlusOne--)
		{
//#pragma omp parallel for
			for (int curState = 0; curState < numStates; curState++)
			{
				T tempSum = 0;
				for (int u = 0; u < numStates; u++)
				{
					tempSum += stateTransEmit[data[jPlusOne]][curState][u] * betas[u][jPlusOne]; // states[curState]->getTransProb(u) * states[u]->getEmitProb[data[jPlusOne]]
				}
				betas[curState][jPlusOne - 1] = tempSum / scaleFactors[jPlusOne - 1];
//				cout << tempSum << "\n";
			}
		}
//		std::cerr << "Backward finished.\n";
		//EM (gammas)
		for (unsigned long j = 0; j < intervals - 1; j++)
		{
//#pragma omp parallel for
			for (int s = 0; s < numStates; s++)
			{
				for (int t = 0; t < numStates; t++)
				{
					gammas[s][t][j] = (alphasOrDeltas[s][j] * stateTransEmit[data[j+1]][s][t] * betas[t][j+1]) / prH;
//					cout << gammas[s][t][j] << "\n";
				}
			}
		}
//		std::cerr << "Gammas finished.\n";
		//deltas
		for (int s = 0; s < numStates; s++)
		{
			alphasOrDeltas[s][intervals - 1] /= prH;
		}
		//gamma row sums
		for (unsigned long j = 0; j < intervals - 1; j++)
		{
//#pragma omp parallel for
			for (int s = 0; s < numStates; s++)
			{
				alphasOrDeltas[s][j] = 0;
				for (int t = 0; t < numStates; t++)
				{
					gammaRowSums[s][t] += gammas[s][t][j];
					alphasOrDeltas[s][j] += gammas[s][t][j];
				}
			}
//		std::cerr << "Deltas and gamma row sums finished.\n";
		}
		//free memory
		for (int i = 0; i < numStates; i++)
		{
			for (int j = 0; j < numStates; j++)
			{
				delete[] gammas[i][j];
			}
			delete[] gammas[i];
			delete[] betas[i];
		}
		delete[] betas;
		delete[] gammas;
		delete[] scaleFactors;
		for (int i = 0; i < numOutputs; i++)
		{
			for (int j = 0; j < numStates; j++)
			{
				delete[] stateTransEmit[i][j];
			}
			delete[] stateTransEmit[i];
		}
		delete[] stateTransEmit;
//		std::cerr << "Memory free finished.\n";
	} // void processChromosome

} //namespace HMM
