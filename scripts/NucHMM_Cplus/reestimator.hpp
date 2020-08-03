#pragma once

#ifdef _OPENMP
#include <omp.h>
#endif

#include "macros.hpp"
#include "stateLD.hpp"

namespace HMM
{
    class reestimator
    {
    public:
        /**
        * Returns a new HMM from the Baum-Welch forward-backward algorithm. Does not spawn threads.
        * @param numStates the number of states in the HMM
        * @param numOutputs the number of outputs of the HMM
        * @param numSequences the number of training observation sequences
        * @param deltas the delta values from the Baum-Welch algorithm, formatted as deltas[seq][state][bin]
        * @param gammaRowSums the sum of the gamma values for each state pair in each observation sequence
        * @param data the observation sequences
        * @param intervals an array of the observation sequence lengths
        * @param startProb will be set to the start probs of the new HMM, must point to appropriately-allocated memory
        * @return the new HMM states
        */
        template <typename T> static State** getNextParams(int numStates, int numOutputs, int numSequences, T*** deltas, T*** gammaRowSums, int** data, unsigned long* intervals, T* startProb);

        /**
        * Returns a new HMM from the Baum-Welch forward-backward algorithm. Uses OpenMP threading.
        * @param numStates the number of states in the HMM
        * @param numOutputs the number of outputs of the HMM
        * @param numSequences the number of training observation sequences
        * @param deltas the delta values from the Baum-Welch algorithm, formatted as deltas[seq][state][bin]
        * @param gammaRowSums the sum of the gamma values for each state pair in each observation sequence
        * @param data the observation sequences
        * @param intervals an array of the observation sequence lengths
        * @param startProb will be set to the start probs of the new HMM, must point to appropriately-allocated memory
        * @return the new HMM states
        */
        template <typename T> static State** getNextParams_threaded(int numStates, int numOutputs, int numSequences, T*** deltas, T*** gammaRowSums, int** data, unsigned long* intervals, T* startProb);
    };

    template <typename T> State** reestimator::getNextParams(int numStates, int numOutputs, int numSequences, T*** deltas, T*** gammaRowSums, int** data, unsigned long* intervals, T* startProb)
    {
        State** HMM = new State*[numStates];
        for (int i = 0; i < numStates; i++)
        {
            HMM[i] = new State(numStates, numOutputs);
        }
        //start probs
        if (startProb != NULL)
        {
            register T denom = 0;
            for (int i = 0; i < numStates; i++)
            {
                for (int seq = 0; seq < numSequences; seq++)
                {
                    denom += deltas[seq][i][0];
                }
            }
            for (int i = 0; i < numStates; i++)
            {
                startProb[i] = 0;
                for (int seq = 0; seq < numSequences; seq++)
                {
                    startProb[i] += deltas[seq][i][0];
                }
                startProb[i] /= denom;
            }
//			std::cerr << "Start prob re-estimation finished.\n";
        }
        //trans probs
        {
            T** gammaDenoms = new T*[numSequences];
            for (int i = 0; i < numSequences; i++)
            {
                gammaDenoms[i] = new T[numStates];
                for (int j = 0; j < numStates; j++)
                {
                    gammaDenoms[i][j] = 0;
//			std::cerr << i << " " << j << "test\n";
                    for (int k = 0; k < numStates; k++)
                    {
                        gammaDenoms[i][j] += gammaRowSums[i][j][k];
                    }
                }
            }
            for (int i = 0; i < numStates; i++)
            {
                for (int j = 0; j < numStates; j++)
                {
                    register T numSeqSum = 0;
                    register T denomSeqSum = 0;
                    for (int seq = 0; seq < numSequences; seq++)
                    {
                        numSeqSum += gammaRowSums[seq][i][j];
                        denomSeqSum += gammaDenoms[seq][i];
                    }
                    HMM[i]->setTransProb(j, numSeqSum / denomSeqSum);
                }
            }
            for (int seq = 0; seq < numSequences; seq++)
            {
                delete[] gammaDenoms[seq];
            }
            delete[] gammaDenoms;
//			std::cerr << "Trans prob re-estimation finished.\n";
        }
        //emit probs
        {
            T* deltaSums = new long double[numStates];
            for (int i = 0; i < numStates; i++)
            {
                deltaSums[i] = 0;
                for (int seq = 0; seq < numSequences; seq++)
                {
                    for (unsigned long j = 0; j < intervals[seq]; j++)
                    {
                        deltaSums[i] += deltas[seq][i][j];
                    }
                }
            }
            for (int i = 0; i < numStates; i++)
            {
                for (int k = 0; k < numOutputs; k++)
                {
                    T sum = 0;
                    for (int seq = 0; seq < numSequences; seq++)
                    {
                        for (unsigned long j = 0; j < intervals[seq]; j++)
                        {
                            if (k == data[seq][j]) sum += deltas[seq][i][j];
                        }
                    }
                    HMM[i]->setEmitProb(k, sum / deltaSums[i]);
                }
            }
            delete[] deltaSums;
//			std::cerr << "Emit prob re-estimation finished.\n";
        }
        return HMM;
    }

    template <typename T> State** reestimator::getNextParams_threaded(int numStates, int numOutputs, int numSequences, T*** deltas, T*** gammaRowSums, int** data, unsigned long* intervals, T* startProb)
    {
        State** HMM = new State*[numStates];
        for (int i = 0; i < numStates; i++)
        {
            HMM[i] = new State(numStates, numOutputs);
        }
        //start probs
        if (startProb != NULL)
        {
            register T denom = 0;
#pragma omp parallel
            {
//std::cout << omp_get_num_threads();
#pragma omp for reduction (+ : denom)
                for (int i = 0; i < numStates; i++)
                {
                    for (int seq = 0; seq < numSequences; seq++)
                    {
                        denom += deltas[seq][i][0];
                    }
                }
#pragma omp for
                for (int i = 0; i < numStates; i++)
                {
                    startProb[i] = 0;
                    for (int seq = 0; seq < numSequences; seq++)
                    {
                        startProb[i] += deltas[seq][i][0];
                    }
                    startProb[i] /= denom;
                }
            }
//			std::cerr << "Start prob re-estimation finished.\n";
        }
        //trans probs
        {
            T** gammaDenoms = new T*[numSequences];
#pragma omp parallel
            {
#pragma omp for
                for (int i = 0; i < numSequences; i++)
                {
                    gammaDenoms[i] = new T[numStates];
                    for (int j = 0; j < numStates; j++)
                    {
                        gammaDenoms[i][j] = 0;
    //			std::cerr << i << " " << j << "test\n";
                        for (int k = 0; k < numStates; k++)
                        {
                            gammaDenoms[i][j] += gammaRowSums[i][j][k];
                        }
                    }
                }
#pragma omp for
                for (int i = 0; i < numStates; i++)
                {
                    for (int j = 0; j < numStates; j++)
                    {
                        register T numSeqSum = 0;
                        register T denomSeqSum = 0;
                        for (int seq = 0; seq < numSequences; seq++)
                        {
                            numSeqSum += gammaRowSums[seq][i][j];
                            denomSeqSum += gammaDenoms[seq][i];
                        }
                        HMM[i]->setTransProb(j, numSeqSum / denomSeqSum);
                    }
                }
            }
            for (int seq = 0; seq < numSequences; seq++)
            {
                delete[] gammaDenoms[seq];
            }
            delete[] gammaDenoms;
//			std::cerr << "Trans prob re-estimation finished.\n";
        }
        //emit probs
        {
            T* deltaSums = new long double[numStates];
#pragma omp parallel
            {
#pragma omp for schedule(dynamic)
                for (int i = 0; i < numStates; i++)
                {
                    deltaSums[i] = 0;
                    for (int seq = 0; seq < numSequences; seq++)
                    {
                        for (unsigned long j = 0; j < intervals[seq]; j++)
                        {
                            deltaSums[i] += deltas[seq][i][j];
                        }
                    }
                }
#pragma omp for schedule(dynamic)
                for (int i = 0; i < numStates; i++)
                {
                    for (int k = 0; k < numOutputs; k++)
                    {
                        T sum = 0;
                        for (int seq = 0; seq < numSequences; seq++)
                        {
                            for (unsigned long j = 0; j < intervals[seq]; j++)
                            {
                                if (k == data[seq][j]) sum += deltas[seq][i][j];
                            }
                        }
                        HMM[i]->setEmitProb(k, sum / deltaSums[i]);
                    }
                }
            }
            delete[] deltaSums;
//			std::cerr << "Emit prob re-estimation finished.\n";
        }
        return HMM;
    }
}

