#pragma once

#include "constants.hpp"

#include <iostream>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include "macros.hpp"

using std::ostream;
using boost::mt19937;
using boost::uniform_real;
using boost::variate_generator;

namespace HMM
{
    class Lnsnum;
    class State;

    class masterUtils
    {
    private:
    	struct DValNode
        {
            double val;
            DValNode* next;
        };

    public:
        /**
        * Calculates a BIC score for an HMM.
        * @param numStates the number of states in the HMM
        * @param numOutputs the number of outputs of the HMM
        * @param prH the PrH values of the observation sequences that the HMM was trained with
        * @param totalIntervals the total length of all the observation sequences that the HMM was trained with
        * @return the BIC score
        */
        static long double _pure calcBIC(int numStates, int numOutputs, Lnsnum* prH, unsigned long totalIntervals);
        /**
        * Calculates a BIC score for an HMM.
        * @param numStates the number of states in the HMM
        * @param numOutputs the number of outputs of the HMM
        * @param logLikelihood the natural log of the total probability of the model over the data
        * @param totalIntervals the total length of all the observation sequences that the HMM was trained with
        * @return the BIC score
        */
        static long double _pure calcBIC(int numStates, int numOutputs, long double logLikelihood, unsigned long totalIntervals);

        /**
        * Sets probabilities in an HMM to zero if it is below a specified threshold (distributing the probability evenly among all other non-zero probabilities of the same state),
        * attempting to guarantee that the model remains fully connected and that all outputs can be emitted by the HMM. Note that this is not absolutely guaranteed as there
        * is a change that for example output X can be emitted only by state 3, that only state 4 can transition to state 3, and that a model would be in state 2 from the previous
        * observation, in which case the Viterbi and forward-backward algorithms fail. This is most likely with HMMs with small numbers of states, and occasionally may require multiple
        * attempts at HMM learning.
        * @param states the states of the HMM
        * @param numStates the number of states in the HMM
        * @param numDataEntries the number of possible probs to consider in one state for the parameter being considered
        * @param State::*getDatum a pointer to a member of HMM::State that gets a desired probability
        * @param State::*setDatum a pointer to a member of HMM::State that sets a desired probability
        * @param minProbToHandle the threshold under which probabilties are zeroa
        * @param numConnectionsToPreserve for each state/output, the minimum number of states that must transition to it or emit it
        */
        static void zeroProbs(State** states, int numStates, int numDataEntries, long double (State::*getDatum) (int) const, void (State::*setDatum) (int, long double), long double minProbToHandle, int numConnectionsToPreserve = 0);

        /**
        * Keeps states or outputs of an HMM from dropping below minProbToHandle, normalizing probs to keep their sum equal to 1. This guarantees that the model remains fully connected,
        * can emit any output in any state and that probabilities do not underflow.
        * @param states the states of the HMM
        * @param numStates the number of states in the HMM
        * @param numOutputs the number of outputs of the HMM
        * @param startProbs the start probs of the HMM
        * @param minProbToHandle the minimum enforced by this function for the HMM probabilities
        */
        static void forceMinimumProbs(State** states, int numStates, int numOutputs, long double* startProbs, long double minProbToHandle);

        //this should be edited to accept any RNG beside mt19937
        /**
        * Generates an HMM with randomly initialized, normalized parameters.
        * @param numStates the desired number of states for the HMM
        * @param numOutputs the desired number of outptus for the HMM states
        * @param rand a Boost mt19937 RNG to use for random initialization (@see boost::variate_generator)
        * @param states will be set to the states of the HMM
        * @param startProb will be set to the start probs of the HMM
        */
        static void generateRandomHMM(int numStates, int numOutputs, variate_generator<mt19937, uniform_real<> >& rand, State**& states, long double*& startProb);

        /**
        * Calculates a table of state difference scores.
        * @param numStatesInWM the number of states in the model
        * @param numOutputs the number of outputs in the model
        * @param states the states of the HMM
        * @return the table
        */
        static long double** getStateDifferenceTable(int numStatesInWM, int numOutputs, State** states);

/*
		//encapsulate the raw C gsl_stats_correlation function
		/ **
		* Calculates the Pearson correlation coefficient. A wrapper around GSL's gsl_states_correlation function.
		* @param x a list of values
		* @param y another list of values
		* @param length the length of the value lists
		* @return the Pearson correlation coefficient
		* /
		static double _pure pearson(double* x, double* y, int length);
*/

        /**
        * Convenience method for freeing the memory occupied by the states of an HMM.
        * @param numStates the number of states in the HMM
        * @param HMM the states in the HMM - will be left invalid
        */
        static void deleteHMMStates(int numStates, State**& HMM);
    };

} //namespace HMM
