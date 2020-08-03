#pragma once

#include <iostream>

#include "macros.hpp"

using std::istream;
using std::ostream;

namespace HMM
{
    class State;

    /**
    * Contains utility methods for loading and saving HMMs.
    */
    class HMMIO
    {
    public:
        /**
        * Outputs an HMM into a standard form.
        * @param os the stream to write to
        * @param states an array of State* that makes up the HMM
        * @param numStates the number of states in the HMM
        * @param numOutputs the number of outputs of the HMM
        * @param startProbs the start probs of the HMM
        */
        static void outputHMM(ostream& os, State** states, int numStates, int numOutputs, long double* startProbs);

        /**
        * Inputs an HMM in the form used by outputHMM.
        * @see outputHMM()
        * @param is the stream to read from
        * @param states will point to an array of State* that will be populated with states
        * @param numStates set to the number of states in the HMM
        * @param numOutputs set to the number of outputs of the HMM
        * @param startProbs will point to the start prob array
        */
        static void inputHMM(istream& is, State**& states, int& numStates, int& numOutputs, long double*& startProbs);

        /**
        * Oputputs an HMM in a human-readable format.
        * @param os the stream to write to
        * @param states an array of states that makes up the HMM
        * @param numStates the number of states in the HMM
        * @param numOutputs the number of outputs of the HMM
        * @param startProbs the start probs of the HMM
        */
        static void printHMM(ostream& os, State** states, int numStates, int numOutputs, long double* startProb);

        /**
        * Oputputs an HMM in a human-readable format.
        * @param os the stream to write to
        * @param states an array of states that makes up the HMM
        * @param numStates the number of states in the HMM
        * @param numOutputs the number of outputs of the HMM
        * @param startProbs the start probs of the HMM
        * @param numHMM a numerical ID for the HMM
        */
        static void printHMM(ostream& os, State** states,int numStates, int numOutputs, long double* startProb, int numHMM);
    };
} //namespace HMM
