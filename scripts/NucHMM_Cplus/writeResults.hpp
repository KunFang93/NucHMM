#pragma once

#include <iostream>
#include <string>

#include "macros.hpp"

using std::ofstream;
using std::string;

namespace HMM
{
    class resultsWriter
    {
    public:
        /**
        * Writes genome-wide results files. A output map will be written ifdef HMM_WRITE_OUTPUT_MAPS, state map ifdef HMM_WRITE_STATE_MAPS
        * @param stateMap a genome-wide map of states
        * @param outputMap a genome-wide map of outputs
        * @param stateFname the filename of the state file to write
        * @param outputFname the filename of the output file to write
        * @param interval the width of the bins
        */
        static void writeResults(int** controlStateMap, int** controlOutputMap, string stateFname, string outputFname, int interval);
		/**
		* Writes genome-wide results files. A output map will be written ifdef HMM_WRITE_OUTPUT_MAPS, state map ifdef HMM_WRITE_STATE_MAPS
		* @param stateMap a genome-wide map of states
		* @param outputMap a genome-wide map of outputs
		* @param stateFname the filename of the state file to write
		* @param outputFname the filename of the output file to write
		* @param numbins the number of the bins
		*/
		static void writeResults_n(int** controlStateMap, int** controlOutputMap, string stateFname, string outputFname, unsigned long* numbins, int** startposition, int** endposition);
        /**
        * Writes a genome-wide map of called marks. Format is [chr#]\t[bin#]\t[data]. Chr# is the chromosome number starting from ZERO (chrX is 22 in humans). Bin# start from zero.
        * Data is the number resulting from the combination of marks expressed in binary form. For example, if there are three possible marks A, B, C, 6 = 0b110 and means that marks
        * B and C are present.
        * @param os the output stream to write to
        * @param data the data map
        * @param interval the width of the bins
        */
        static void writePrecompData(ofstream& os, int** data, int interval);
    };
} //namespace HMM
