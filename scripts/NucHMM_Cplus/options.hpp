#pragma once

#include <map>
#include <string>

using std::map;
using std::string;

//default paraeters
#define DEFAULT_THRESHOLD 1.0e-6l
#define DEFAULT_BIN_WIDTH 1000
#define DEFAULT_ITERATIONS 300

namespace HMM
{
    class options
    {
        public:
        /**
        * Get the global options instance.
        * @return the options instance
        */
        static options* getOptions();

        /**
        * Delete the global options instance. Call before program exit.
        * @precondition options must be valid.
        */
        static void cleanupOptions();

        /**
        * Sets the chrLengths to the lengths of the hg18 chromosomes.
        */
        static void setHg19ChrLengthArray();

        virtual ~options();

        bool bicEveryIteration;
        bool zeroProbs;
        bool enforceMinProb;
        int chromosomes;
        int numStates;
        int numOutputs;
        int binWidth;
        int iterations;
        int numDataSets;
        long double minProb;
        char** fnames;
        bool usingCustomInitialHMM;
        string initialHMMFile;
        bool outputHMMFile;
        string finalHMMFile;
        bool customGenome;
		bool numbin;
        string customGenomeFile;
		string numbinFile;
        map<string, int> chrNames;
        map<int, string> chrIDs;
        unsigned long* chrLengths;
		unsigned long* chrBins;
        long double minBICChange;
        bool printAllProbs;
        bool valid;

        private:
        options();
        options(const options& rhs){}
        options& operator=(const options& rhs){return *this;}
        static options* instance;
    };
} //namespace HMM
