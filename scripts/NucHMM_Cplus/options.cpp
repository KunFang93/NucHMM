#include "options.hpp"

#include <sstream>
#include <string>

#include "constants.hpp"

using std::ostringstream;
using std::string;

namespace HMM
{
    options* options::instance = NULL;

    options::options()
    : bicEveryIteration(false)
    , zeroProbs(false)
    , enforceMinProb(true)
    , chromosomes(0)
    , numStates(0)
    , numOutputs(0)
    , binWidth(DEFAULT_BIN_WIDTH)
    , iterations(DEFAULT_ITERATIONS)
    , numDataSets(0)
    , minProb(DEFAULT_THRESHOLD)
    , fnames(NULL)
    , usingCustomInitialHMM(false)
    , initialHMMFile()
    , finalHMMFile()
    , customGenome(false)
    , customGenomeFile()
    , chrNames()
    , chrIDs()
    , chrLengths(NULL)
	, chrBins(NULL)
    , minBICChange(-1.0l)
    , printAllProbs(false)
    , valid(false) //must be made valid
    {}

    //assuming this is consistent!
    options::~options()
    {
        if (fnames != NULL)
        {
            for (int i = 0; i < numDataSets; i++)
            {
                delete[] fnames[i];
            }
            delete[] fnames;
        }
        if (chrLengths != NULL) delete[] chrLengths;
    }

    options* options::getOptions()
    {
        if (instance == NULL)
        {
            instance = new options();
        }
        return instance;
    }

    void options::cleanupOptions()
    {
//        assert(Options::instance->valid);
        delete instance;
        instance = NULL;
    }

	void options::setHg19ChrLengthArray()
	{
        if (instance == NULL)
        {
            instance = new options();
        }
        instance->chromosomes = CHROMOSOMES_HG19;
#ifndef USE_Y_CHROMOSOME
	    instance->chromosomes--; //Y chromosome excluded by default
#endif
        unsigned long* arr = new unsigned long[instance->chromosomes];
	    for (int i = 0; i < instance->chromosomes; i++)
	    {
	        arr[i] = hg19ChrLengths[i];
	    }
	    instance->chrLengths = arr;
	    //based off of similar code in BELT
        for (int i = 0; i < 22; i++)
        {
            ostringstream oss;
            oss << "chr" << (i + 1);
            string chr = oss.str();
            instance->chrNames[chr] = i;
	    instance->chrIDs[i] = chr;
        }
        string chrX("chrX");
        instance->chrNames[chrX] = 22;
        instance->chrIDs[22] = chrX;
#ifdef USE_Y_CHROMOSOME
        string chrY("chrY");
        instance->chrNames[chrY] = 23;
        instance->chrIDs[23] = chrY;
#endif
	}
} //namespace HMM
