#include "commandLineParser.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#include "constants.hpp"
#include "options.hpp"

using std::cerr;
using std::endl;
using std::fabs;
using std::strcpy;
using std::string;
using std::strlen;

namespace HMM
{
    options* commandLineParser::parseCommandLine(int argc, char** argv)
    {
        options* options = options::getOptions();
        if (argc == 1)
        {
            cerr << "NucHMM-learn version " << 1.0 << endl;
            cerr << "Usage: NucHMM-learn <options> numStates numOutputs filename1 filename2...\n";
            cerr << "Options:\n";
            cerr << "-b: calculate and report BIC after every iteration\n";
            cerr << "-z: zero probabilities when they drop below the min probability threshold\n";
            cerr << "\tConflicts with -F\n";
            cerr << "-F: do not keep probabilities from dropping below the min probability threshold\n";
            cerr << "\tConflicts with -z\n";
            cerr << "-H [filename]: supply initial HMM parameters. Accepts a .rawhmm file.\n";
            cerr << "\tIf used, the number of states and number of outputs specified on the command line will be overridden by the HMM file.\n";
            cerr << "-h [filename]: output the final HMM in .rawhmm format. Generally recommended.\n";
            cerr << "-m [double]: minimum transition, emission or start probability (default: 1.0e-6)\n";
            cerr << "-i [int]: maximum number of iterations (default: 300)\n";
            cerr << "\tUse -1 to train continuously (usually used with -d)\n";
            cerr << "-w [int]: bin width (default: 1000)\n";
            cerr << "-d [double]: stop training if |decrease of BIC| is below this value\n";
            cerr << "\tImplies -b.\n";
            cerr << "-g [filename]: use a custom genome file (default genome is hg18). This file must have been specified to the precompiler as well.\n";
			cerr << "-v [filename]: the number of bins for each chromosome.\n";
            cerr << "-p: print transition, emission, and start probabilities of the HMM after each iteration.\n";
            return options;
        }
        typedef commandLineParser::STATE state_t;
        state_t state;
        state = AWAITING_OPTION;
        int locToStartRequiredOptions = -1;
        int* intToFill = NULL;
        long double* longDoubleToFill = NULL;
        string* stringToFill = NULL;
        for (int i = 1; i < argc; i++)
        {
            switch(state)
            {
                case AWAITING_OPTION:
                {
                    if (argv[i][0] == '-')
                    {
                        //this will work even if argv[i] == "-"; it will grab the terminating \0
                        char option = argv[i][1];
                        switch(option)
                        {
                            case 'b':
                            {
                                options->bicEveryIteration = true;
                            } break;
                            case 'z':
                            {
                                if (!options->enforceMinProb)
                                {
                                    cerr << "Error: -F and -z cannot be combined.\n";
                                    return options;
                                }
                                options->zeroProbs = true;
                            } break;
                            case 'm':
                            {
                                longDoubleToFill = &options->minProb;
                                state = RECEIVING_PARAM_LONG_DOUBLE;
                            } break;
                            case 'i':
                            {
                                intToFill = &options->iterations;
                                state = RECEIVING_PARAM_INT;
                            } break;
                            case 'w':
                            {
                                intToFill = &options->binWidth;
                                state = RECEIVING_PARAM_INT;
                            } break;
                            case 'F':
                            {
                                if (options->zeroProbs)
                                {
                                    cerr << "Error: -F and -z cannot be combined.\n";
                                    return options;
                                }
                                options->enforceMinProb = false;
                            } break;
                            case 'H':
                            {
                                options->usingCustomInitialHMM = true;
                                stringToFill = &options->initialHMMFile;
                                state = RECEIVING_PARAM_STRING;
                            } break;
                            case 'h':
                            {
                                options->outputHMMFile = true;
                                stringToFill = &options->finalHMMFile;
                                state = RECEIVING_PARAM_STRING;
                            } break;
                            case 'g':
                            {
                                options->customGenome = true;
                                stringToFill = &options->customGenomeFile;
                                state = RECEIVING_PARAM_STRING;
                            } break;
							case 'v':
							{
								options->numbin = true;
								stringToFill = &options->numbinFile;
								state = RECEIVING_PARAM_STRING;
							} break;
                            case 'd':
                            {
                                options->bicEveryIteration = true;
                                longDoubleToFill = &options->minBICChange;
                                state = RECEIVING_PARAM_LONG_DOUBLE;
                            } break;
                            case 'p':
                            {
                                options->printAllProbs = true;
                            } break;
                            default:
                            {
                                cerr << "Error: invalid option " << argv[i] << "\n";
                                return options;
                            } break;
                        }
                    }
                    else
                    {
                        state = NO_OPTIONS_REMAIN;
                        locToStartRequiredOptions = i;
                        i = argc;
                    }
                } break;
                case RECEIVING_PARAM_LONG_DOUBLE:
                {
                    *longDoubleToFill = atof(argv[i]);
                    state = AWAITING_OPTION;
                } break;
                case RECEIVING_PARAM_INT:
                {
                    *intToFill = atoi(argv[i]);
                    state = AWAITING_OPTION;
                } break;
                case RECEIVING_PARAM_STRING:
                {
                    *stringToFill = argv[i];
                    state = AWAITING_OPTION;
                } break;
                default:
                {
                    cerr << "Internal error: invalid or unexpected state reached in command line parser.\n";
                    return options;
                } break;
            }
        }
        if(locToStartRequiredOptions == -1 || locToStartRequiredOptions + 2 >= argc)
        {
            //locToStartRequiredOptions remains -1 if the parser has run out of options before expecting to start reading the following parameters
            int missingOptions = (locToStartRequiredOptions < 0)? -1 : (locToStartRequiredOptions + 2 - argc);
            switch(missingOptions)
            {
                //this deliberately uses fall-through
                case (-1):
                case (2):
                {
                    cerr << "Error: please supply number of states in desired HMM.\n";
                }
                case (1):
                {
                    cerr << "Error: please supply number of outputs in data file(s).\n";
                }
                case (0):
                {
                    cerr << "Error: missing data file(s).\n";
                } break;
                default:
                {
                    cerr << "Error: missing at least one of: numStates numOutputs dataFile(s)\n";
                } break;
            }
            return options;
        }
        options->numStates = atoi(argv[locToStartRequiredOptions]);
        options->numOutputs = atoi(argv[locToStartRequiredOptions + 1]);
        int beginFnames = locToStartRequiredOptions + 2;
        int numDataSets =  argc - beginFnames;
        options->fnames = new char*[numDataSets];
        options->numDataSets = numDataSets;
        for (int i = beginFnames; i < argc; i++)
        {
            options->fnames[i - beginFnames] = new char[strlen(argv[i]) + 1];
            strcpy(options->fnames[i - beginFnames], argv[i]);
        }
        options->valid = true;
        //sanity checks
        if (options->minBICChange == 0)
        {
            cerr << "Error: minimum BIC change cannot be 0.\n";
            options->valid = false;
        }
        if (options->minBICChange != -1) options->minBICChange = fabs(options->minBICChange);
        if (options->iterations != -1 && options->iterations < 1)
        {
            cerr << "Error: number of iterations less than 1 and not -1.\n";
            options->valid = false;
        }
        if (options->binWidth < 1)
        {
            cerr << "Error: bin size less than 1.\n";
            options->valid = false;
        }
        if (options->minProb <= 0.0l || options->minProb >= 1.0l / (long double) options->numStates)
        {
            cerr << "Error: minimum probability must be between 0 and 1 / (number of states) (exclusive).\n";
            options->valid = false;
        }
        return options;
    }
}
