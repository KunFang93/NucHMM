#pragma once

#define VERSION "1.0"

//preprocessor tags
//write state maps in final step
#define HMM_WRITE_STATE_MAP
//write output maps in final step
#define HMM_WRITE_OUTPUT_MAP
//learn the HMM over the Y chromosome (only if hg18 default is used)
#define USE_Y_CHROMOSOME
//exclude a read if its strand is unknown (default: assume + strand)
//#define EXCLUDE_UNKNOWN_STRAND
//use a continuous (gamma-based) instead of discrete (factorial-based) Poisson distribution to calculate read thresholds for calling marks
//#define CONTINUOUS_MARK_THRESHOLD
//warn on invalid chromsome name in preprocessing (such as chrM)
//#define WARN_INVALID_CHROMOSOME

//only use one (or none) of these!
//shift reads in their 5' to 3' direction by binSize / 2
//not supported if CONTINUOUS_MARK_THRESHOLD is also specified
//#define SHIFT_READS_ERNST_KELLIS
//shift reads by half their average fragment size as estimated by BELT
//not supported if CONTINUOUS_MARK_THRESHOLD is also specified
#define SHIFT_READS_BELT
//default estimated fragment length if BELT algorithm is unsucessful
#define BELT_DEFAULT_FRAG_LEN 150

namespace HMM
{
    //when CONTINUOUS_MARK_THRESHOLD is defined, how many iterations to search for k such that P == thresholdToCallMark
	const int numThresholdSearchIterations = 15;

//hg18
#define CHROMOSOMES_HG19 24

	const unsigned long hg19ChrLengths[] = {
		249250621,
		243199373,
		198022430,
		191154276,
		180915260,
		171115067,
		159138663,
		146364022,
		141213431,
		135534747,
		135006516,
		133851895,
        115169878,
        107349540,
        102531392,
        90354753,
        81195210,
        78077248,
        59128983,
        63025520,
        48129895,
        51304566,
        155270560,
        59373566,
        16571
	};

} //namespace HMM
