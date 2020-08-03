#pragma once

#include <cmath>
#include <iostream>
#include <string>

#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include "constants.hpp"
#include "macros.hpp"
#include "options.hpp"

using std::ifstream;
using std::ostream;
using std::pow;
using std::string;
using boost::math::factorial;
using boost::math::tgamma;

namespace HMM
{
    class dataLoad
    {
    public:
        /**
        * Loads a bed file of ChIP-seq/MBD-seq reads, binning each read according to its midpoint by default.
        * @param fname the path to the bed file
        * @param interval the width of the bins to put the reads into
        * @param outs the output stream to write output to
        * @param useFivePrime if true, bin each read according to its 5' end (accounting for strand).
        * @param readShift the amount to shift each read in the 5' to 3' direction. Use -1 for bin size / 2.
        * @return the number of reads in each bin, formatted as data[chromosome number][bin]
        */
        static int** readFile(const string& fname, int interval, ostream& outs, int readShift = -1, bool useFivePrime = false);

        /**
        * Counts the total number of reads in a bed file.
        * @param reads the bed file, formatted as data[chromosome number][bin]
        * @param interval the width of the bins
        * @return the total number of reads
        */
        static unsigned long totalReads(int** reads, int interval);

        /**
        * Calculates the average number of marks per bin.
        * @param totalReads the total number of reads
        * @param interval the width of the bins
        * @return the average number of marks
        */
        static long double getMeanMarks(unsigned long totalReads, int interval);

        /**
        * Calculates the threshold to call a mark present or absent based on the discrete Poisson distribution.
        * @param meanMarks the average number of marks
        * @param thresholdToCallMark the maximum p-value under which a mark is to be considered present
        * @return the threshold
        */
        static int getThreshold(long double meanMarks, double thresholdToCallMark);

        /**
        * Calculates the threshold to call a mark present or absent based on an extension of the Poisson distribution to continuous k.
        * @param meanMarks the average number of marks
        * @param thresholdToCallMark the maximum p-value under which a mark is to be considered present
        * @return the threshold
        */
        static long double getContinuousThreshold(long double meanMarks, double thresholdToCallMark);

        /**
        * Gets the distribution of marks calculated from a data file.
        * @param reads the bed file, formatted as data[chromosome number][bin]
        * @param threshold the number-of-reads-per-bin cutoff for calling a mark
        * @param interval the width of the bins
        * @param numBinsMarked optional, will be set to the number of bins marked
        * @return a genome-wide map of marks, formatted as marks[chromosome number][bin]
        */
        static bool** getMarks(int** reads, int threshold, int interval, int* numBinsMarked = NULL);
        /**
        * Gets the distribution of marks calculated from a data file.
        * @param reads the bed file, formatted as data[chromosome number][bin]
        * @param threshold the number-of-reads-per-bin cutoff for calling a mark
        * @param interval the width of the bins
        * @param numBinsMarked optional, will be set to the number of bins marked
        * @return a genome-wide map of marks, formatted as marks[chromosome number][bin]
        */
        static bool** getMarks(double** reads, double threshold, int interval, int* numBinsMarked = NULL);

        /**
        * Convenience method for allocating enough memory to hold a genome-wide data set.
        * @param T the type to allocate the memory with
        * @param interval the width of the bins in base pairs
        * @return the memory, allocated as data[chromosome number][bin]
        */
        template <typename T> _malloc static T** initDataArray(int interval);

        /**
        * Convenience method for allocating enough memory to hold multiple genome-wide data sets
        * @param T the type to allocate the memory with
        * @param dataSets the number of data sets
        * @param interval the width of the bins in base pairs
        * @return the memory, allocated as data[chromosome number + (num chromosomes * data set)][bin]
        */
		template <typename T> _malloc static T** initDataArray_n(unsigned long* numbins);

		/**
		* Convenience method for allocating enough memory to hold multiple genome-wide data sets
		* @param T the type to allocate the memory with
		* @param dataSets the number of data sets
		* @param numbin the num of the bins for each chromosome
		* @return the memory, allocated as data[chromosome number + (num chromosomes * data set)][bin]
		*/
        template <typename T> _malloc static T** initMultipleDataSetDataArray_n(int dataSets, unsigned long* numbins);
		/**
		* Inserts a map of marks for a dataset into a combined data array.
		* @param marks the map of marks
		* @param data the combined array, must already be allocated, existing data will be clobbered
		* @param bit the bit the change
		* @param numbin the num of the bins for each chromosome
		*/

		template <typename T> _malloc static T** initMultipleDataSetDataArray(int dataSets, int interval);

		/**
		* Inserts a map of marks for a dataset into a combined data array.
		* @param marks the map of marks
		* @param data the combined array, must already be allocated, existing data will be clobbered
		* @param bit the bit the change
		* @param interval the width of the bins
		*/

        static void putMarksInData(bool** marks, int** data, int bit, int interval);

        /**
        * Loads a precompiled data file.
        * @param is the input stream
        * @param data the already-allocated array to read the data into, will be overwritten
        */
        static void loadPrecompData(ifstream& is, int** data);

		/**
		* Loads a precompiled data file.
		* @param is the input stream
		* @param data the already-allocated array to read the data into, will be overwritten
		*/
		static void loadPrecompData_out(ifstream& is, int** data, int** start, int** end);

        /**
        * Calculates the Poisson CDF.
        * @param k number of events
        * @param lambda expected number of events
        * @return P(N >= k)
        */
        inline static long double _pure poisson(int k, long double lambda);

        /**
        * Calculates the Poisson CDF, extended to a continuous function k > 0 by removing
        * the floor functions and changing the factorial in the denominator into a gamma.
        * @param k number of events
        * @param lambda expected number of events
        * @return P(N >= k)
        */
        inline static long double _pure continuousPoisson(long double k, long double lambda);

        /**
        * Loads a bed file of ChIP-seq/MBD-seq reads, binning each read accoding to the coverage of the bin by the read.
        * @param fname the path to the bed file
        * @param totalReads will be set to the total number of reads in the file
        * @param outs the output stream to write output to
        * @param interval the width of the bins to put the reads into
        * @return the number of reads in each bin, formatted as data[chromosome number][bin]
        */
        static double** readFileDistributeReads(const string& fname, unsigned long& totalReads, ostream& outs, int interval);

        /**
        * Passes the specified chromosome lengths to the options singleton.
        * @param input the chromosome lengths file
        */
        static void getChrlenFromListFile(ifstream& input);
		
		/**
		* Passes the specified chromosome bin number to the options singleton.
		* @param input the chromosome lengths file
		*/
		static void getChrBinFromListFile(ifstream& input);

    private:
        /**
        * Calculates the mark presence threshold by recursive binary search.
        * @see HMM::dataLoad::getContinuousThreshold(...)
        * @param meanMarks the average number of marks
        * @param lower the current lower bound of the search
        * @param upper the current upper bound of the search
        * @param iterationNum the current recursion depth
        * @param curThreshold is set to the closest threshold found
        * @param thresholdToCallMark the maximum p-value under which a mark is to be considered present
        */
        static long double continuousThresholdRec(long double& meanMarks, long double lower, long double upper, int iterationNum, long double& curThreshold, double& thresholdToCallMark);
    };

	template <typename T> _malloc T** dataLoad::initDataArray(int interval)
	{
		T** arr;
		arr = new T*[options::getOptions()->chromosomes];
		for (int i = 0; i < options::getOptions()->chromosomes; i++)
		{
			unsigned long numBins = (options::getOptions()->chrLengths[i] / interval) + 1;
			arr[i] = new T[numBins];
		}
		return arr;
	}

	template <typename T> _malloc T** dataLoad::initDataArray_n(unsigned long* numbins)
	{
		T** arr;
		arr = new T*[options::getOptions()->chromosomes];
		for (int i = 0; i < options::getOptions()->chromosomes; i++)
		{
			unsigned long numBins = numbins[i % options::getOptions()->chromosomes];;
			arr[i] = new T[numBins];
		}
		return arr;
	}

	template <typename T> _malloc T** dataLoad::initMultipleDataSetDataArray(int dataSets, int interval)
	{
		T** arr;
		arr = new T*[options::getOptions()->chromosomes * dataSets];
		for (int i = 0; i < options::getOptions()->chromosomes * dataSets; i++)
		{
			unsigned long numBins = (options::getOptions()->chrLengths[i % options::getOptions()->chromosomes] / interval) + 1;
			arr[i] = new T[numBins];
		}
		return arr;
	}

	template <typename T> _malloc T** dataLoad::initMultipleDataSetDataArray_n(int dataSets, unsigned long* numbins)
	{
		T** arr;
		arr = new T*[options::getOptions()->chromosomes * dataSets];
		for (int i = 0; i < options::getOptions()->chromosomes * dataSets; i++)
		{
			unsigned long numBins = numbins[i % options::getOptions()->chromosomes];
			arr[i] = new T[numBins];
		}
		return arr;
	}

    inline long double _pure dataLoad::poisson(int k, long double lambda)
    {
//        return (pow(lambda,k) * exp(-lambda)) / (long double) factorial<long double>(k); //Poisson distribution
        return 1.0 - (tgamma<long double>(floorl(k + 1.0), lambda) / (long double) factorial<long double>(floorl(k))); //what wikipedia says is the discrete Poisson cdf
    }

    inline long double _pure dataLoad::continuousPoisson(long double k, long double lambda)
    {
//        return (pow(lambda,k) * exp(-lambda)) / (long double) tgamma<long double>(k + 1.0); //continuous Poisson distribution
        return 1.0 - (tgamma<long double>(k + 1.0, lambda) / tgamma<long double>(k + 1.0)); //what I hope is the continuous Poisson cdf, testing it with WolframAlpha seems to hold up
    }
} //namespace HMM
