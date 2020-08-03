#if !(defined(unix) || defined(__unix__) || defined(__unix))
#error This program only supports POSIX environments.
#endif

#include <cstdlib>
#include <iostream>
#include <string>

#include "constants.hpp"
#include "macros.hpp"

//#define TEST_DRIVER

#ifdef TEST_DRIVER
#include "test.hpp"
#else
#include "commandLineParser.hpp"
#include "master.hpp"
#include "options.hpp"
#endif

#ifdef TEST_DRIVER
using HMM::testLearn;
#else
using HMM::commandLineParser;
using HMM::master;
using HMM::options;
#endif

using std::cout;
using std::endl;
using std::string;

int main(int argc, char* argv[])
{
#ifdef TEST_DRIVER
    if (argc == 1)
    {
        cout << "Test driver for HMM_MPI_LEARN. Usage: HMM_MPI_LEARN #iterations" << endl;
        return(EXIT_FAILURE);
    }
    int iterations = atoi(argv[1]);
    testLearn t;
    return t.test(iterations);
#else
    commandLineParser::parseCommandLine(argc, argv);
    if (options::getOptions()->valid == false) return(EXIT_FAILURE);
    master m;
    int resultCode = m.run(cout);
    options::cleanupOptions();
    return resultCode;
#endif
}
