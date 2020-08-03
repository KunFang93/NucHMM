#pragma once

#include <iostream>
#include <string>

using std::ostream;
using std::string;

namespace HMM
{
    class master
    {
    public:
    // double _thresholdToCallMark = 1.0e-4
        /**
        * Creates a new master object. Uses the options singleton.
        */
        explicit master();

        /**
        * Runs this master on some data. Uses the options singleton.
        * @param outs the ostream to write text output to
        * @return a standard UNIX result code
        */
        int run(ostream& outs);

    private:
        int numStates;
        int numOutputs;
        //bin width
        int interval;
    };
}
