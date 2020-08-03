#include <cstring>
#include <string>

#include "options.hpp"

using std::string;

namespace HMM
{
    class commandLineParser
    {
    public:
        /**
        * Parse command-line options.
        * @param argc the number of arguments counting the first (the path to the executable)
        * @param argv an array of C strings containing the arguments
        * @return a commandLineParser::Options struct containing the options
        */
        static options* parseCommandLine(int argc, char** argv);
    protected:
    private:
        //state as in state machine, not related to HMM state
        enum STATE {AWAITING_OPTION, RECEIVING_PARAM_INT, RECEIVING_PARAM_LONG_DOUBLE, RECEIVING_PARAM_STRING, NO_OPTIONS_REMAIN};
    };
}
