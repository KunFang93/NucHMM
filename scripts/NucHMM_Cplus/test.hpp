#pragma once

//#include "macros.hpp"

namespace HMM
{
    /**
    * A test driver for the various HMM-related algorithms. Randomly initializes a model, performs Baum-Welch training,
    * test the Poisson distribution code and other things.
    * It is important that the BIC score decreases with each training iteration!
    */
    class testLearn
    {
    public:
        /**
        * Creates a new TestLearn object.
        */
        explicit testLearn();

        /**
        * Runs the test code.
        * @param iterations the number of iterations for the Baum-Welch training test
        */
        int test(int iterations = 100);
    };
}
