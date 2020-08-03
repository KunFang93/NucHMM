#pragma once

#include "stateTemplate.hpp"

namespace HMM
{
    /**
    * Specialization of StateTemplate<> for long double.
    */
    class StateLD : public StateTemplate<long double>
    {
    public:
        StateLD(int numStates, int numOutputs) : StateTemplate<long double>(numStates, numOutputs)
        {};

        StateLD(const StateLD& rhs) : StateTemplate<long double>(rhs)
        {};
    };

//I know this shouldn't really exist, but I really don't feel like hunting down every reference to State
    /**
    * @deprecated
    * Specialization of StateTemplate<> for long double.
    */
    class State : public StateTemplate<long double>
    {
    public:
        State(int numStates, int numOutputs) : StateTemplate<long double>(numStates, numOutputs)
        {};

        State(const State& rhs) : StateTemplate<long double>(rhs)
        {};
    };
}
