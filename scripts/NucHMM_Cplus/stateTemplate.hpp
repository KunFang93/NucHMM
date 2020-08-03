#pragma once

#include <cmath>

#include <sys/time.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

//deliberately not including constants.hpp for genericity

using std::pow;
using boost::mt19937;
using boost::uniform_real;
using boost::variate_generator;

namespace HMM
{

    /**
    * A state in an HMM
    * @param T the data type of this state's trans and emit probs, must support the usual operator overloads
    */
    template <typename T> class StateTemplate
    {
    public:
        /**
        * Creates a new StateTemplate
        * @param numStates the number of states in the HMM this StateTemplate is part of
        * @param the number of outputs of the state
        */
        StateTemplate(int numStates, int numOutputs);
        /**
        * Copy constructor for StateTemplate
        * @param rhs the other StateTemplate
        */
        StateTemplate(const StateTemplate<T>& rhs); //we have to be sure we're making deep copies

        virtual ~StateTemplate();

        /**
        * Initializes this state with random (boost mt19937) transition and emission parameters
        */
        void initRandom();

        /**
        * Returns the transition probability of this state to another state in the model
        * @param id the number of the other state
        * @return the trans prob
        */
        T getTransProb(int id) const;
        /**
        * Returns the emission probability of this state to an output
        * @param id the number of the output
        * @return the emit prob
        */
        T getEmitProb(int id) const;

        /**
        * Sets a transition probability. Does NOT normalize this State's trans probs.
        * @param id the number of the state to set the transition to
        * @param prob the prob to set it to
        */
        void setTransProb(int id, T prob);
        /**
        * Sets an emission probability. Does NOT normalize this State's emit probs.
        * @param id the number of the output to set the emission to
        * @param prob the prob to set it to
        */
        void setEmitProb(int id, T prob);

        /**
        * Returns the number of states in the HMM that this State belongs to, as reported to this State's constructor.
        * @return the number of states
        */
        int getNumStates() const;
        /**
        * Returns the number of outputs of this state, as reported by this State's constructor.
        * @return the number of outputs
        */
        int getNumOutputs() const;

    private:
        int totalStates;
        int totalOutputs;
        T* __restrict__ transProb;
        T* __restrict__ emitProb;
        static mt19937 gen;
        static variate_generator<mt19937, uniform_real<> > rand;
        //get a microsecond-granular value
        static unsigned int getSeed();
    };

	template <typename T> inline T StateTemplate<T>::getTransProb(int id) const
	{
		return transProb[id];
	}

	template <typename T> inline T StateTemplate<T>::getEmitProb(int id) const
	{
		return emitProb[id];
	}

	template <typename T> inline void StateTemplate<T>::setTransProb(int id, T prob)
	{
		transProb[id] = prob;
	}

	template <typename T> inline void StateTemplate<T>::setEmitProb(int id, T prob)
	{
		emitProb[id] = prob;
	}

    template <typename T> inline int StateTemplate<T>::getNumStates() const {return totalStates;}

    template <typename T> inline int StateTemplate<T>::getNumOutputs() const {return totalOutputs;}

	template <typename T> mt19937 StateTemplate<T>::gen( static_cast<unsigned int>(StateTemplate::getSeed()));
	template <typename T> variate_generator<mt19937, uniform_real<> > StateTemplate<T>::rand(gen, uniform_real<double>(0.0, 1.0));

	template <typename T> StateTemplate<T>::StateTemplate(int numStates, int numOutputs)
	: totalStates (numStates)
	, totalOutputs (numOutputs)
	, transProb (new T[numStates])
	, emitProb (new T[numOutputs])
	{}

    template <typename T> StateTemplate<T>::StateTemplate(const StateTemplate<T>& rhs)
    : totalStates (rhs.getNumStates())
    , totalOutputs (rhs.getNumOutputs())
    {
        transProb = new T[totalStates];
        emitProb = new T[totalOutputs];
        for (int i = 0; i < totalStates; i++)
        {
            transProb[i] = rhs.getTransProb(i);
        }
        for (int i = 0; i < totalStates; i++)
        {
            emitProb[i] = rhs.getEmitProb(i);
        }
    }

	template <typename T> StateTemplate<T>::~StateTemplate()
	{
		delete[] transProb;
		delete[] emitProb;
	}

    template <typename T> unsigned int StateTemplate<T>::getSeed()
    {
        timeval curTime;
        gettimeofday(&curTime, NULL);
        return static_cast<unsigned int>(curTime.tv_usec);
    }

	template <typename T> void StateTemplate<T>::initRandom()
	{
		T count = 0;
		for (int i = 0; i < totalStates; i++)
		{
			transProb[i] = rand();
			count += transProb[i];
		}
		for (int i = 0; i < totalStates; i++)
		{
			transProb[i] /= count;
		}
		count = 0;
		for (int i = 0; i < totalOutputs; i++)
		{
			emitProb[i] = rand();
			count += emitProb[i];
		}
		for (int i = 0; i < totalOutputs; i++)
		{
			emitProb[i] /= count;
		}
	}
}
