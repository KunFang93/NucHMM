#pragma once

// requires mpfr.h
// #define LNSNUM_ENABLE_ADD_BIG

#include <cmath>
#include <iostream>

#ifdef LNSNUM_ENABLE_ADD_BIG
#include <mpfr.h>
#endif

using std::exp;
using std::log;
using std::ostream;

namespace HMM
{
    /**
    * A number that stores the natural log of its value (and a zero bit). Intented for use with extremely large or extremely small numbers, or in cases where fast division is
    * more important than high precision (division is implemented as subtraction of logs). Precision decreases with increasing value distance from e.
    * @warning Note that addition and subtraction can fail if the magnitudes of the numbers are significantly different (due to over/underflow of intermediate numbers)!
    * Use Lnsnum::addBig if this becomes a problem (though it is slow and requires MPFR), or Lnsnum::balancedAdd for lists of Lnsnum (though this uses recursion
    * of depth lg(n) and can potentially blow the stack)
    */
	class Lnsnum
	{
		public:
			Lnsnum();
			Lnsnum(long double n);

			Lnsnum& operator=(const Lnsnum& rhs);
			Lnsnum& operator+=(const Lnsnum& rhs);
			Lnsnum& operator-=(const Lnsnum& rhs);
			Lnsnum& operator*=(const Lnsnum& rhs);
			Lnsnum& operator/=(const Lnsnum& rhs);
			const Lnsnum operator+(const Lnsnum& rhs) const;
			const Lnsnum operator-(const Lnsnum& rhs) const;
			const Lnsnum operator*(const Lnsnum& rhs) const;
			const Lnsnum operator/(const Lnsnum& rhs) const;
			Lnsnum& operator=(const long double& rhs);

            //these comparisons deliberately do not check if |a-b|<epsilon - that's left up to the user
			bool operator==(const Lnsnum& rhs) const;
			bool operator!=(const Lnsnum& rhs) const;
			bool operator==(const long double& rhs) const;
			bool operator!=(const long double& rhs) const;
			friend bool operator<(const Lnsnum& lhs, const Lnsnum& rhs);
			friend bool operator>(const Lnsnum& lhs, const Lnsnum& rhs);
			friend bool operator<=(const Lnsnum& lhs, const Lnsnum& rhs);
			friend bool operator>=(const Lnsnum& lhs, const Lnsnum& rhs);

			friend ostream& operator <<(ostream& os, const Lnsnum& rhs);

			static Lnsnum pow(const Lnsnum& base, const Lnsnum& exponent);
			/**
			* Adds a list of Lnsnum.
			* @warning this uses recursion of depth lg(listSize), can potentially blow the stack with very long lists
			* @warning the Lnsnum in the list should be of similar magnitude, else the addition can potentially fail as intermediate numbers over/underflow
			* @param list the list of Lnsnum
			* @param listSize the number of Lnsnum in the list
			* @return the sum
			*/
			static Lnsnum balancedAdd(const Lnsnum* list, const int listSize);

#ifdef LNSNUM_ENABLE_ADD_BIG
            /**
            * Adds two Lnsnum of significantly different magnitudes. Much slower than operator+ and requires MPFR, but avoids most if not all intermediate over/underflow issues.
            * @param a one Lnsnum
            * @param b another Lnsnum
            * @return their sum
            */
			static Lnsnum addBig(const Lnsnum& a, const Lnsnum& b);
#endif
            /**
            * Returns the natural log of this Lnsnum's value (actually, the internal log representation of the value)
            * @return log(this)
            */
			long double getLog() const;
			/**
			* Returns the value of this Lnsnum as a long double.
			* @warning this can easily over/underflow, check the log (see getLog()) against log(LDBL_MAX) and log(LDBL_MIN) to guard against this
			*/
			long double getVal() const;
			/**
			* Returns true if the value of this Lnsnum is zero, otherwise returns false
			* @return (*this == 0)
			*/
			bool isZero() const;
			/**
			* Sets the Lnsnum to e^n if it is nonzero.
			* @warning if this is zero, setLog will have no effect! If you are unsure, you can use Lnsnum l = 1 before l.setLog(n).
			* @param n the natural log of the new value
			*/
			void setLog(long double n);
		private:
			long double num;
			bool zero;
			static Lnsnum balancedAddHelper(const Lnsnum* list, const int begin, const int end);

#ifdef LNSNUM_ENABLE_ADD_BIG
			static mpfr_rnd_t rndMode;
			static mpfr_prec_t defaultPrec;
#endif
	}; //class Lnsnum

	inline Lnsnum& Lnsnum::operator+=(const Lnsnum &rhs)
	{
		if (rhs.isZero()) {}
		else if (zero)
		{
			num = rhs.getLog();
			zero = false;
		}
		else
		{
			num = rhs.getLog() + log(exp(num - rhs.getLog()) + 1);
		}
		return *this;
	}

	inline Lnsnum& Lnsnum::operator-=(const Lnsnum &rhs)
	{
		if (rhs.isZero()) {}
		else if (zero) {/* implement if negative numbers are necessary */}
		else if (num == rhs.getLog())
		{
			num = 0.0;
			zero = true;
		}
		else num = rhs.getLog() + log(exp(num - rhs.getLog()) - 1);
		return *this;
	}

	inline Lnsnum& Lnsnum::operator*=(const Lnsnum &rhs)
	{
		if (zero || rhs.isZero())
		{
			zero = true;
			num = 0.0;
		}
		else num += rhs.getLog();
		return *this;
	}

	inline Lnsnum& Lnsnum::operator/=(const Lnsnum &rhs)
	{
		if (zero || rhs == 0.0)
		{
			zero = true;
			num = 0.0;
		}
		else num -= rhs.getLog();
		return *this;
	}

	inline const Lnsnum Lnsnum::operator+(const Lnsnum &rhs) const
	{
		Lnsnum result = *this;
		result += rhs;
		return result;
	}

	inline const Lnsnum Lnsnum::operator-(const Lnsnum &rhs) const
	{
		Lnsnum result = *this;
		result -= rhs;
		return result;
	}

	inline const Lnsnum Lnsnum::operator*(const Lnsnum &rhs) const
	{
		Lnsnum result = *this;
		result *= rhs;
		return result;
	}

	inline const Lnsnum Lnsnum::operator/(const Lnsnum &rhs) const
	{
		Lnsnum result = *this;
		result /= rhs;
		return result;
	}

	inline Lnsnum& Lnsnum::operator=(const Lnsnum &rhs)
	{
		if (this != &rhs)
		{
			if (rhs.isZero()) zero = true;
			else zero = false;
			num = rhs.getLog();
		}
		return *this;
	}

	inline Lnsnum& Lnsnum::operator=(const long double &rhs)
	{
		if (rhs == 0.0)
		{
			zero = true;
			num = 0.0;
		}
		else
		{
			zero = false;
			num = log(rhs);
		}
		return *this;
	}
	inline long double Lnsnum::getLog() const {return num;}

	inline long double Lnsnum::getVal() const
	{
		if (zero) return 0.0;
		else return exp(num);
	}

	inline bool Lnsnum::isZero() const {return zero;}

} // namespace HMM
