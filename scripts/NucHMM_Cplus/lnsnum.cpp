#include <cmath>
#include <iostream>

#ifdef LNSNUM_ENABLE_ADD_BIG
#include <mpfr.h>
#endif

#include "lnsnum.hpp"

using std::ceil;
using std::exp;
using std::floor;
using std::log;
using std::ostream;

namespace HMM
{
#ifdef LNSNUM_ENABLE_ADD_BIG
	mpfr_rnd_t Lnsnum::rndMode = mpfr_get_default_rounding_mode();
	mpfr_prec_t Lnsnum::defaultPrec = mpfr_get_default_prec();
#endif
	Lnsnum::Lnsnum()
	: num (0)
	, zero (true)
	{}

	Lnsnum::Lnsnum(long double n)
	{
		if (n == 0.0)
		{
			num = 0;
			zero = true;
		}
		else
		{
			num = log(n);
			zero = false;
		}
	}

	bool Lnsnum::operator==(const Lnsnum &rhs) const
	{
		if (zero != rhs.isZero()) return false;
		else if (zero && rhs.isZero()) return true;
		else if (num == rhs.getLog()) return true;
		else return false;
	}

	bool Lnsnum::operator!=(const Lnsnum &rhs) const
	{
		return !(*this == rhs);
	}

	bool Lnsnum::operator==(const long double &rhs) const
	{
		if (zero != (rhs == 0.0)) return false;
		else if (zero && rhs == 0.0) return true;
		else if (num == log(rhs)) return true;
		else return false;
	}

	bool Lnsnum::operator!=(const long double &rhs) const
	{
		return !(*this == rhs);
	}

	ostream& operator <<(ostream &os, const Lnsnum &rhs)
	{
		os << rhs.getVal();
		return os;
	}

	Lnsnum operator+(long double lhs, Lnsnum &rhs)
	{
		Lnsnum lhsLnsnum(lhs);
		return lhsLnsnum + rhs;
	}

	Lnsnum operator+(Lnsnum &lhs, long double &rhs)
	{
		return (rhs + lhs);
	}

	Lnsnum operator-(long double lhs, Lnsnum &rhs)
	{
		Lnsnum lhsLnsnum(lhs);
		return lhsLnsnum - rhs;
	}

	Lnsnum operator-(Lnsnum &lhs, long double rhs)
	{
		Lnsnum rhsLnsnum(rhs);
		return lhs - rhsLnsnum;
	}

	Lnsnum operator*(long double lhs, Lnsnum &rhs)
	{
		Lnsnum lhsLnsnum(lhs);
		return lhsLnsnum * rhs;
	}

	Lnsnum operator*(Lnsnum &lhs, long double &rhs)
	{
		return (rhs * lhs);
	}

	Lnsnum operator/(long double lhs, Lnsnum &rhs)
	{
		Lnsnum lhsLnsnum(lhs);
		return lhsLnsnum / rhs;
	}

	Lnsnum operator/(Lnsnum &lhs, long double rhs)
	{
		Lnsnum rhsLnsnum(rhs);
		return lhs / rhsLnsnum;
	}

	bool operator<(const Lnsnum &lhs, const Lnsnum &rhs)
	{
		if (rhs.isZero()) return false;
		else if (lhs.isZero() && !rhs.isZero()) return true;
		else return (lhs.getLog() < rhs.getLog());
	}
	bool operator>(const Lnsnum &lhs, const Lnsnum &rhs)
	{
		if (lhs.isZero()) return false;
		else if (!lhs.isZero() && rhs.isZero()) return true;
		else return (lhs.getLog() > rhs.getLog());
	}
	bool operator<=(const Lnsnum &lhs, const Lnsnum &rhs)
	{
		if (lhs.isZero() && !rhs.isZero()) return true;
		else if (!lhs.isZero() && rhs.isZero()) return false;
		else if (lhs.isZero() && rhs.isZero()) return true;
		else return (lhs.getLog() <= rhs.getLog());
	}
	bool operator>=(const Lnsnum &lhs, const Lnsnum &rhs)
	{
		if (lhs.isZero() && !rhs.isZero()) return false;
		else if (!lhs.isZero() && rhs.isZero()) return true;
		else if (lhs.isZero() && rhs.isZero()) return true;
		else return (lhs.getLog() >= rhs.getLog());
	}
	bool operator<(const Lnsnum &lhs, const long double &rhs)
	{
		if (lhs.isZero()) return (0 < rhs);
		else return (lhs.getLog() < log(rhs));
	}
	bool operator>(const Lnsnum &lhs, const long double &rhs)
	{
		if (lhs.isZero()) return (0 > rhs);
		else return (lhs.getLog() > log(rhs));
	}
	bool operator<=(const Lnsnum &lhs, const long double &rhs)
	{
		if (lhs.isZero()) return (0 <= rhs);
		else return (lhs.getLog() <= log(rhs));
	}
	bool operator>=(const Lnsnum &lhs, const long double &rhs)
	{
		if (lhs.isZero()) return (0 >= rhs);
		else return (lhs.getLog() >= log(rhs));
	}
	bool operator<(const long double &lhs, const Lnsnum &rhs)
		{return (rhs > lhs);}
	bool operator>(const long double &lhs, const Lnsnum &rhs)
		{return (rhs < lhs);}
	bool operator<=(const long double &lhs, const Lnsnum &rhs)
		{return (rhs >= lhs);}
	bool operator>=(const long double &lhs, const Lnsnum &rhs)
		{return (rhs <= lhs);}

	Lnsnum Lnsnum::pow(const Lnsnum &base, const Lnsnum &exponent)
	{
		if (exponent.isZero())
		{
			Lnsnum ans(1);
			return ans;
		}
		else if (base.isZero())
		{
			Lnsnum ans(0);
			return ans;
		}
		else if (base.getLog() < 0.0)
		{
			long double temp = log(-base.getLog());
			temp += exponent.getLog();
			long double temp2 = -exp(temp);
			Lnsnum ans(1);
			ans.setLog(temp2);
			return ans;
		}
		else
		{
			long double temp = log(base.getLog());
			temp += exponent.getLog();
			long double temp2 = exp(temp);
			Lnsnum ans(1);
			ans.setLog(temp2);
			return ans;
		}
	}

#ifdef LNSNUM_ENABLE_ADD_BIG
	Lnsnum Lnsnum::addBig(const Lnsnum& a, const Lnsnum& b)
	{
		if (b.isZero()) {return a;}
		else if (a.isZero()) {return b;}
		else
		{
			long double result;
			mpfr_t l;
			mpfr_t r;
			mpfr_t dif;
			mpfr_t ans;
			mpfr_t log;
			mpfr_inits2(defaultPrec, l, r, dif, ans, log, (mpfr_ptr) NULL);
			mpfr_set_d(l, a.getLog(), rndMode);
			mpfr_set_d(r, b.getLog(), rndMode);

			mpfr_sub(dif, l, r, rndMode);
			mpfr_exp(ans, dif, rndMode);
			mpfr_add_ui(dif, ans, (unsigned long int) 1, rndMode);
			mpfr_log(log, dif, rndMode);
			mpfr_add(ans, log, r, rndMode);
			result = mpfr_get_d(ans, rndMode);
			mpfr_clears(l, r, dif, ans, log, (mpfr_ptr) NULL);
			Lnsnum res(1);
			res.setLog(result);
			return res;
		}
	}
#endif

	void Lnsnum::setLog(long double n)
	{
		num = n;
	}

	Lnsnum Lnsnum::balancedAdd(const Lnsnum* list, const int listSize)
	{
		if (listSize == 0) return 0;
		else return Lnsnum::balancedAddHelper(list, 0, listSize - 1);
	}

	Lnsnum Lnsnum::balancedAddHelper(const Lnsnum* list, const int begin, const int end)
	{
		if (begin == end) return list[begin];
		else
		{
			int midpoint = (end + begin) / 2;
			return Lnsnum::balancedAddHelper(list, begin, midpoint) + Lnsnum::balancedAddHelper(list, midpoint + 1, end);
		}
	}

} // namespace HMM
