#pragma once

#if __GNUC__  >= 3
# define _noinline     __attribute_noinline__
# define _pure         __attribute_pure__
# define _const        __attribute_const__
# define _noreturn     __attribute_noreturn__
# define _malloc       __attribute_malloc__
# define _must_check   __attribute_warn_unused_result__
# define _deprecated   __attribute_deprecated__
# define _used         __attribute_used__
# define _align(x)     __attribute_aligned__ (x)
# define _align_max    __attribute_aligned__
# define branch_prediction_unlikely(expr) __builtin_expect(!!(expr), 0)
# define branch_prediction_likely(expr) __builtin_expect(!!(expr), 1)
#else
# define _pure         /* no pure */
# define _const        /* no const */
# define _noreturn     /* no noreturn */
# define _malloc       /* no malloc */
# define _must_check   /* no warn_unused_result */
# define _deprecated   /* no deprecated */
# define _used         /* no used */
# define _align(x)     /* no aligned */
# define _align_max    /* no align_max */
# define branch_prediction_unlikely(expr) (expr)
# define branch_prediction_likely(expr) (expr)
#endif
//http://lwn.net/Articles/255364/
