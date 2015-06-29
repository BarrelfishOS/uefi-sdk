//===-- lib/truncdfsf2.c - quad -> double conversion --------------*- C -*-===//
//
//                     The LLVM Compiler Infrastructure
//
// This file is dual licensed under the MIT and the University of Illinois Open
// Source Licenses. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#include <Base.h>

#include "../Arm/Llvm_int_lib.h"

typedef long double src_t;
typedef __uint128_t src_rep_t;
#define SRC_REP_C (__uint128_t)
static const int srcSigBits = 112;

typedef double dst_t;
typedef UINT64 dst_rep_t;
#define DST_REP_C (UINT64)
static const int dstSigBits = 52;

static inline src_rep_t srcToRep(src_t x) {
    const union { src_t f; src_rep_t i; } rep = {.f = x};
    return rep.i;
}

static inline dst_t dstFromRep(dst_rep_t x) {
    const union { dst_t f; dst_rep_t i; } rep = {.i = x};
    return rep.f;
}

static inline dst_t __truncXfYf2__(src_t a) {
    // Various constants whose values follow from the type parameters.
    // Any reasonable optimizer will fold and propagate all of these.
    const int srcBits = sizeof(src_t)*CHAR_BIT;
    const int srcExpBits = srcBits - srcSigBits - 1;
    const int srcInfExp = (1 << srcExpBits) - 1;
    const int srcExpBias = srcInfExp >> 1;

    const src_rep_t srcMinNormal = SRC_REP_C(1) << srcSigBits;
    const src_rep_t srcSignificandMask = srcMinNormal - 1;
    const src_rep_t srcInfinity = (src_rep_t)srcInfExp << srcSigBits;
    const src_rep_t srcSignMask = SRC_REP_C(1) << (srcSigBits + srcExpBits);
    const src_rep_t srcAbsMask = srcSignMask - 1;
    const src_rep_t roundMask = (SRC_REP_C(1) << (srcSigBits - dstSigBits)) - 1;
    const src_rep_t halfway = SRC_REP_C(1) << (srcSigBits - dstSigBits - 1);
    const src_rep_t srcQNaN = SRC_REP_C(1) << (srcSigBits - 1);
    const src_rep_t srcNaNCode = srcQNaN - 1;

    const int dstBits = sizeof(dst_t)*CHAR_BIT;
    const int dstExpBits = dstBits - dstSigBits - 1;
    const int dstInfExp = (1 << dstExpBits) - 1;
    const int dstExpBias = dstInfExp >> 1;

    const int underflowExponent = srcExpBias + 1 - dstExpBias;
    const int overflowExponent = srcExpBias + dstInfExp - dstExpBias;
    const src_rep_t underflow = (src_rep_t)underflowExponent << srcSigBits;
    const src_rep_t overflow = (src_rep_t)overflowExponent << srcSigBits;

    const dst_rep_t dstQNaN = DST_REP_C(1) << (dstSigBits - 1);
    const dst_rep_t dstNaNCode = dstQNaN - 1;

    // Break a into a sign and representation of the absolute value
    const src_rep_t aRep = srcToRep(a);
    const src_rep_t aAbs = aRep & srcAbsMask;
    const src_rep_t sign = aRep & srcSignMask;
    dst_rep_t absResult;

    if (aAbs - underflow < aAbs - overflow) {
        // The exponent of a is within the range of normal numbers in the
        // destination format.  We can convert by simply right-shifting with
        // rounding and adjusting the exponent.
        absResult = aAbs >> (srcSigBits - dstSigBits);
        absResult -= (dst_rep_t)(srcExpBias - dstExpBias) << dstSigBits;

        const src_rep_t roundBits = aAbs & roundMask;
        // Round to nearest
        if (roundBits > halfway)
            absResult++;
        // Ties to even
        else if (roundBits == halfway)
            absResult += absResult & 1;
    }
    else if (aAbs > srcInfinity) {
        // a is NaN.
        // Conjure the result by beginning with infinity, setting the qNaN
        // bit and inserting the (truncated) trailing NaN field.
        absResult = (dst_rep_t)dstInfExp << dstSigBits;
        absResult |= dstQNaN;
        absResult |= ((aAbs & srcNaNCode) >> (srcSigBits - dstSigBits)) & dstNaNCode;
    }
    else if (aAbs >= overflow) {
        // a overflows to infinity.
        absResult = (dst_rep_t)dstInfExp << dstSigBits;
    }
    else {
        // a underflows on conversion to the destination type or is an exact
        // zero.  The result may be a denormal or zero.  Extract the exponent
        // to get the shift amount for the denormalization.
        const int aExp = aAbs >> srcSigBits;
        const int shift = srcExpBias - dstExpBias - aExp + 1;

        const src_rep_t significand = (aRep & srcSignificandMask) | srcMinNormal;

        // Right shift by the denormalization amount with sticky.
        if (shift > srcSigBits) {
            absResult = 0;
        } else {
            const BOOLEAN sticky = significand << (srcBits - shift);
            src_rep_t denormalizedSignificand = significand >> shift | sticky;
            absResult = denormalizedSignificand >> (srcSigBits - dstSigBits);
            const src_rep_t roundBits = denormalizedSignificand & roundMask;
            // Round to nearest
            if (roundBits > halfway)
                absResult++;
            // Ties to even
            else if (roundBits == halfway)
                absResult += absResult & 1;
        }
    }

    // Apply the signbit to (dst_t)abs(a).
    const dst_rep_t result = absResult | sign >> (srcBits - dstBits);
    return dstFromRep(result);
}

double __trunctfdf2(long double a) {
    return __truncXfYf2__(a);
}
