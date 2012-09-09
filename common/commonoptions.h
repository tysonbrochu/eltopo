// ---------------------------------------------------------
//
//  options.h
//  Tyson Brochu 2008
//
//  Constants and macro defines
//
// ---------------------------------------------------------

#ifndef COMMONOPTIONS_H
#define COMMONOPTIONS_H

#include <cstddef>

// ---------------------------------------------------------
// Global constants
// ---------------------------------------------------------

// Types

#ifdef _MSC_VER
typedef ptrdiff_t ptrdiff_t;
#endif

const double UNINITIALIZED_DOUBLE = 0x0F;
const double BIG_DOUBLE = 1e30;
const size_t UNINITIALIZED_SIZE_T = static_cast<size_t> (~0);

#endif

