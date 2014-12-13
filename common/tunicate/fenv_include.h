#ifndef FENVINCLUDE_H
#define FENVINCLUDE_H

#ifndef _MSC_VER

//if not Windows, include the standard fenv.v.h
#include <fenv.h>

#else 

//on Windows, cook up the functions we need
#include <float.h>
#include <stdio.h>
#include <errno.h>

#pragma fenv_access(on)


//define the functions we need. I hate Microsoft! Just implement C99 already, c'mon.
#define FE_DOWNWARD   _RC_DOWN
#define FE_UPWARD     _RC_UP
#define FE_TONEAREST  _RC_NEAR
#define FE_TOWARDZERO _RC_CHOP

inline void fesetround(unsigned int choice) {
  _controlfp(choice, _MCW_RC);
}

inline int fegetround() {
  unsigned int control_word = _controlfp(0,0);
  int roundmode = (control_word & _MCW_RC);
  return roundmode;
}

#endif

#endif
