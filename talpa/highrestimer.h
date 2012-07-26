//
//  highrestimer.h
//  eltopo3d_project
//
//  Created by Tyson Brochu on 11-12-06.
//

#ifndef eltopo3d_project_highrestimer_h
#define eltopo3d_project_highrestimer_h

#include <mach/clock.h>
#include <mach/mach_time.h>
#include <mach/mach.h>
#include <CoreServices/CoreServices.h>


inline uint64_t get_absolute_time()
{
    return mach_absolute_time();
}


inline uint64_t convert_interval_to_nanoseconds( uint64_t time_interval )
{
    Nanoseconds elapsed_nano = AbsoluteToNanoseconds( *(AbsoluteTime *) &time_interval ); 
    return * (uint64_t *) &elapsed_nano;
}

inline double convert_interval_to_seconds( uint64_t time_interval )
{
    Nanoseconds elapsed_nano = AbsoluteToNanoseconds( *(AbsoluteTime *) &time_interval ); 
    return 1e-9 * double(*(uint64_t*)&elapsed_nano);
}

#endif
