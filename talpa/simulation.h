// ---------------------------------------------------------
//
//  simulation.h
//  Tyson Brochu 2008
//
//  Time keeping and running status.
//
// ---------------------------------------------------------

#ifndef EL_TOPO_SIMULATION_H
#define EL_TOPO_SIMULATION_H

// ---------------------------------------------------------
// Nested includes
// ---------------------------------------------------------

#include <options.h>

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

// ---------------------------------------------------------
//  Class definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
///  Simulation parameters
///
// ---------------------------------------------------------

class Simulation
{
public:
    
    /// simulation time step size
    ///
    double m_dt;
    
    /// end time
    ///
    double m_max_t;
    
    /// current simulation time
    ///
    double m_curr_t;
    
    /// whether we're currently running an entire simulation
    ///
    bool m_running;
    
    /// whether we're currently running a single timestep of the simulation
    ///
    bool m_currently_advancing_simulation;
    
    
    /// Specify a target time step and optional end time
    ///
    Simulation( double in_dt, double in_max_t = BIG_DOUBLE ) :
    m_dt(in_dt),
    m_max_t(in_max_t),
    m_curr_t(0),
    m_running(false),
    m_currently_advancing_simulation(false)
    {
        assert( m_dt > 0 );
        assert( m_max_t > 0 );
    }
    
    bool done_simulation( )
    {
        return m_curr_t > m_max_t - 1e-7;
    }
    
    
    
};


#endif



