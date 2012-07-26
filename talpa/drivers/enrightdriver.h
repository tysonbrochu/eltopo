// ---------------------------------------------------------
//
//  enright.h
//  Tyson Brochu 2008
//
//  Mesh driver for motion specified by an analytic velocity field, using RK4 time integration.
//
// ---------------------------------------------------------

#ifndef EL_TOPO_ENRIGHT_DRIVER_H
#define EL_TOPO_ENRIGHT_DRIVER_H

// ---------------------------------------------------------
// Nested includes
// ---------------------------------------------------------

#include <cmath>
#include <meshdriver.h>
#include <runstats.h>
#include <surftrack.h>

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

// ---------------------------------------------------------
//  Class definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Driver for the Enright test.
///
// ---------------------------------------------------------

class EnrightDriver : public MeshDriver
{
    
public:
    
    EnrightDriver() : m_initial_volume(UNINITIALIZED_DOUBLE) {}
    
    inline void initialize( const SurfTrack& );
    
    // set velocity at (t,x)
    inline void enright_velocity( double t, const Vec3d& pos, Vec3d& out );
    
    inline virtual void set_predicted_vertex_positions( const SurfTrack& surf, std::vector<Vec3d>& predicted_positions, double current_t, double& adaptive_dt );
    
    inline void compute_error( const SurfTrack& surf, double current_t );
    
    // Volume enclosed by surface at t = 0.  Compare to end time to get error.
    double m_initial_volume;
    
};

// ---------------------------------------------------------
//  Inline functions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Initialization.  Just grab the initial volume for error computation later.
///
// ---------------------------------------------------------

inline void EnrightDriver::initialize( const SurfTrack& surf )
{
    m_initial_volume = surf.get_volume();
}


// ---------------------------------------------------------
///
/// Set velocity at location pos, time t.
///
// ---------------------------------------------------------

inline void EnrightDriver::enright_velocity( double t, const Vec3d& pos, Vec3d& out )
{
    double x = pos[0]; 
    double y = pos[1]; 
    double z = pos[2];
    
    out = Vec3d( 2.0 * std::sin(M_PI*x) * std::sin(M_PI*x) * std::sin(2.0*M_PI*y) * std::sin(2.0*M_PI*z),
                -std::sin(2.0*M_PI*x) * std::sin(M_PI*y)*std::sin(M_PI*y) * std::sin(2.0*M_PI*z),
                -std::sin(2.0*M_PI*x) * std::sin(2.0*M_PI*y) * std::sin(M_PI*z) * std::sin(M_PI*z) );
    
    out *= sin(M_PI * t * 2 / 3);    // modulate with a period of 3
    
}


// ---------------------------------------------------------
///
/// Set velocity vectors on all vertices according to enright_velocity.  
/// Uses RK4 to get x(t+1) from x(t), then returns v = (x(t+1) - x(t)) / dt.
///
// ---------------------------------------------------------

inline void EnrightDriver::set_predicted_vertex_positions( const SurfTrack& surf, std::vector<Vec3d>& predicted_positions, double current_t, double& dt )
{   
    const std::vector<Vec3d>& x = surf.get_positions();  
    
    predicted_positions.resize( x.size() );
    
    for(unsigned int i=0; i < x.size(); ++i) 
    {
        Vec3d v;
        
        // RK4
        // -----------
        // k1 = dt * f( t, x );
        enright_velocity( current_t, x[i], v );
        Vec3d k1 = dt * v;
        
        // k2 = dt * f( t + 0.5*dt, x + 0.5*k1 );
        enright_velocity( current_t + 0.5*dt, x[i] + 0.5 * k1, v );
        Vec3d k2 = dt * v;
        
        // k3 = dt * f( t + 0.5*dt, x + 0.5*k2 );
        enright_velocity( current_t + 0.5*dt, x[i] + 0.5 * k2, v );
        Vec3d k3 = dt * v;
        
        // k4 = dt * f( t + dt, x + k3 );
        enright_velocity( current_t + dt, x[i] + k3, v );
        Vec3d k4 = dt * v;
        
        predicted_positions[i] = x[i] + 1./6. * ( k1 + k4 ) + 1./3. * ( k2 + k3 );      
    }
    
}


// ---------------------------------------------------------
///
/// Compute and output change in volume so far.
///
// ---------------------------------------------------------

inline void EnrightDriver::compute_error( const SurfTrack& surf, double current_t )
{
    double delta_volume = surf.get_volume() - m_initial_volume;
    double percent_error = delta_volume / m_initial_volume * 100;
    
    static unsigned int curr_frame = 0;
    extern RunStats g_stats;
    g_stats.set_double( "last_t", current_t );   
    g_stats.set_double( "delta_volume", delta_volume );
    g_stats.set_double( "percent_error", percent_error );
    ++curr_frame;
}


#endif
