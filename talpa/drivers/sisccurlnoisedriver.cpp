// ---------------------------------------------------------
//
//  sisccurlnoisedriver.cpp
//  Tyson Brochu 2008
//
//  Mesh driver for divergence-free, pseuo-random velocity
//
// ---------------------------------------------------------

// ---------------------------------------------------------
// Includes
// ---------------------------------------------------------

#include <sisccurlnoisedriver.h>

#include <runstats.h>
#include <surftrack.h>

// ---------------------------------------------------------
// Global externs
// ---------------------------------------------------------

// ---------------------------------------------------------
// Local constants, typedefs, macros
// ---------------------------------------------------------

// ---------------------------------------------------------
// Static function definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
// Member function definitions
// ---------------------------------------------------------


// ---------------------------------------------------------
///
/// Constructor
///
// ---------------------------------------------------------

SISCCurlNoiseDriver::SISCCurlNoiseDriver() : 
   CurlNoise3( ),
   noise_lengthscale(1),
   noise_gain(1),
   noise(),
   m_initial_volume(UNINITIALIZED_DOUBLE)
{
    noise_lengthscale[0]=1.5;
    noise_gain[0]=1.3;  
}


// ---------------------------------------------------------
///
/// Initializer.  Get the initial volume enclosed by the surface for comparison later on.
///
// ---------------------------------------------------------

void SISCCurlNoiseDriver::initialize( const SurfTrack& surf )
{
    m_initial_volume = surf.get_volume();
}

// ---------------------------------------------------------
///
/// 3D vector field which defines the velocity via its curl.
///
// ---------------------------------------------------------

Vec3d SISCCurlNoiseDriver::potential(double x, double y, double z) const
{
    Vec3d psi(0,0,0);
    double height_factor=0.5;
    
    static const Vec3d centre( 0.0, 1.0, 0.0 );
    static double radius = 4.0;
    
    for(unsigned int i=0; i<noise_lengthscale.size(); ++i)
    {
        double sx=x/noise_lengthscale[i];
        double sy=y/noise_lengthscale[i];
        double sz=z/noise_lengthscale[i];
        
        Vec3d psi_i( 0.f, 0.f, noise2(sx,sy,sz));
        
        double dist = mag( Vec3d(x,y,z) - centre );      
        double scale = max( (radius - dist)/radius, 0.0 );
        psi_i *= scale;
        
        psi+=height_factor*noise_gain[i]*psi_i;
    }
    
    return psi;
}


// ---------------------------------------------------------
///
/// For each vertex on the mesh, take the curl of potential to get velocity vector at that vertex location.
/// Uses RK4 to get x(t+1) from x(t), then returns linear trajectory v = (x(t+1) - x(t)) / dt.
///
// ---------------------------------------------------------

void SISCCurlNoiseDriver::set_predicted_vertex_positions( const SurfTrack& surf, std::vector<Vec3d>& predicted_positions, double /*current_t*/, double& adaptive_dt )
{
    const std::vector<Vec3d>& positions = surf.get_positions();
    Vec3d v, midx;
    
    // Update velocities
    for ( unsigned int i = 0; i < positions.size(); ++i )
    {    
        
        const Vec3d& p = positions[i];
        
        if ( surf.m_mesh.m_vertex_to_triangle_map[i].empty() )
        {
            predicted_positions[i] = p;
            continue;
        }
        
        // RK4
        // -----------
        // k1 = dt * f( t, x );
        get_velocity( p, v );
        Vec3d k1 = adaptive_dt * v;
        
        // k2 = dt * f( t + 0.5*dt, x + 0.5*k1 );
        advance_time( 0.5*adaptive_dt );
        get_velocity( p + 0.5 * k1, v );
        Vec3d k2 = adaptive_dt * v;
        
        //k3 = dt * f( t + 0.5*dt, x + 0.5*k2 );
        get_velocity( p + 0.5 * k2, v );
        Vec3d k3 = adaptive_dt * v;
        
        //k4 = dt * f( t + dt, x + k3 );
        advance_time( 0.5*adaptive_dt );
        get_velocity( p + 0.5 * k3, v );
        Vec3d k4 = adaptive_dt * v;
        
        predicted_positions[i] = p + 1.0/6.0 * ( k1 + k4 ) + 1.0/3.0 * ( k2 + k3 );
        
        advance_time( -adaptive_dt );
    }
    
    advance_time( adaptive_dt );
}


// ---------------------------------------------------------
///
/// Compute and output difference between current volume and initial volume
///
// ---------------------------------------------------------

void SISCCurlNoiseDriver::compute_error( const SurfTrack& surf, double current_t )
{
    double delta_volume = surf.get_volume() - m_initial_volume;
    double percent_error = delta_volume / m_initial_volume * 100;
    
    static unsigned int curr_frame = 0;
    extern RunStats g_stats;  
    g_stats.set_double( "last_t", current_t );   
    g_stats.set_double( "last_delta_volume", delta_volume );
    g_stats.set_double( "last_percent_error", percent_error );
    ++curr_frame;
    
}


