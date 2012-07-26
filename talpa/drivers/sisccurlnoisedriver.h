// ---------------------------------------------------------
//
//  sisccurlnoisedriver.h
//  Tyson Brochu 2008
//
//  Mesh driver for divergence-free, pseuo-random velocity
//
// ---------------------------------------------------------

#ifndef EL_TOPO_SISCCURLNOISEDRIVER_H
#define EL_TOPO_SISCCURLNOISEDRIVER_H

// ---------------------------------------------------------
// Nested includes
// ---------------------------------------------------------

#include <curlnoise.h>
#include <meshdriver.h>
#include <noise.h>

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

class SurfTrack;

// ---------------------------------------------------------
//  Class definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Driver for the curl noise example in the SISC paper.
///
// ---------------------------------------------------------

struct SISCCurlNoiseDriver: public CurlNoise3, public MeshDriver
{
    
    /// Default constructor
    ///
    SISCCurlNoiseDriver(void);
    
    /// Initializer
    ///
    void initialize( const SurfTrack& surf );
    
    //
    // Virtual functions 
    //
    
    virtual bool seed_particles(std::vector<Vec3d>& /*x*/, double /*dt*/ ) const { return true; }   
    void advance_time(double /*dt*/) {}
    
    /// z-component of noise vector (x and y are 0)
    ///
    double noise2(double x, double y, double z) const { return noise(z-203.994, x+169.47, y-205.31); }
    
    /// Take the curl of this function to get velocity
    ///
    Vec3d potential(double x, double y, double z) const;
    
    /// Set velocity on each mesh vertex
    ///
    virtual void set_predicted_vertex_positions( const SurfTrack& surf, std::vector<Vec3d>& predicted_positions, double current_t, double& adaptive_dt );
    
    /// Compute and output error (difference between current volume and initial volume)
    ///
    void compute_error( const SurfTrack& surf, double current_t );
    
    /// Parameters for each scale of noise
    ///
    std::vector<double> noise_lengthscale, noise_gain;
    
    /// Noise generator
    /// 
    FlowNoise3 noise;
    
    /// Volume enclosed by surface at t = 0.  Compare to volume at current time to get error.
    ///
    double m_initial_volume;
    
};

// ---------------------------------------------------------
//  Inline functions
// ---------------------------------------------------------


#endif
