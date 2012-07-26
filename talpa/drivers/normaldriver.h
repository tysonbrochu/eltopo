// ---------------------------------------------------------
//
//  normaldriver.h
//  Tyson Brochu 2008
//
//  Mesh driver for motion in the normal direction using vertex normals.
//  NOTE: This was implemented only for comparison against FaceOff.  Not recommended for general use.
//
// ---------------------------------------------------------

#ifndef EL_TOPO_NORMALDRIVER_H
#define EL_TOPO_NORMALDRIVER_H

// ---------------------------------------------------------
// Nested includes
// ---------------------------------------------------------

#include <meshdriver.h>
#include <mat.h>
#include <vec.h> 
#include <vector>

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

class NonDestructiveTriMesh;
class SurfTrack;

// ---------------------------------------------------------
//  Class definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Motion in the normal direction using vertex normals 
/// (area-weighted average of incident triangle normals).
///
/// This driver starts with two spheres, moves in the positive
/// normal direction, then in the negative normal direction.
/// Error computations are done vs. an analytical solution.
///
// ---------------------------------------------------------

class NormalDriver : public MeshDriver
{
public:
    
    NormalDriver( double in_speed, 
                 const Vec3d& in_sphere_a_centre,  
                 const Vec3d& in_sphere_b_centre, 
                 double in_max_radius, 
                 double in_interior_radius ) : 
    speed(in_speed),
    sphere_a_centre(in_sphere_a_centre),
    sphere_b_centre(in_sphere_b_centre),   
    max_radius(in_max_radius),
    interior_radius(in_interior_radius)
    {}
    
    void initialize( const SurfTrack& ) {}
    
    // Assign a velocity vector to each mesh vertex
    // 
    void set_predicted_vertex_positions( const SurfTrack& surf, std::vector<Vec3d>& predicted_positions, double current_t, double& adaptive_dt );
    
    /// Compute the distance from vertices to analytic surface, weighted by associated vertex area
    /// 
    double compute_l1_error( const SurfTrack& surf );
    
    /// Compute the maximum distance from vertices to analytic surface
    /// 
    double compute_inf_error( const SurfTrack& surf );
    
    /// Compute and output both L1 and L_inf errors
    ///
    void compute_error( const SurfTrack& surf, double current_t );
    
    /// Speed of normal motion
    double speed;
    
    //
    // The following members are used to compute the analytic solution
    //
    
    /// Geometry to compare against
    Vec3d sphere_a_centre, sphere_b_centre;
    
    /// radius of spheres when motion switches from positive to negative
    double max_radius;         
    
    /// difference between maximum radius and final radius
    double interior_radius;    
    
};


// ---------------------------------------------------------
//  Inline functions
// ---------------------------------------------------------

#endif


