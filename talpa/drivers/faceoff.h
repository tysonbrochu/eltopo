// ---------------------------------------------------------
//
//  faceoff.h
//  Tyson Brochu 2008
//
//  Mesh driver for motion in the normal direction using the faceoff method (entropy solution).
//
// ---------------------------------------------------------

#ifndef EL_TOPO_FACEOFF_H
#define EL_TOPO_FACEOFF_H

// ---------------------------------------------------------
// Nested includes
// ---------------------------------------------------------

#include <meshdriver.h>
#include <mat.h>
#include <vec.h> 
#include <vector>

// TEMP
#include <meshoperator.h>

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

class NonDestructiveTriMesh;

// ---------------------------------------------------------
//  Class definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Motion in the normal direction using Face Offsetting [Jiao 2007].
///
/// Specifically, this driver starts with two spheres, moves in the 
/// positive normal direction, then in the negative normal direction.
/// Error computations are done vs. an analytical solution to this problem.
///
// ---------------------------------------------------------

class FaceOffDriver : public MeshDriver
{
public:
    
    /// Constructor 
    /// 
    FaceOffDriver( double in_speed, const Vec3d& in_sphere_a_centre, const Vec3d& in_sphere_b_centre, double in_max_radius, double in_interior_radius ) : 
    speed(in_speed) ,
    sphere_a_centre(in_sphere_a_centre),
    sphere_b_centre(in_sphere_b_centre),   
    max_radius(in_max_radius),
    interior_radius(in_interior_radius)
    {}
    
    void initialize( const SurfTrack& );
    
    /// Get the quadric metric tensor at a vertex from the given incident triangles
    ///   
    void compute_quadric_metric_tensor( const std::vector<Vec3d>& triangle_normals, 
                                       const std::vector<double>& triangle_areas, 
                                       const std::vector<size_t>& incident_triangles,
                                       Mat33d& quadric_metric_tensor );
    
    /// Return intersection point between a set of planes in the least-squares sense
    ///
    void intersection_point( const std::vector<Vec3d>& triangle_normals, 
                            const std::vector<double>& triangle_plane_distances,
                            const std::vector<double>& triangle_areas, 
                            const std::vector<size_t>& incident_triangles,
                            Vec3d& out);
    
    /// Assign a velocity vector to each mesh vertex
    /// 
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


