// ---------------------------------------------------------
//
//  meancurvature.h
//  Tyson Brochu 2008
//
//  Mesh driver for motion in the normal direction scaled by mean curvature
//
// ---------------------------------------------------------

#ifndef EL_TOPO_MEANCURVATURE_H
#define EL_TOPO_MEANCURVATURE_H

// ---------------------------------------------------------
// Nested includes
// ---------------------------------------------------------

#include <meshdriver.h>
#include <vec.h>

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

// ---------------------------------------------------------
//  Class definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Motion in the normal direction scaled by mean curvature.
///
// ---------------------------------------------------------

class MeanCurvatureDriver : public MeshDriver
{
public:
    
    /// Constructor specifying speed of flow and final signed distance field to compare against.
    /// 
    MeanCurvatureDriver( double in_curvature_multiplier, 
                        const Array3d& in_final_signed_distance, 
                        const Vec3d& in_final_domain_low, 
                        double in_final_domain_dx );
    
    /// Compute the area of a triangle associated with the specified vertex.  (See [Meyer et al. 2002].)
    ///
    static double mixed_area( size_t vertex_index, size_t triangle_index, const SurfTrack& surf );
    
    /// Compute mean curvature times normal at a vertex
    ///
    static void vertex_mean_curvature_normal( size_t vertex_index, const SurfTrack& surf, Vec3d& out );
    
    /// Compute MC * normal at a vertex and return the sum of weights used (for computing the time step restriction)
    ///
    static void vertex_mean_curvature_normal( size_t vertex_index, const SurfTrack& surf, Vec3d& out, double& weight_sum );
    
    /// Set velocities on each mesh vertex
    ///
    void set_predicted_vertex_positions( const SurfTrack& surf, std::vector<Vec3d>& predicted_positions, double current_t, double& adaptive_dt );
    
    /// Compute the distance from vertices to analytic surface, weighted by associated vertex area
    /// 
    double compute_l1_error( const SurfTrack& surf ) const;
    
    /// Compute the maximum distance from vertices to analytic surface
    /// 
    double compute_inf_error( const SurfTrack& surf ) const;
    
    /// Compute and output both L1 and L_inf errors
    ///
    void compute_error( const SurfTrack& surf, double current_t );
    
    /// Speed of motion
    ///
    double curvature_multiplier;
    
    //
    // For error computation: what the signed distance function "should" be at the end time
    //
    
    Array3d final_signed_distance;
    Vec3d final_domain_low;
    double final_domain_dx;
    
};

// ---------------------------------------------------------
//  Inline functions
// ---------------------------------------------------------


#endif
