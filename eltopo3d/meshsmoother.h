// ---------------------------------------------------------
//
//  meshsmoother.h
//  Tyson Brochu 2011
//  
//  Functions related to the tangent-space mesh smoothing operation.
//
// ---------------------------------------------------------


#ifndef EL_TOPO_MESHSMOOTHER_H
#define EL_TOPO_MESHSMOOTHER_H

// ---------------------------------------------------------
//  Nested includes
// ---------------------------------------------------------

#include <meshoperator.h>
#include <vector>

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

class SurfTrack;
template<unsigned int N, class T> struct Vec;
typedef Vec<3,double> Vec3d;
typedef Vec<2,size_t> Vec2st;
typedef Vec<3,size_t> Vec3st;

// ---------------------------------------------------------
//  Class definitions
// ---------------------------------------------------------

class MeshSmoother : public MeshOperator
{
    
public:
    
    MeshSmoother( SurfTrack& surf ) :
      MeshOperator( surf )
    {}
    
    /// NULL-space smoothing of all vertices
    ///
    void process_mesh();
    
    /// Compute the maximum timestep which does not invert a triangle normal
    ///
    static double compute_max_timestep_quadratic_solve( const std::vector<Vec3st>& tris, 
                                                       const std::vector<Vec3d>& positions, 
                                                       const std::vector<Vec3d>& displacements, 
                                                       bool verbose );   
    
    /// Find a new vertex location using NULL-space smoothing
    ///
    void null_space_smooth_vertex(size_t v, 
                                  const std::vector<double>& triangle_areas, 
                                  const std::vector<Vec3d>& triangle_normals, 
                                  const std::vector<Vec3d>& triangle_centroids, 
                                  Vec3d& displacement ) const;      
        
};

#endif

