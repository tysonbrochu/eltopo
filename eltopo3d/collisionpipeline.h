// ---------------------------------------------------------
//
//  collisionpipeline.h
//  Tyson Brochu 2011
//  
//  Encapsulates all collision detection and resolution functions.
//
// ---------------------------------------------------------

#ifndef EL_TOPO_COLLISIONPIPELINE_H
#define EL_TOPO_COLLISIONPIPELINE_H

#include <deque>
#include <options.h>
#include <vec.h>


bool check_edge_triangle_intersection_by_index(size_t edge_a, 
                                               size_t edge_b, 
                                               size_t triangle_a, 
                                               size_t triangle_b, 
                                               size_t triangle_c, 
                                               const std::vector<Vec3d>& m_positions, 
                                               bool verbose );

inline bool check_triangle_triangle_intersection(Vec3st triangle_a, 
                                                 Vec3st triangle_b, 
                                                 const std::vector<Vec3d>& positions );


class BroadPhase;
class DynamicSurface;
struct ImpactZone;

// A potentially colliding pair of primitives.  Each pair is a triple of size_ts:
//  elements 0 and 1 are the indices of the primitives involved.
//  element 2 specifies if the potential collision is point-triangle or edge-edge
typedef std::deque<Vec3st> CollisionCandidateSet;

// --------------------------------------------------------
///
/// A collision between a triangle and a vertex or between two edges
///
// --------------------------------------------------------

struct Collision
{
    
    Collision() :
    m_is_edge_edge( false ),
    m_vertex_indices( Vec4st(static_cast<size_t>(~0)) ),
    m_normal( Vec3d(UNINITIALIZED_DOUBLE) ),
    m_barycentric_coordinates( Vec4d(UNINITIALIZED_DOUBLE) ),
    m_relative_displacement( UNINITIALIZED_DOUBLE )
    {}   
    
    Collision( bool in_is_edge_edge, const Vec4st& in_vertex_indices, const Vec3d& in_normal, const Vec4d& in_barycentric_coordinates, double in_relative_displacement ) :
    m_is_edge_edge( in_is_edge_edge ),
    m_vertex_indices( in_vertex_indices ),
    m_normal( in_normal ),
    m_barycentric_coordinates( in_barycentric_coordinates ),
    m_relative_displacement( in_relative_displacement )
    {
        if ( !m_is_edge_edge ) { assert( m_barycentric_coordinates[0] == 1.0 ); }
    }
    
    // One or more vertices is shared between this Collision and other
    inline bool overlap_vertices( const Collision& other ) const;
    
    // ALL vertices are shared between this Collision and other
    inline bool same_vertices( const Collision& other ) const;
    
    // Are the two elements both edges
    bool m_is_edge_edge;
    
    // Which vertices are involved in the collision
    Vec4st m_vertex_indices;
    
    // Collision normal
    Vec3d m_normal;
    
    // Barycentric coordinates of the point of intersection
    Vec4d m_barycentric_coordinates;
    
    // Magnitude of relative motion over the timestep
    double m_relative_displacement;
    
};


// --------------------------------------------------------
///
/// Edge-triangle intersection info.
///
// --------------------------------------------------------

struct Intersection
{
    
    Intersection(size_t edge_index, 
                 size_t triangle_index ) :
    m_edge_index( edge_index ),
    m_triangle_index( triangle_index )
    {}
    
    size_t m_edge_index;
    size_t m_triangle_index;
    
};

// --------------------------------------------------------
///
/// The results of processing a group of collision candidates.
///
// --------------------------------------------------------

struct ProcessCollisionStatus
{
    ProcessCollisionStatus() :
    collision_found( false ),
    overflow(false),
    all_candidates_processed( false )
    {}
    
    /// Whether one or more collisions was found.
    bool collision_found;
    
    /// Whether the number of collision candidates overflowed the candidate container.
    bool overflow;
    
    /// Whether all collision candidates were processed, or if the processing was terminated early.
    /// This is not necessarily equivalent to (!overflow): processing might stop early without overflow.
    bool all_candidates_processed;
    
};


// --------------------------------------------------------
///
/// Encapsulates all collision detection and resolution.
///
// --------------------------------------------------------

class CollisionPipeline
{
    
public:
    
    CollisionPipeline( DynamicSurface& surface,
                      BroadPhase& broadphase,
                      double in_friction_coefficient );
    

    // ---------------------------------------------------------
    // Proximity detection and resolution

    /// Repulsion forces
    ///
    void handle_proximities( double dt );
    
    // ---------------------------------------------------------
    // Collision detection and resolution

    /// Sequential impulses
    ///
    bool handle_collisions( double dt );
    
    /// Get all collisions at once
    ///   
    bool detect_collisions( std::vector<Collision>& collisions );
    
    /// Get collisions involving vertices in the impact zones
    /// 
    bool detect_new_collisions( const std::vector<ImpactZone> impact_zones, 
                               std::vector<Collision>& collisions );
    
    /// Get any collisions involving an edge and a triangle
    ///
    void detect_collisions( size_t edge_index, size_t triangle_index, std::vector<Collision>& collisions );
    
    bool check_if_collision_persists( const Collision& collision );    
    
    // ---------------------------------------------------------
    // Intersection detection 
    
    /// Get the set of triangles intersected by the given segment. 
    /// Also returns the barycentric coordinate along the segment of each intersection.
    void get_triangle_intersections(const Vec3d& segment_point_a, 
                                    const Vec3d& segment_point_b,
                                    std::vector<double>& hit_ss,
                                    std::vector<size_t>& hit_triangles ) const;
    
    /// Count the number of triangles intersected by the given segment.
    size_t get_number_of_triangle_intersections( const Vec3d& segment_point_a, 
                                                const Vec3d& segment_point_b ) const;
    
    /// Using exact interseciton testing, count the number of triangle intersected by the given segment.
    size_t get_number_of_triangle_intersections_exact( const Vec3d& segment_point_a, 
                                                      const Vec3d& segment_point_b ) const;
    
    /// Test the given triangle against all other triangles in the mesh for intersection.
    bool check_triangle_vs_all_triangles_for_intersection( size_t tri_index );
    bool check_triangle_vs_all_triangles_for_intersection( const Vec3st& tri );
    
    /// Get all self-intersections in the surface
    void get_intersections( bool degeneracy_counts_as_intersection, 
                           bool use_new_positions, 
                           std::vector<Intersection>& intersections );
    
    /// Look for self-intersections, but stop when the first one is found
    void get_first_intersection( bool degeneracy_counts_as_intersection, 
                                bool use_new_positions, 
                                Intersection& intersections );
    
    /// Fire an assert if the mesh contains a self-intersection. Uses m_positions as the vertex locations.
    void assert_mesh_is_intersection_free( bool degeneracy_counts_as_intersection );              
    
    /// Using m_newpositions as the vertex locations, fire an assert if the mesh contains a self-intersection.
    void assert_predicted_mesh_is_intersection_free( bool degeneracy_counts_as_intersection ); 

    
    double m_friction_coefficient;
    
private: 
    
    friend class DynamicSurface;
    friend class EdgeCollapser;
    
    void apply_edge_edge_impulse( const Collision& collision, double impulse_magnitude, double dt );
    
    void apply_triangle_point_impulse( const Collision& collision, double impulse_magnitude, double dt );
    
    void apply_impulse( const Vec4d& alphas, 
                       const Vec4st& vertex_indices, 
                       double impulse_magnitude, 
                       const Vec3d& normal,
                       double dt );
    
    void add_point_candidates( size_t v,
                              bool return_solid,
                              bool return_dynamic,
                              CollisionCandidateSet& collision_candidates );
    
    void add_edge_candidates( size_t e,
                             bool return_solid,
                             bool return_dynamic,
                             CollisionCandidateSet& collision_candidates );
    
    void add_triangle_candidates( size_t t,
                                 bool return_solid,
                                 bool return_dynamic,
                                 CollisionCandidateSet& collision_candidates );
    
    void add_point_update_candidates( size_t v, CollisionCandidateSet& collision_candidates );
    
    
    bool detect_segment_segment_collision( const Vec3st& candidate, Collision& collision );
    
    bool detect_point_triangle_collision( const Vec3st& candidate, Collision& collision );
    
    /// Test the candidates for proximity and apply impulses
    ///
    void process_proximity_candidates( double dt,
                                      CollisionCandidateSet& candidates );
    
    void dynamic_point_vs_solid_triangle_proximities(double dt);
    void dynamic_triangle_vs_all_point_proximities(double dt);
    void dynamic_edge_vs_all_edge_proximities(double dt);  
    
    
    /// Test the candidates and fix any collisions with impulses
    ///
    void process_collision_candidates( double dt,
                                      CollisionCandidateSet& candidates,
                                      bool add_to_new_candidates,
                                      CollisionCandidateSet& new_candidates,
                                      ProcessCollisionStatus& status );
    
    /// Test the candidates and return collision info
    ///
    void test_collision_candidates( CollisionCandidateSet& candidates,
                                   std::vector<Collision>& collisions,
                                   ProcessCollisionStatus& status );
    
    bool any_collision( CollisionCandidateSet& candidates, Collision& collision );
    
    void dynamic_point_vs_solid_triangle_collisions( double dt,
                                                    bool collect_candidates,
                                                    CollisionCandidateSet& update_collision_candidates,
                                                    ProcessCollisionStatus& status );
    
    void dynamic_triangle_vs_all_point_collisions( double dt,
                                                  bool collect_candidates,
                                                  CollisionCandidateSet& update_collision_candidates,
                                                  ProcessCollisionStatus& status );
    
    void dynamic_edge_vs_all_edge_collisions( double dt,
                                             bool collect_candidates,
                                             CollisionCandidateSet& update_collision_candidates,
                                             ProcessCollisionStatus& status );
    
    DynamicSurface& m_surface;
    BroadPhase& m_broad_phase;
    
    
    
};


// --------------------------------------------------------
///
/// Determine if another collision has any vertices in common with this collision.
///
// --------------------------------------------------------

inline bool Collision::overlap_vertices( const Collision& other ) const
{
    for ( unsigned short i = 0; i < 4; ++i )
    {
        if ( m_vertex_indices[i] == other.m_vertex_indices[0] || 
            m_vertex_indices[i] == other.m_vertex_indices[1] || 
            m_vertex_indices[i] == other.m_vertex_indices[2] || 
            m_vertex_indices[i] == other.m_vertex_indices[3] )
        {
            return true;
        }
    }
    
    return false;
}

// --------------------------------------------------------
///
/// Determine if another collision has all the same vertices as this collision.
///
// --------------------------------------------------------

inline bool Collision::same_vertices( const Collision& other ) const
{
    bool found[4];
    for ( unsigned short i = 0; i < 4; ++i )
    {
        if ( m_vertex_indices[i] == other.m_vertex_indices[0] || 
            m_vertex_indices[i] == other.m_vertex_indices[1] || 
            m_vertex_indices[i] == other.m_vertex_indices[2] || 
            m_vertex_indices[i] == other.m_vertex_indices[3] )
        {
            found[i] = true;
        }
        else
        {
            found[i] = false;
        }
    }
    
    return ( found[0] && found[1] && found[2] && found[3] );
}


// --------------------------------------------------------
///
/// Returns true if the an edge from one of the triangles intersects the other triangle.
/// NOTE: Using this routine will produce duplicate checks.  Better to use check_edge_triangle_intersection where possible.
///
// --------------------------------------------------------

inline bool check_triangle_triangle_intersection(Vec3st triangle_a, 
                                                 Vec3st triangle_b, 
                                                 const std::vector<Vec3d>& positions )
{
    if ( triangle_a[0] == triangle_a[1] || triangle_b[0] == triangle_b[1] )    
    { 
        return false; 
    }
    
    if ( check_edge_triangle_intersection_by_index( triangle_a[0], triangle_a[1], 
                                                   triangle_b[0], triangle_b[1], triangle_b[2], 
                                                   positions, false ) )
    {
        return true;
    }
    
    if ( check_edge_triangle_intersection_by_index( triangle_a[1], triangle_a[2], 
                                                   triangle_b[0], triangle_b[1], triangle_b[2], 
                                                   positions, false ) )
    {
        return true;
    }
    
    if ( check_edge_triangle_intersection_by_index( triangle_a[2], triangle_a[0], 
                                                   triangle_b[0], triangle_b[1], triangle_b[2], 
                                                   positions, false ) )
    {
        return true;
    }
    
    if ( check_edge_triangle_intersection_by_index( triangle_b[0], triangle_b[1], 
                                                   triangle_a[0], triangle_a[1], triangle_a[2], 
                                                   positions, false ) )
    {
        return true;
    }
    
    if ( check_edge_triangle_intersection_by_index( triangle_b[1], triangle_b[2], 
                                                   triangle_a[0], triangle_a[1], triangle_a[2], 
                                                   positions, false ) )
    {
        return true;
    }
    
    if ( check_edge_triangle_intersection_by_index( triangle_b[2], triangle_b[0], 
                                                   triangle_a[0], triangle_a[1], triangle_a[2], 
                                                   positions, false ) )
    {
		return true;
    }
    
    return false;
}

#endif
