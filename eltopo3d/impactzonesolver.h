// ---------------------------------------------------------
//
//  impactzonesolver.h
//  Tyson Brochu 2011
//  
//  Encapsulates two impact zone solvers: inelastic impact zones, and rigid impact zones.
//
// ---------------------------------------------------------


#ifndef EL_TOPO_IMPACTZONES_H
#define EL_TOPO_IMPACTZONES_H

// ---------------------------------------------------------
//  Nested includes
// ---------------------------------------------------------

#include <collisionpipeline.h>
#include <vector>

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

class DynamicSurface;

// ---------------------------------------------------------
//  Class definitions
// ---------------------------------------------------------

// --------------------------------------------------------
///
/// Used in the simultaneous handling of collisions: a set of connected elements which are in collision
///
// --------------------------------------------------------

struct ImpactZone
{
    ImpactZone() :
    m_collisions(),
    m_all_solved( false )
    {}
    
    // Get the set of all vertices in this impact zone
    void get_all_vertices( std::vector<size_t>& vertices ) const;
    
    // Whether this ImpactZones shares vertices with other
    bool share_vertices( const ImpactZone& other ) const;
    
    // Set of collisions with connected vertices
    std::vector<Collision> m_collisions;  
    
    // Whether all collisions in this zone have been solved (i.e. no longer colliding)
    bool m_all_solved;
    
};


// --------------------------------------------------------
///
/// Process collisions using inelastic impact zone solver and rigid impact zone solver.
///
// --------------------------------------------------------

class ImpactZoneSolver
{
    
public:
    
    ImpactZoneSolver( DynamicSurface& surface );
    
    bool inelastic_impact_zones(double dt);
    
    bool rigid_impact_zones(double dt);
    
protected:
    
    /// iteratively run collision detection and inelastic projection on an active set of collisions
    ///
    bool iterated_inelastic_projection( ImpactZone& iz, double dt );
    
    /// attempt to set normal velocity to zero for all collisions in the impact zone
    ///
    bool inelastic_projection( const ImpactZone& iz );
    
    /// Compute the best-fit single rigid motion for a set of vertices.
    ///
    bool calculate_rigid_motion(double dt, std::vector<size_t>& vs);
    
    DynamicSurface& m_surface;
    
    /// For rigid impact zones, treat solid vertices as having high but not infinite mass.  Use this value for the mass.
    ///
    const double m_rigid_zone_infinite_mass;     
    
};


// ---------------------------------------------------------
//  Inline functions
// ---------------------------------------------------------

// --------------------------------------------------------
///
/// Extract the set of all vertices in all collisions in an ImpactZone
///
// --------------------------------------------------------

inline void ImpactZone::get_all_vertices( std::vector<size_t>& vertices ) const
{
    vertices.clear();
    for ( size_t i = 0; i < m_collisions.size(); ++i )
    {
        add_unique( vertices, m_collisions[i].m_vertex_indices[0] );
        add_unique( vertices, m_collisions[i].m_vertex_indices[1] );
        add_unique( vertices, m_collisions[i].m_vertex_indices[2] );
        add_unique( vertices, m_collisions[i].m_vertex_indices[3] );
    }
}


// --------------------------------------------------------
///
/// Determine whether another ImpactZone shares any vertices with this ImpactZone
///
// --------------------------------------------------------

inline bool ImpactZone::share_vertices( const ImpactZone& other ) const
{
    for ( size_t i = 0; i < m_collisions.size(); ++i )
    {
        for ( size_t j = 0; j < other.m_collisions.size(); ++j )
        {
            if ( m_collisions[i].overlap_vertices( other.m_collisions[j] ) )
            {
                return true;
            }
        }
    }
    
    return false;
}


#endif


