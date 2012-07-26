// ---------------------------------------------------------
//
//  dynamicsurface.cpp
//  Tyson Brochu 2008
//  
//  A triangle mesh with associated vertex locations and  masses.  Query functions for getting geometry info.
//
// ---------------------------------------------------------

// ---------------------------------------------------------
// Includes
// ---------------------------------------------------------

#include <dynamicsurface.h>

#include <broadphasegrid.h>
#include <cassert>
#include <ccd_wrapper.h>
#include <collisionqueries.h>
#include <ctime>
#include <impactzonesolver.h>
#include <lapack_wrapper.h>
#include <mat.h>
#include <queue>
#include <runstats.h>
#include <vec.h>
#include <vector>
#include <wallclocktime.h>


// ---------------------------------------------------------
// Local constants, typedefs, macros
// ---------------------------------------------------------

// ---------------------------------------------------------
//  Extern globals
// ---------------------------------------------------------

// ---------------------------------------------------------
// Static function definitions
// ---------------------------------------------------------


// ---------------------------------------------------------
// Member function definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// DynamicSurface constructor.  Copy triangles and vertex locations.
///
// ---------------------------------------------------------

DynamicSurface::DynamicSurface( const std::vector<Vec3d>& vertex_positions, 
                               const std::vector<Vec3st>& triangles,
                               const std::vector<double>& masses,
                               double in_proximity_epsilon,
                               double in_friction_coefficient,
                               bool in_collision_safety,
                               bool in_verbose ) :
    m_proximity_epsilon( in_proximity_epsilon ),
    m_verbose( in_verbose ),   
    m_collision_safety( in_collision_safety ),
    m_masses( masses ), 
    m_mesh(), 
    m_broad_phase( new BroadPhaseGrid() ),
    m_collision_pipeline( *this, *m_broad_phase, in_friction_coefficient ),    // allocated and initialized in the constructor body
    m_aabb_padding( max( in_proximity_epsilon, 1e-4 ) ),
    pm_positions(vertex_positions), 
    pm_newpositions(vertex_positions),
    m_velocities(0)
{
    
    if ( m_verbose )
    {
        std::cout << "constructing dynamic surface" << std::endl;
    }
    
    // if masses not provided, set all to 1.0
    if ( m_masses.size() == 0 )
    {
        m_masses.resize( get_num_vertices(), 1.0 );
    }
    else 
    {
        for ( size_t i = 0; i < m_masses.size(); ++i )
        {
            if ( m_masses[i] == 0.0 )
            {
                std::cerr << "Zero-mass vertex detected.  Please use #INF for infinite-mass vertices." << std::endl;
                assert(0);
            }
            
        }
        
    }
    
    m_mesh.set_num_vertices( get_num_vertices() );   
    m_mesh.replace_all_triangles( triangles );
    
    
    if ( m_verbose )
    {
        std::cout << "constructed dynamic surface" << std::endl;
    }
    
}


// ---------------------------------------------------------
///
/// Destructor. Frees memory allocated by DynamicSurface for the broad phase and collision pipeline objects.
///
// ---------------------------------------------------------

DynamicSurface::~DynamicSurface()
{
    delete m_broad_phase;
}


// ---------------------------------------------------------
///
/// Compute the unsigned distance to the surface.
///
// ---------------------------------------------------------

double DynamicSurface::distance_to_surface( const Vec3d& p, size_t& closest_triangle ) const
{
    
    double padding = m_aabb_padding;
    double min_distance = BIG_DOUBLE;
    
    while ( min_distance == BIG_DOUBLE )
    {
        
        Vec3d xmin( p - Vec3d( padding ) );
        Vec3d xmax( p + Vec3d( padding ) );
        
        std::vector<size_t> nearby_triangles;   
        
        m_broad_phase->get_potential_triangle_collisions( xmin, xmax, true, true, nearby_triangles );
        
        for ( size_t j = 0; j < nearby_triangles.size(); ++j )
        {
            const Vec3st& tri = m_mesh.get_triangle( nearby_triangles[j] );
            
            if ( tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2] ) { continue; }
            
            double curr_distance;
            check_point_triangle_proximity( p, get_position(tri[0]), get_position(tri[1]), get_position(tri[2]), curr_distance );
            if ( curr_distance < padding )
            {   
                min_distance = min( min_distance, curr_distance );
                closest_triangle = nearby_triangles[j];
            }
        }
        
        padding *= 2.0;
        
    }
    
    return min_distance;
    
}


// --------------------------------------------------------
///
/// Break up the triangle mesh into connected components, determine surface IDs for all vertices.
///
// --------------------------------------------------------

void DynamicSurface::partition_surfaces( std::vector<size_t>& surface_ids, std::vector< std::vector< size_t> >& surfaces ) const
{
    
    static const size_t UNASSIGNED = (size_t) ~0;
    
    surfaces.clear();
    
    surface_ids.clear();
    surface_ids.resize( get_num_vertices(), UNASSIGNED );
    
    size_t curr_surface = 0;
    
    while ( true )
    { 
        size_t next_unassigned_vertex;
        for ( next_unassigned_vertex = 0; next_unassigned_vertex < surface_ids.size(); ++next_unassigned_vertex )
        {
            if ( m_mesh.m_vertex_to_edge_map[next_unassigned_vertex].empty() ) { continue; }
            
            if ( surface_ids[next_unassigned_vertex] == UNASSIGNED )
            {
                break;
            }
        }
        
        if ( next_unassigned_vertex == surface_ids.size() )
        {
            break;
        }
        
        std::queue<size_t> open;
        open.push( next_unassigned_vertex );
        
        std::vector<size_t> surface_vertices;
        
        while ( false == open.empty() )
        {
            size_t vertex_index = open.front();
            open.pop();
            
            if ( m_mesh.m_vertex_to_edge_map[vertex_index].empty() ) { continue; }
            
            if ( surface_ids[vertex_index] != UNASSIGNED )
            {
                assert( surface_ids[vertex_index] == curr_surface );
                continue;
            }
            
            surface_ids[vertex_index] = curr_surface;
            surface_vertices.push_back( vertex_index );
            
            const std::vector<size_t>& incident_edges = m_mesh.m_vertex_to_edge_map[vertex_index];
            
            for( size_t i = 0; i < incident_edges.size(); ++i )
            {
                size_t adjacent_vertex = m_mesh.m_edges[ incident_edges[i] ][0];
                if ( adjacent_vertex == vertex_index ) { adjacent_vertex = m_mesh.m_edges[ incident_edges[i] ][1]; }
                
                if ( surface_ids[adjacent_vertex] == UNASSIGNED )
                {
                    open.push( adjacent_vertex );
                }
                else
                {
                    assert( surface_ids[adjacent_vertex] == curr_surface );
                }
                
            } 
        }
        
        surfaces.push_back( surface_vertices );
        
        ++curr_surface;
        
    }
    
    //
    // assert all vertices are assigned and share volume IDs with their neighbours
    //
    
    for ( size_t i = 0; i < surface_ids.size(); ++i )
    {
        if ( m_mesh.m_vertex_to_edge_map[i].empty() ) { continue; }
        
        assert( surface_ids[i] != UNASSIGNED );
        
        const std::vector<size_t>& incident_edges = m_mesh.m_vertex_to_edge_map[i];    
        for( size_t j = 0; j < incident_edges.size(); ++j )
        {
            size_t adjacent_vertex = m_mesh.m_edges[ incident_edges[j] ][0];
            if ( adjacent_vertex == i ) { adjacent_vertex = m_mesh.m_edges[ incident_edges[j] ][1]; }
            assert( surface_ids[adjacent_vertex] == surface_ids[i] );         
        } 
        
    }
    
}

// --------------------------------------------------------
///
/// Compute all vertex normals (unweighted average).
///
// --------------------------------------------------------

void DynamicSurface::get_all_vertex_normals( std::vector<Vec3d>& normals ) const
{
    normals.resize( get_num_vertices() );
    for ( size_t i = 0; i < get_num_vertices(); ++i )
    {
        normals[i] = get_vertex_normal(i);
    }
}




// ---------------------------------------------------------
///
/// Compute rank of the quadric metric tensor at a vertex
///
// ---------------------------------------------------------

unsigned int DynamicSurface::classify_vertex( size_t v ) const
{     
    if ( m_mesh.m_vertex_to_triangle_map[v].empty() )     { return 0; }
    
    const std::vector<size_t>& incident_triangles = m_mesh.m_vertex_to_triangle_map[v];
    
    
    Mat33d A(0,0,0,0,0,0,0,0,0);
    
    for ( size_t i = 0; i < incident_triangles.size(); ++i )
    {
        size_t triangle_index = incident_triangles[i];
        Vec3d normal = get_triangle_normal(triangle_index);
        double w = get_triangle_area(triangle_index);
        
        A(0,0) += normal[0] * w * normal[0];
        A(1,0) += normal[1] * w * normal[0];
        A(2,0) += normal[2] * w * normal[0];
        
        A(0,1) += normal[0] * w * normal[1];
        A(1,1) += normal[1] * w * normal[1];
        A(2,1) += normal[2] * w * normal[1];
        
        A(0,2) += normal[0] * w * normal[2];
        A(1,2) += normal[1] * w * normal[2];
        A(2,2) += normal[2] * w * normal[2];
    }
    
    // get eigen decomposition
    double eigenvalues[3];
    double work[9];
    int info = ~0, n = 3, lwork = 9;
    LAPACK::get_eigen_decomposition( &n, A.a, &n, eigenvalues, work, &lwork, &info );
    
    if ( info != 0 )
    {
        if ( m_verbose )
        {
            std::cout << "Eigen decomp failed.  Incident triangles: " << std::endl;
            for ( size_t i = 0; i < incident_triangles.size(); ++i )
            {
                size_t triangle_index = incident_triangles[i];
                Vec3d normal = get_triangle_normal(triangle_index);
                double w = get_triangle_area(triangle_index);
                
                std::cout << "normal: ( " << normal << " )    ";  
                std::cout << "area: " << w << std::endl;
            }
        }
        return 4;
    }
    
    // compute rank of primary space
    unsigned int rank = 0;
    for ( unsigned int i = 0; i < 3; ++i )
    {
        if ( eigenvalues[i] > G_EIGENVALUE_RANK_RATIO * eigenvalues[2] )
        {
            ++rank;
        }
    }
    
    return rank;
    
}


// ---------------------------------------------------------
///
/// Advance mesh by one time step 
///
// ---------------------------------------------------------

void DynamicSurface::integrate( double desired_dt, double& actual_dt )
{     
    
    static const bool DEGEN_DOES_NOT_COUNT = false;   
    static const bool USE_NEW_POSITIONS = true;
    
    if ( m_verbose ) 
    {
        std::cout << "---------------------- El Topo: integration and collision handling --------------------" << std::endl;
    }
    
    double curr_dt = desired_dt;
    bool success = false;
        
    while ( !success )
    {
        
        const std::vector<Vec3d> saved_predicted_positions = get_newpositions();
        
        if ( m_collision_safety )
        {        
            std::vector<Intersection> intersections;
            m_collision_pipeline.get_intersections( DEGEN_DOES_NOT_COUNT, false, intersections );
            if ( !intersections.empty() )
            {
                std::cout << "Intersection in mesh before integration" << std::endl;
                assert(false);
            }
        }
        
        m_velocities.resize( get_num_vertices() );
        for(size_t i = 0; i < get_num_vertices(); i++)
        {
            m_velocities[i] = ( get_newposition(i) - get_position(i) ) / curr_dt;  
        }
        
        // Handle proximities
        
        if ( m_collision_safety )
        {
            m_collision_pipeline.handle_proximities( curr_dt );
        }
        
        if ( m_collision_safety )
        {        
            
            // Handle continuous collisions
            bool all_collisions_handled = false;
            
            all_collisions_handled = m_collision_pipeline.handle_collisions( curr_dt );
            
            // failsafe impact zones 
            
            ImpactZoneSolver impactZoneSolver( *this );
            
            bool solver_ok = all_collisions_handled;
            
            if ( !solver_ok )
            {
                //if ( m_verbose ) 
                { std::cout << "IIZ" << std::endl; }
                solver_ok = impactZoneSolver.inelastic_impact_zones( curr_dt );            
            }
            
            if ( !solver_ok )
            {
                //if ( m_verbose ) 
                { std::cout << "RIZ" << std::endl; }
                // punt to rigid impact zones
                solver_ok = impactZoneSolver.rigid_impact_zones( curr_dt );
            }  
            
            if ( !solver_ok )
            {
                // back up and try again:
                
                curr_dt = 0.5 * curr_dt;
                for ( size_t i = 0; i < get_num_vertices(); ++i )
                {
                    set_newposition(i, get_position(i) + 0.5 * (saved_predicted_positions[i] - get_position(i)) ) ;
                }
                
                continue;      
            }
            
            
            // verify intersection-free predicted mesh
            std::vector<Intersection> intersections;
            m_collision_pipeline.get_intersections( DEGEN_DOES_NOT_COUNT, USE_NEW_POSITIONS, intersections );
            
            if ( !intersections.empty() )
            {
                std::cout << "Intersection in predicted mesh." << std::endl;
                
                //            extern bool simplex_verbose;
                //            simplex_verbose = true;
                //
                //            std::vector<Intersection> curr_intersections;
                //            get_intersections( true, false, curr_intersections );
                //            
                //            simplex_verbose = false;
                //            
                //            std::cout << "intersections at time step beginning: " << curr_intersections.size() << std::endl;
                //            
                //            std::vector<Collision> collisions;
                //            bool detect_finished = m_collision_pipeline.detect_collisions( collisions );
                //            
                //            assert( detect_finished );
                //            std::cout << "\n ======== num collisions: " << collisions.size() << "\n" << std::endl;
                //            
                //            if ( collisions.size() == 0 )
                //            {
                //               for ( size_t i = 0; i < intersections.size(); ++i )
                //               {
                //                  size_t e = intersections[i].m_edge_index;
                //                  size_t t = intersections[i].m_triangle_index;
                //                  const Vec2st& edge = m_mesh.m_edges[e];
                //                  const Vec3st& triangle = m_mesh.get_triangle(t);
                //                                       
                //                  const Vec3d& e0 = m_newpositions[edge[0]];
                //                  const Vec3d& e1 = m_newpositions[edge[1]];
                //                  const Vec3d& t0 = m_newpositions[triangle[0]];
                //                  const Vec3d& t1 = m_newpositions[triangle[1]];
                //                  const Vec3d& t2 = m_newpositions[triangle[2]];
                //                     
                //                  simplex_verbose = true;
                //
                //                  segment_triangle_intersection( e0, edge[0], 
                //                                                 e1, edge[1],
                //                                                 t0, triangle[0], 
                //                                                 t1, triangle[1], 
                //                                                 t2, triangle[2], 
                //                                                 false, true );
                //                  
                //                  bool is_degenerate;
                //                  //bool non_exact_result = triangle_intersects_segment( t0, t1, t2, e0, e1, 1e-7, true, is_degenerate );
                //                  //std::cout << "approximate test result: " << non_exact_result << ", is_degenerate: " << is_degenerate << std::endl;
                //                  
                //                  std::vector<Collision> check_collisions;
                //                  m_collision_pipeline.detect_collisions( e, t, check_collisions );
                //                  std::cout << "Collisions for this edge-triangle pair: " << check_collisions.size() << std::endl;
                //                  
                //                  std::cout << "triangle area: " << get_triangle_area( t ) << ", edge length: " << get_edge_length( e ) << std::endl;
                //                  
                //                  simplex_verbose = false;
                //                  
                //
                //               }
                //               
                //               std::cout << "Intersection in predicted mesh but no collisions detected." << std::endl;
                //               
                //               write_binary_file( m_mesh, m_positions, m_masses, 0.0, "/Users/tyson/scratch/tbrochu/collisiondebug/current.bin" );
                //               write_binary_file( m_mesh, m_newpositions, m_masses, curr_dt, "/Users/tyson/scratch/tbrochu/collisiondebug/predicted.bin" );
                //               
                //               assert( false );
                //            
                //            }
                
                
                if ( all_collisions_handled )
                {
                    std::cout << "Intersection in predicted mesh but handle collisions returned ok." << std::endl;
                    assert( false );
                }
                
                //if ( m_verbose )
                {
                    std::cout << "Intersection in predicted mesh, cutting timestep." << std::endl;
                }
                
                // back up and try again:
                
                curr_dt = 0.5 * curr_dt;
                for ( size_t i = 0; i < get_num_vertices(); ++i )
                {
                    set_newposition( i, get_position(i) + 0.5 * ( saved_predicted_positions[i] - get_position(i) ) );
                }
                
                continue;      
                
            }                 
            
        }
        
        // Set m_positions
        set_positions_to_newpositions();
        
        if ( m_collision_safety )
        {
            m_collision_pipeline.assert_mesh_is_intersection_free( DEGEN_DOES_NOT_COUNT );
        }
        
        actual_dt = curr_dt;
        
        success = true;
        
    }
    
    static unsigned int step = 0;
    ++step;
    
}

// ---------------------------------------------------------
///
/// Construct static acceleration structure
///
// ---------------------------------------------------------

void DynamicSurface::rebuild_static_broad_phase()
{
    assert( m_collision_safety );
    m_broad_phase->update_broad_phase( *this, false );
}

// ---------------------------------------------------------
///
/// Construct continuous acceleration structure
///
// ---------------------------------------------------------

void DynamicSurface::rebuild_continuous_broad_phase()
{
    assert( m_collision_safety );
    m_broad_phase->update_broad_phase( *this, true );
}


// ---------------------------------------------------------
///
/// Update the broadphase elements incident to the given vertex
///
// ---------------------------------------------------------

void DynamicSurface::update_static_broad_phase( size_t vertex_index )
{
    const std::vector<size_t>& incident_tris = m_mesh.m_vertex_to_triangle_map[ vertex_index ];
    const std::vector<size_t>& incident_edges = m_mesh.m_vertex_to_edge_map[ vertex_index ];
    
    Vec3d low, high;
    vertex_static_bounds( vertex_index, low, high );
    m_broad_phase->update_vertex( vertex_index, low, high, vertex_is_solid(vertex_index) );
    
    for ( size_t t = 0; t < incident_tris.size(); ++t )
    {
        triangle_static_bounds( incident_tris[t], low, high );
        m_broad_phase->update_triangle( incident_tris[t], low, high, triangle_is_solid(incident_tris[t]) );
    }
    
    for ( size_t e = 0; e < incident_edges.size(); ++e )
    {
        edge_static_bounds( incident_edges[e], low, high );
        m_broad_phase->update_edge( incident_edges[e], low, high, edge_is_solid(incident_edges[e]) );
    }
    
}


// ---------------------------------------------------------
///
/// Update the broadphase elements incident to the given vertex, using current and predicted vertex positions
///
// ---------------------------------------------------------

void DynamicSurface::update_continuous_broad_phase( size_t vertex_index )
{
    assert( m_collision_safety );
    
    const std::vector<size_t>& incident_tris = m_mesh.m_vertex_to_triangle_map[ vertex_index ];
    const std::vector<size_t>& incident_edges = m_mesh.m_vertex_to_edge_map[ vertex_index ];
    
    Vec3d low, high;
    vertex_continuous_bounds( vertex_index, low, high );
    m_broad_phase->update_vertex( vertex_index, low, high, vertex_is_solid(vertex_index) );
    
    for ( size_t t = 0; t < incident_tris.size(); ++t )
    {
        triangle_continuous_bounds( incident_tris[t], low, high );
        m_broad_phase->update_triangle( incident_tris[t], low, high, triangle_is_solid(incident_tris[t]) );
    }
    
    for ( size_t e = 0; e < incident_edges.size(); ++e )
    {
        edge_continuous_bounds( incident_edges[e], low, high );
        m_broad_phase->update_edge( incident_edges[e], low, high, edge_is_solid(incident_edges[e]) );
    }
}


// ---------------------------------------------------------
///
/// Compute the (padded) AABB of a vertex
///
// ---------------------------------------------------------

void DynamicSurface::vertex_static_bounds(size_t v, Vec3d &xmin, Vec3d &xmax) const
{
    if ( m_mesh.m_vertex_to_triangle_map[v].empty() )
    {
        xmin = Vec3d(m_aabb_padding);
        xmax = -Vec3d(m_aabb_padding);
    }
    else
    {
        xmin = get_position(v) - Vec3d(m_aabb_padding);
        xmax = get_position(v) + Vec3d(m_aabb_padding);
    }
}

// ---------------------------------------------------------
///
/// Compute the AABB of an edge
///
// ---------------------------------------------------------

void DynamicSurface::edge_static_bounds(size_t e, Vec3d &xmin, Vec3d &xmax) const
{
    const Vec2st& edge = m_mesh.m_edges[e];
    if ( edge[0] == edge[1] )
    {
        xmin = Vec3d(m_aabb_padding);
        xmax = -Vec3d(m_aabb_padding); 
    }
    else
    {            
        minmax( get_position(edge[0]), get_position(edge[1]), xmin, xmax);
        xmin -= Vec3d(m_aabb_padding);
        xmax += Vec3d(m_aabb_padding);
    }
}

// ---------------------------------------------------------
///
/// Compute the AABB of a triangle
///
// ---------------------------------------------------------

void DynamicSurface::triangle_static_bounds(size_t t, Vec3d &xmin, Vec3d &xmax) const
{
    const Vec3st& tri = m_mesh.get_triangle(t);  
    if ( tri[0] == tri[1] )
    {
        xmin = Vec3d(m_aabb_padding);
        xmax = -Vec3d(m_aabb_padding);
    }
    else
    {      
        minmax(get_position(tri[0]), get_position(tri[1]), get_position(tri[2]), xmin, xmax);
        xmin -= Vec3d(m_aabb_padding);
        xmax += Vec3d(m_aabb_padding);
    }
}

// ---------------------------------------------------------
///
/// Compute the AABB of a continuous vertex
///
// ---------------------------------------------------------

void DynamicSurface::vertex_continuous_bounds(size_t v, Vec3d &xmin, Vec3d &xmax) const
{
    if ( m_mesh.m_vertex_to_triangle_map[v].empty() )
    {
        xmin = Vec3d(m_aabb_padding);
        xmax = -Vec3d(m_aabb_padding);
    }
    else
    {
        minmax( get_position(v), get_newposition(v), xmin, xmax);
        xmin -= Vec3d(m_aabb_padding);
        xmax += Vec3d(m_aabb_padding);
    }
}

// ---------------------------------------------------------
///
/// Compute the AABB of a continuous edge
///
// ---------------------------------------------------------

void DynamicSurface::edge_continuous_bounds(size_t e, Vec3d &xmin, Vec3d &xmax) const
{
    const Vec2st& edge = m_mesh.m_edges[e];   
    if ( edge[0] == edge[1] )
    {
        xmin = Vec3d(m_aabb_padding);
        xmax = -Vec3d(m_aabb_padding);
    }
    else
    {      
        minmax(get_position(edge[0]), get_newposition(edge[0]), 
               get_position(edge[1]), get_newposition(edge[1]), 
               xmin, xmax);
        xmin -= Vec3d(m_aabb_padding);
        xmax += Vec3d(m_aabb_padding);
    }
}

// ---------------------------------------------------------
///
/// Compute the AABB of a continuous triangle
///
// ---------------------------------------------------------

void DynamicSurface::triangle_continuous_bounds(size_t t, Vec3d &xmin, Vec3d &xmax) const
{
    const Vec3st& tri = m_mesh.get_triangle(t);
    if ( tri[0] == tri[1] )
    {
        xmin = Vec3d(m_aabb_padding);
        xmax = -Vec3d(m_aabb_padding);
    }
    else
    {
        minmax(get_position(tri[0]), get_newposition(tri[0]), 
               get_position(tri[1]), get_newposition(tri[1]), 
               get_position(tri[2]), get_newposition(tri[2]), 
               xmin, xmax);
        
        xmin -= Vec3d(m_aabb_padding);
        xmax += Vec3d(m_aabb_padding);
    }
}


// --------------------------------------------------------

static bool aabbs_intersect( const Vec3d& a_xmin, const Vec3d& a_xmax, const Vec3d& b_xmin, const Vec3d& b_xmax )
{
    if ( (a_xmin[0] <= b_xmax[0] && a_xmin[1] <= b_xmax[1] && a_xmin[2] <= b_xmax[2]) &&
        (a_xmax[0] >= b_xmin[0] && a_xmax[1] >= b_xmin[1] && a_xmax[2] >= b_xmin[2]) )
    {
        return true;
    }
    
    return false;
}


// --------------------------------------------------------

void DynamicSurface::check_static_broad_phase_is_up_to_date() const
{
    
    // Verify by running against the n^2 broad phase
    
    //
    // vertex vs triangle
    //
    
    for ( size_t i = 0; i < get_num_vertices(); ++i )
    {
        if ( m_mesh.vertex_is_deleted(i) ) { continue; }
        
        // First, accelerated broad phase
        
        Vec3d aabb_low, aabb_high;
        vertex_static_bounds( i, aabb_low, aabb_high );
        
        std::vector<size_t> overlapping_triangles;
        m_broad_phase->get_potential_triangle_collisions(aabb_low, aabb_high, true, true, overlapping_triangles); 
        
        // filter deleted triangles
        for ( int k = 0; k < (int)overlapping_triangles.size(); ++k )
        {
            if ( m_mesh.triangle_is_deleted( overlapping_triangles[k] ) )
            {
                overlapping_triangles.erase( overlapping_triangles.begin() + k );
                --k;
            }
        }
        
        // Second, brute force check
        
        std::vector<size_t> brute_force_overlapping_triangles;
        
        for ( size_t j = 0; j < m_mesh.num_triangles(); ++j )
        {
            if ( m_mesh.triangle_is_deleted(j) ) { continue; }
            
            Vec3d tri_aabb_low, tri_aabb_high;
            triangle_static_bounds( j, tri_aabb_low, tri_aabb_high );
            
            if ( aabbs_intersect( aabb_low, aabb_high, tri_aabb_low, tri_aabb_high ) )
            {
                brute_force_overlapping_triangles.push_back( j );            
            }
        }
        
        assert( overlapping_triangles.size() == brute_force_overlapping_triangles.size() );
        
        std::sort(overlapping_triangles.begin(), overlapping_triangles.end());
        std::sort(brute_force_overlapping_triangles.begin(), brute_force_overlapping_triangles.end());
        
        for ( size_t k = 0; k < overlapping_triangles.size(); ++k )
        {
            assert( overlapping_triangles[k] == brute_force_overlapping_triangles[k] );
        }
        
    }
    
    //
    // edge vs edge
    //
    
    for ( size_t i = 0; i < m_mesh.m_edges.size(); ++i )
    {
        if ( m_mesh.edge_is_deleted(i) ) { continue; }
        
        // Accelerated broad phase
        
        Vec3d aabb_low, aabb_high;
        edge_static_bounds( i, aabb_low, aabb_high );
        
        std::vector<size_t> overlapping_edges;
        m_broad_phase->get_potential_edge_collisions( aabb_low, aabb_high, true, true, overlapping_edges );
        
        // filter deleted edges
        for ( int k = 0; k < (int)overlapping_edges.size(); ++k )
        {
            if ( m_mesh.edge_is_deleted( overlapping_edges[k] ) )
            {
                overlapping_edges.erase( overlapping_edges.begin() + k );
                --k;
            }
        }
        
        // Brute force
        std::vector<size_t> brute_force_overlapping_edges;
        for ( size_t j = 0; j < m_mesh.m_edges.size(); ++j )
        {
            if ( m_mesh.edge_is_deleted(j) ) { continue; }
            
            Vec3d edge_aabb_low, edge_aabb_high;
            edge_static_bounds( j, edge_aabb_low, edge_aabb_high );
            
            if ( aabbs_intersect( aabb_low, aabb_high, edge_aabb_low, edge_aabb_high ) )
            {
                brute_force_overlapping_edges.push_back( j );
            }
        }
        
        if ( overlapping_edges.size() != brute_force_overlapping_edges.size() )
        {
            
            std::cout << "edge " << i << ": " << m_mesh.m_edges[i] << std::endl;
            std::cout << "overlapping_edges.size(): " << overlapping_edges.size() << std::endl;
            for ( size_t k = 0; k < overlapping_edges.size(); ++k )
            {
                std::cout << k << ": " << overlapping_edges[k] << std::endl;
            }
            
            std::cout << "brute_force_overlapping_edges.size(): " << brute_force_overlapping_edges.size() << std::endl;
            for ( size_t k = 0; k < brute_force_overlapping_edges.size(); ++k )
            {
                std::cout << k << ": " << brute_force_overlapping_edges[k] << std::endl;
            }
            
        }
        
        assert( overlapping_edges.size() == brute_force_overlapping_edges.size() );
        
        std::sort( overlapping_edges.begin(), overlapping_edges.end() );
        std::sort( brute_force_overlapping_edges.begin(), brute_force_overlapping_edges.end() );
        
        for ( size_t k = 0; k < overlapping_edges.size(); ++k )
        {
            assert( overlapping_edges[k] == brute_force_overlapping_edges[k] );
        }
    }
    
    //
    // triangle vs vertex
    //
    
    for ( size_t i = 0; i < m_mesh.num_triangles(); ++i )
    {
        if ( m_mesh.triangle_is_deleted(i) ) { continue; }
        
        // Accelerated broad phase
        
        Vec3d aabb_low, aabb_high;
        triangle_static_bounds( i, aabb_low, aabb_high );
        
        std::vector<size_t> overlapping_vertices;
        m_broad_phase->get_potential_vertex_collisions( aabb_low, aabb_high, true, true, overlapping_vertices );
        
        // filter deleted vertices
        for ( int k = 0; k < (int)overlapping_vertices.size(); ++k )
        {
            if ( m_mesh.vertex_is_deleted( overlapping_vertices[k] ) )
            {
                overlapping_vertices.erase( overlapping_vertices.begin() + k );
                --k;
            }
        }
        
        // Brute force
        std::vector<size_t> brute_force_overlapping_vertices;
        for ( size_t j = 0; j < get_num_vertices(); ++j )
        {
            if ( m_mesh.vertex_is_deleted(j) ) { continue; }
            
            Vec3d vertex_aabb_low, vertex_aabb_high;
            vertex_static_bounds( j, vertex_aabb_low, vertex_aabb_high );
            
            if ( aabbs_intersect( aabb_low, aabb_high, vertex_aabb_low, vertex_aabb_high ) )
            {
                brute_force_overlapping_vertices.push_back( j );
            }
        }
        
        if ( overlapping_vertices.size() != brute_force_overlapping_vertices.size() )
        {
            std::cout << "triangle " << i << ": " << m_mesh.get_triangle(i) << std::endl;
            std::cout << "overlapping_vertices.size(): " << overlapping_vertices.size() << std::endl;
            for ( size_t k = 0; k < overlapping_vertices.size(); ++k )
            {
                std::cout << k << ": " << overlapping_vertices[k] << " --- ";
                std::cout << "is deleted: " << m_mesh.vertex_is_deleted( overlapping_vertices[k] ) << std::endl;
            }
            
            std::cout << "brute_force_overlapping_vertices.size(): " << brute_force_overlapping_vertices.size() << std::endl;
            for ( size_t k = 0; k < brute_force_overlapping_vertices.size(); ++k )
            {
                std::cout << k << ": " << brute_force_overlapping_vertices[k] << " --- ";
                std::cout << "is deleted: " << m_mesh.vertex_is_deleted( brute_force_overlapping_vertices[k] ) << std::endl;
            }
        }
        
        assert( overlapping_vertices.size() == brute_force_overlapping_vertices.size() );
        
        std::sort( overlapping_vertices.begin(), overlapping_vertices.end() );
        std::sort( brute_force_overlapping_vertices.begin(), brute_force_overlapping_vertices.end() );
        
        for ( size_t k = 0; k < overlapping_vertices.size(); ++k )
        {
            assert( overlapping_vertices[k] == brute_force_overlapping_vertices[k] );
        }
        
    }
    
    
}



// --------------------------------------------------------

void DynamicSurface::check_continuous_broad_phase_is_up_to_date() const
{
    
    // Verify by running against the n^2 broad phase
    
    //
    // vertex vs triangle
    //
    
    for ( size_t i = 0; i < get_num_vertices(); ++i )
    {
        if ( m_mesh.vertex_is_deleted(i) ) { continue; }
        
        // First, accelerated broad phase
        
        Vec3d aabb_low, aabb_high;
        vertex_continuous_bounds( i, aabb_low, aabb_high );
        
        std::vector<size_t> overlapping_triangles;
        m_broad_phase->get_potential_triangle_collisions(aabb_low, aabb_high, true, true, overlapping_triangles); 
        
        // filter deleted triangles
        for ( int k = 0; k < (int)overlapping_triangles.size(); ++k )
        {
            if ( m_mesh.triangle_is_deleted( overlapping_triangles[k] ) )
            {
                overlapping_triangles.erase( overlapping_triangles.begin() + k );
                --k;
            }
        }
        
        // Second, brute force check
        
        std::vector<size_t> brute_force_overlapping_triangles;
        
        for ( size_t j = 0; j < m_mesh.num_triangles(); ++j )
        {
            if ( m_mesh.triangle_is_deleted(j) ) { continue; }
            
            Vec3d tri_aabb_low, tri_aabb_high;
            triangle_continuous_bounds( j, tri_aabb_low, tri_aabb_high );
            
            if ( aabbs_intersect( aabb_low, aabb_high, tri_aabb_low, tri_aabb_high ) )
            {
                brute_force_overlapping_triangles.push_back( j );            
            }
        }
        
        assert( overlapping_triangles.size() == brute_force_overlapping_triangles.size() );
        
        std::sort(overlapping_triangles.begin(), overlapping_triangles.end());
        std::sort(brute_force_overlapping_triangles.begin(), brute_force_overlapping_triangles.end());
        
        for ( size_t k = 0; k < overlapping_triangles.size(); ++k )
        {
            assert( overlapping_triangles[k] == brute_force_overlapping_triangles[k] );
        }
        
    }
    
    //
    // edge vs edge
    //
    
    for ( size_t i = 0; i < m_mesh.m_edges.size(); ++i )
    {
        if ( m_mesh.edge_is_deleted(i) ) { continue; }
        
        // Accelerated broad phase
        
        Vec3d aabb_low, aabb_high;
        edge_continuous_bounds( i, aabb_low, aabb_high );
        
        std::vector<size_t> overlapping_edges;
        m_broad_phase->get_potential_edge_collisions( aabb_low, aabb_high, true, true, overlapping_edges );
        
        // filter deleted edges
        for ( int k = 0; k < (int)overlapping_edges.size(); ++k )
        {
            if ( m_mesh.edge_is_deleted( overlapping_edges[k] ) )
            {
                overlapping_edges.erase( overlapping_edges.begin() + k );
                --k;
            }
        }
        
        // Brute force
        std::vector<size_t> brute_force_overlapping_edges;
        for ( size_t j = 0; j < m_mesh.m_edges.size(); ++j )
        {
            if ( m_mesh.edge_is_deleted(j) ) { continue; }
            
            Vec3d edge_aabb_low, edge_aabb_high;
            edge_continuous_bounds( j, edge_aabb_low, edge_aabb_high );
            
            if ( aabbs_intersect( aabb_low, aabb_high, edge_aabb_low, edge_aabb_high ) )
            {
                brute_force_overlapping_edges.push_back( j );
            }
        }
        
        if ( overlapping_edges.size() != brute_force_overlapping_edges.size() )
        {
            
            std::cout << "edge " << i << ": " << m_mesh.m_edges[i] << std::endl;
            std::cout << "overlapping_edges.size(): " << overlapping_edges.size() << std::endl;
            for ( size_t k = 0; k < overlapping_edges.size(); ++k )
            {
                std::cout << k << ": " << overlapping_edges[k] << std::endl;
            }
            
            std::cout << "brute_force_overlapping_edges.size(): " << brute_force_overlapping_edges.size() << std::endl;
            for ( size_t k = 0; k < brute_force_overlapping_edges.size(); ++k )
            {
                std::cout << k << ": " << brute_force_overlapping_edges[k] << std::endl;
            }
            
        }
        
        assert( overlapping_edges.size() == brute_force_overlapping_edges.size() );
        
        std::sort( overlapping_edges.begin(), overlapping_edges.end() );
        std::sort( brute_force_overlapping_edges.begin(), brute_force_overlapping_edges.end() );
        
        for ( size_t k = 0; k < overlapping_edges.size(); ++k )
        {
            assert( overlapping_edges[k] == brute_force_overlapping_edges[k] );
        }
    }
    
    //
    // triangle vs vertex
    //
    
    for ( size_t i = 0; i < m_mesh.num_triangles(); ++i )
    {
        if ( m_mesh.triangle_is_deleted(i) ) { continue; }
        
        // Accelerated broad phase
        
        Vec3d aabb_low, aabb_high;
        triangle_continuous_bounds( i, aabb_low, aabb_high );
        
        std::vector<size_t> overlapping_vertices;
        m_broad_phase->get_potential_vertex_collisions( aabb_low, aabb_high, true, true, overlapping_vertices );
        
        // filter deleted vertices
        for ( int k = 0; k < (int)overlapping_vertices.size(); ++k )
        {
            if ( m_mesh.vertex_is_deleted( overlapping_vertices[k] ) )
            {
                overlapping_vertices.erase( overlapping_vertices.begin() + k );
                --k;
            }
        }
        
        // Brute force
        std::vector<size_t> brute_force_overlapping_vertices;
        for ( size_t j = 0; j < get_num_vertices(); ++j )
        {
            if ( m_mesh.vertex_is_deleted(j) ) { continue; }
            
            Vec3d vertex_aabb_low, vertex_aabb_high;
            vertex_continuous_bounds( j, vertex_aabb_low, vertex_aabb_high );
            
            if ( aabbs_intersect( aabb_low, aabb_high, vertex_aabb_low, vertex_aabb_high ) )
            {
                brute_force_overlapping_vertices.push_back( j );
            }
        }
        
        if ( overlapping_vertices.size() != brute_force_overlapping_vertices.size() )
        {
            std::cout << "triangle " << i << ": " << m_mesh.get_triangle(i) << std::endl;
            std::cout << "overlapping_vertices.size(): " << overlapping_vertices.size() << std::endl;
            for ( size_t k = 0; k < overlapping_vertices.size(); ++k )
            {
                std::cout << k << ": " << overlapping_vertices[k] << " --- ";
                std::cout << "is deleted: " << m_mesh.vertex_is_deleted( overlapping_vertices[k] ) << std::endl;
            }
            
            std::cout << "brute_force_overlapping_vertices.size(): " << brute_force_overlapping_vertices.size() << std::endl;
            for ( size_t k = 0; k < brute_force_overlapping_vertices.size(); ++k )
            {
                std::cout << k << ": " << brute_force_overlapping_vertices[k] << " --- ";
                std::cout << "is deleted: " << m_mesh.vertex_is_deleted( brute_force_overlapping_vertices[k] ) << std::endl;
                
                Vec3d lo, hi;
                bool is_solid = vertex_is_solid( brute_force_overlapping_vertices[k] );
                m_broad_phase->get_vertex_aabb( brute_force_overlapping_vertices[k], is_solid, lo, hi );
                std::cout << "AABB: " << lo << " - " << hi << std::endl;
                std::cout << "x: " << pm_positions[brute_force_overlapping_vertices[k]] << ", new_x: " << pm_newpositions[brute_force_overlapping_vertices[k]] << std::endl;
                
                bool query_overlaps_broadphase_aabb = aabbs_intersect( aabb_low, aabb_high, lo, hi );
                std::cout << "query_overlaps_broadphase_aabb: " << query_overlaps_broadphase_aabb << std::endl;
                
                
                BroadPhaseGrid* grid_bf = static_cast<BroadPhaseGrid*>(m_broad_phase);
                
                const std::vector<Vec3st>& cells = grid_bf->m_dynamic_vertex_grid.m_elementidxs[ brute_force_overlapping_vertices[k] ];
                std::cout << "cells: " << std::endl;
                for ( size_t m = 0; m < cells.size(); ++m )
                {
                    std::cout << cells[m] << std::endl;
                }
                
            }
        }
        
        assert( overlapping_vertices.size() == brute_force_overlapping_vertices.size() );
        
        std::sort( overlapping_vertices.begin(), overlapping_vertices.end() );
        std::sort( brute_force_overlapping_vertices.begin(), brute_force_overlapping_vertices.end() );
        
        for ( size_t k = 0; k < overlapping_vertices.size(); ++k )
        {
            assert( overlapping_vertices[k] == brute_force_overlapping_vertices[k] );
        }
        
    }
    
}







