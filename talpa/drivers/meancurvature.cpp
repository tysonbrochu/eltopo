// ---------------------------------------------------------
//
//  meancurvature.cpp
//  Tyson Brochu 2008
//
//  Mesh driver for motion in the normal direction scaled by mean curvature
//
// ---------------------------------------------------------

// ---------------------------------------------------------
// Includes
// ---------------------------------------------------------

#include <meancurvature.h>

#include <array3_utils.h>
#include <fstream>
#include <iomesh.h>
#include <runstats.h>
#include <surftrack.h>


// ---------------------------------------------------------
// Global externs
// ---------------------------------------------------------

// ---------------------------------------------------------
// Local constants, typedefs, macros
// ---------------------------------------------------------

// ---------------------------------------------------------
// Static and nonmember function definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
// Member function definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Constructor
///
// ---------------------------------------------------------

MeanCurvatureDriver::MeanCurvatureDriver( double in_curvature_multiplier, 
                                         const Array3d& in_final_signed_distance, 
                                         const Vec3d& in_domain_low, 
                                         double in_domain_dx ) : 
   curvature_multiplier( in_curvature_multiplier ),
   final_signed_distance( in_final_signed_distance ),
   final_domain_low( in_domain_low ),
   final_domain_dx( in_domain_dx )
{}

// ---------------------------------------------------------
///
/// Compute the area of a triangle associated with the specified vertex.  (See [Meyer et al. 2002].)
///
// ---------------------------------------------------------

double MeanCurvatureDriver::mixed_area( size_t vertex_index, size_t triangle_index, const SurfTrack& surf )
{
    const Vec3st& tri = surf.m_mesh.get_triangle(triangle_index);
    
    Vec2st opposite_edge;
    if ( vertex_index == tri[0] )
    {
        opposite_edge = Vec2st( tri[1], tri[2] );
    }
    else if ( vertex_index == tri[1] )
    {
        opposite_edge = Vec2st( tri[2], tri[0] );
    }
    else
    {
        opposite_edge = Vec2st( tri[0], tri[1] );
    }
    
    const Vec3d& a = surf.get_position(vertex_index);
    const Vec3d& b = surf.get_position(opposite_edge[0]);
    const Vec3d& c = surf.get_position(opposite_edge[1]);
    
    bool obtuse_triangle = ( ( dot(b-a, c-a) < 0.0 ) || ( dot(a-b, c-b) < 0.0 ) || ( dot(a-c, b-c) < 0.0 ) );
    
    if ( obtuse_triangle )
    {
        //std::cout << "obtuse_triangle " << triangle_index << ": " << tri << std::endl;
        
        if ( dot(b-a, c-a) < 0.0 )
        {
            // obtuse at a
            return 0.5 * surf.get_triangle_area( triangle_index );
        }
        else
        {
            // obtuse somewhere else
            return 0.25 * surf.get_triangle_area( triangle_index );
        }
    }
    else
    {
        // not obtuse, use voronoi area
        
        double cross_c = mag( cross( a-c, b-c ) );      
        double cot_c = dot( a-c, b-c) / cross_c;      
        
        double cross_b = mag( cross( a-b, c-b ) );      
        double cot_b = dot( a-b, c-b) / cross_b;      
        
        return 1.0 / 8.0 * (mag2(b-a) * cot_c + mag2(c-a) * cot_b);
    }
    
}



// ---------------------------------------------------------
///
/// Compute mean curvature times normal at a vertex
///
// ---------------------------------------------------------

void MeanCurvatureDriver::vertex_mean_curvature_normal( size_t vertex_index, const SurfTrack& surf, Vec3d& out )
{
    double dummy;
    vertex_mean_curvature_normal( vertex_index, surf, out, dummy );
}

// ---------------------------------------------------------
///
/// Compute mean curvature times normal at a vertex and return the sum of weights used (for computing the time step restriction)
///
// ---------------------------------------------------------

void MeanCurvatureDriver::vertex_mean_curvature_normal( size_t vertex_index, const SurfTrack& surf, Vec3d& out, double& weight_sum )
{
    Vec3d mean_curvature_normal( 0, 0, 0 );
    weight_sum = 0;
    
    double edge_length_sum = 0.0;
    
    for ( size_t i = 0; i < surf.m_mesh.m_vertex_to_edge_map[vertex_index].size(); ++i )
    {
        size_t e = surf.m_mesh.m_vertex_to_edge_map[vertex_index][i];
        const Vec2st& curr_edge = surf.m_mesh.m_edges[e];
        Vec3d edge_vector;
        if ( curr_edge[0] == vertex_index )
        {
            edge_vector = surf.get_position( curr_edge[1] ) - surf.get_position( vertex_index );
        }
        else
        {
            assert( curr_edge[1] == vertex_index );
            edge_vector = surf.get_position( curr_edge[0] ) - surf.get_position( vertex_index );
        }
        
        edge_length_sum += mag( edge_vector );
        
        if ( surf.m_mesh.m_edge_to_triangle_map[e].size() != 2 )
        {
            // TODO: properly handle more than 2 incident triangles
            continue;
        }
        
        size_t tri0 = surf.m_mesh.m_edge_to_triangle_map[e][0];
        size_t tri1 = surf.m_mesh.m_edge_to_triangle_map[e][1];
        
        size_t third_vertex_0 = surf.m_mesh.get_third_vertex( curr_edge[0], curr_edge[1], surf.m_mesh.get_triangle(tri0) );
        size_t third_vertex_1 = surf.m_mesh.get_third_vertex( curr_edge[0], curr_edge[1], surf.m_mesh.get_triangle(tri1) );
        
        Vec3d v00 = surf.get_position( curr_edge[0] ) - surf.get_position( third_vertex_0 );
        Vec3d v10 = surf.get_position( curr_edge[1] ) - surf.get_position( third_vertex_0 );
        
        double cross_0 = mag( cross( v00, v10 ) );
        if ( cross_0 < 1e-10 )
        {
            continue;
        }
        double cot_0 = dot(v00, v10) / cross_0;
        
        Vec3d v01 = surf.get_position( curr_edge[0] ) - surf.get_position( third_vertex_1 );
        Vec3d v11 = surf.get_position( curr_edge[1] ) - surf.get_position( third_vertex_1 );
        
        double cross_1 = mag( cross( v01, v11 ) );
        if ( cross_1 < 1e-10 )
        {
            continue;
        }
        
        double cot_1 = dot(v01, v11) / cross_1;
        
        double weight = cot_0 + cot_1;
        weight_sum += weight;
        
        mean_curvature_normal += weight * edge_vector;
        
    }
    
    double vertex_area = 0.0;
    for ( size_t i = 0; i < surf.m_mesh.m_vertex_to_triangle_map[vertex_index].size(); ++i )
    {
        vertex_area += mixed_area( vertex_index, surf.m_mesh.m_vertex_to_triangle_map[vertex_index][i], surf );
    }
    
    double coeff = 1.0 / (2.0 * vertex_area);
    
    weight_sum *= coeff;
    
    out = coeff * mean_curvature_normal;
    
}


// ---------------------------------------------------------
///
/// Set velocities on each mesh vertex
///
// ---------------------------------------------------------

void MeanCurvatureDriver::set_predicted_vertex_positions( const SurfTrack& surf, std::vector<Vec3d>& predicted_positions, double /*current_t*/, double& adaptive_dt )
{
    
    size_t n = surf.get_num_vertices();
    predicted_positions.resize(n);
    
    double t_limit = BIG_DOUBLE;
    
    std::vector<Vec3d> velocities(n);
    
    for ( size_t i = 0; i < n; ++i )
    {
        if ( surf.m_mesh.m_vertex_to_triangle_map[i].size() > 0 )
        {
            double weight_sum;
            vertex_mean_curvature_normal( i, surf, velocities[i], weight_sum );
            
            velocities[i] *= curvature_multiplier;
            
            if ( weight_sum == 0 )
            {
                velocities[i] = Vec3d(0,0,0);
                continue;
            }
            
            t_limit = min( t_limit, 1.0 / (curvature_multiplier * weight_sum) );
            
        }
        else
        {
            velocities[i] = Vec3d(0,0,0);
        }
    }
    
    adaptive_dt = min( adaptive_dt, t_limit );   
    
    for ( size_t i = 0; i < n; ++i )
    {
        predicted_positions[i] = surf.get_position(i) + adaptive_dt * velocities[i];
    }
    
}


// ---------------------------------------------------------
///
/// Compute L1 error against the final-time signed distance field
///
// ---------------------------------------------------------

double MeanCurvatureDriver::compute_l1_error( const SurfTrack& surf ) const
{
    
    double total_error = 0.0;
    double total_area = 0.0;
    
    for ( size_t i = 0; i < surf.get_num_vertices(); ++i )
    {
        if ( surf.m_mesh.m_vertex_to_triangle_map[i].empty() )
        {
            continue;
        }
        
        const Vec3d& p = surf.get_position(i);
        double dist = cubic_interpolate_value( (1.0/final_domain_dx) * (p - final_domain_low), final_signed_distance );
        
        double area = 0;
        for ( size_t j = 0; j < surf.m_mesh.m_vertex_to_triangle_map[i].size(); ++j )
        {
            area += surf.get_triangle_area( surf.m_mesh.m_vertex_to_triangle_map[i][j] );
        }
        area /= 3;
        
        total_error += dist * area;
        total_area += area;
    }
    
    total_error /= total_area;
    
    return total_error;
    
}

// ---------------------------------------------------------
///
/// Compute max error against the final-time signed distance field
///
// ---------------------------------------------------------

double MeanCurvatureDriver::compute_inf_error( const SurfTrack& surf ) const
{   
    double max_error = -1.0;
    
    for ( size_t i = 0; i < surf.get_num_vertices(); ++i )
    {
        if ( surf.m_mesh.m_vertex_to_triangle_map[i].empty() )
        {
            continue;
        }
        
        const Vec3d& p = surf.get_position(i);
        double dist = cubic_interpolate_value( (1.0/final_domain_dx) * (p - final_domain_low), final_signed_distance );
        
        max_error = max ( max_error, dist );
    }
    
    return max_error;
    
}


// ---------------------------------------------------------
///
/// Compute and output current errors measures and dump complete error history to a file
///
// ---------------------------------------------------------

void MeanCurvatureDriver::compute_error( const SurfTrack& surf, double current_t )
{      
    double inf_error = compute_inf_error(surf);
    double l1_error = compute_l1_error(surf);   
    
    static unsigned int curr_frame = 0;
    extern RunStats g_stats;   
    g_stats.set_double( "last_t", current_t );   
    g_stats.set_double( "last_inf_error", inf_error );
    g_stats.set_double( "last_l1_error", l1_error );
    g_stats.update_min_double( "min_inf_error", fabs(inf_error) );
    g_stats.update_min_double( "min_l1_error", fabs(l1_error) );
    ++curr_frame;
    
}


