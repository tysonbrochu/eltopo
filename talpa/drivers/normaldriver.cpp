// ---------------------------------------------------------
//
//  normaldriver.cpp
//  Tyson Brochu 2008
//
//  Mesh driver for motion in the normal direction using vertex 
//  normals (area-weighted average of incident triangle normals).
//
// ---------------------------------------------------------

// ---------------------------------------------------------
// Includes
// ---------------------------------------------------------

#include <normaldriver.h>

#include <faceoff.h>
#include <geometryinit.h>
#include <meshsmoother.h>
#include <nondestructivetrimesh.h>
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
/// Assign the normal at each vertex as the vertex's velocity
///
// ---------------------------------------------------------


void NormalDriver::set_predicted_vertex_positions( const SurfTrack& surf, std::vector<Vec3d>& predicted_positions, double current_t, double& adaptive_dt )
{
    std::vector<Vec3d> displacements( surf.get_num_vertices(), Vec3d(0,0,0) );
    std::vector<Vec3d> velocities( surf.get_num_vertices(), Vec3d(0,0,0) );
    
    for ( unsigned int i = 0; i < surf.get_num_vertices(); ++i )
    {
        if ( surf.m_mesh.m_vertex_to_triangle_map[i].empty() ) 
        { 
            displacements[i] = Vec3d(0,0,0);
            continue;
        }
        
        Vec3d normal(0,0,0);
        double sum_areas = 0.0;
        for ( unsigned int j = 0; j < surf.m_mesh.m_vertex_to_triangle_map[i].size(); ++j )
        {
            double area = surf.get_triangle_area( surf.m_mesh.m_vertex_to_triangle_map[i][j] );
            normal += surf.get_triangle_normal( surf.m_mesh.m_vertex_to_triangle_map[i][j] ) * area;
            sum_areas += area;
        }
        //normal /= sum_areas;
        normal /= mag(normal);
        
        double switch_speed = (current_t >= 1.0) ? -speed : speed;
        velocities[i] = switch_speed * normal;
        displacements[i] = adaptive_dt * velocities[i];
    }
    
    double capped_dt = MeshSmoother::compute_max_timestep_quadratic_solve( surf.m_mesh.get_triangles(), surf.get_positions(), displacements, false );
    
    adaptive_dt = min( adaptive_dt, capped_dt );
    
    for ( unsigned int i = 0; i < surf.get_num_vertices(); ++i )
    {
        predicted_positions[i] = surf.get_position(i) + adaptive_dt * velocities[i];
    }
    
}

// ---------------------------------------------------------
///
/// Compute the L1 error between the current mesh state and the analytical final state
///
// ---------------------------------------------------------

double NormalDriver::compute_l1_error( const SurfTrack& surf )
{
    
    double total_error = 0.0;
    double total_area = 0.0;
    
    for ( unsigned int i = 0; i < surf.get_num_vertices(); ++i )
    {
        if ( surf.m_mesh.m_vertex_to_triangle_map[i].empty() )
        {
            // ignore deleted vertices
            continue;
        }
        
        double dist = fabs( signed_distance_entropy( surf.get_position(i), sphere_a_centre, sphere_b_centre, max_radius, interior_radius ) );
        
        double area = 0;
        for ( unsigned int j = 0; j < surf.m_mesh.m_vertex_to_triangle_map[i].size(); ++j )
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
/// Compute the L_inf error between the current mesh state and the analytical final state
///
// ---------------------------------------------------------

double NormalDriver::compute_inf_error( const SurfTrack& surf )
{
    
    double max_error = -1.0;
    
    for ( unsigned int i = 0; i < surf.get_num_vertices(); ++i )
    {
        if ( surf.m_mesh.m_vertex_to_triangle_map[i].empty() )
        {
            // ignore deleted vertices
            continue;
        }
        
        double dist = signed_distance_entropy( surf.get_position(i), sphere_a_centre, sphere_b_centre, max_radius, interior_radius );
        max_error = max( max_error, fabs(dist) );
    }
    
    return max_error;
    
}


// ---------------------------------------------------------
///
/// Compute and output L1 and L_inf error measures
///
// ---------------------------------------------------------

void NormalDriver::compute_error( const SurfTrack& surf, double current_t )
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




