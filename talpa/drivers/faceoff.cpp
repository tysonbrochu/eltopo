// ---------------------------------------------------------
//
//  faceoff.cpp
//  Tyson Brochu 2008
//
//  Mesh driver for motion in the normal direction using the faceoff method (entropy solution).
//
// ---------------------------------------------------------

// ---------------------------------------------------------
// Includes
// ---------------------------------------------------------

#include <faceoff.h>

#include <cassert>
#include <geometryinit.h>
#include <lapack_wrapper.h>
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
///
/// Get the eigen decomposition of a 3x3 matrix and return results in descending order of eigenvalue magnitude
///
// ---------------------------------------------------------

namespace {
    
    void eigen_decomposition_descending( const Mat33d& mat, double* eigenvalues, double* eigenvectors )
    {
        eigenvectors[0] = mat(0,0);
        eigenvectors[1] = mat(1,0);
        eigenvectors[2] = mat(2,0);
        eigenvectors[3] = mat(0,1);
        eigenvectors[4] = mat(1,1);
        eigenvectors[5] = mat(2,1); 
        eigenvectors[6] = mat(0,2);
        eigenvectors[7] = mat(1,2);
        eigenvectors[8] = mat(2,2);   
        
        double work[9];
        int info = ~0, n = 3, lwork = 9;
        LAPACK::get_eigen_decomposition( &n,  eigenvectors, &n, eigenvalues, work, &lwork, &info );
        
        assert( info == 0 );
        
        // Now put them in descending order
        swap( eigenvalues[0], eigenvalues[2] );
        
        swap( eigenvectors[0], eigenvectors[6] );
        swap( eigenvectors[1], eigenvectors[7] );
        swap( eigenvectors[2], eigenvectors[8] );
    }
    
}

// ---------------------------------------------------------
// Member function definitions
// ---------------------------------------------------------

void FaceOffDriver::initialize( const SurfTrack& )
{
}

// ---------------------------------------------------------
///
/// Get the quadric metric tensor at a vertex from the given incident triangles
///
// ---------------------------------------------------------

void FaceOffDriver::compute_quadric_metric_tensor( const std::vector<Vec3d>& triangle_normals, 
                                                  const std::vector<double>& triangle_areas, 
                                                  const std::vector<size_t>& incident_triangles,
                                                  Mat33d& quadric_metric_tensor ) 
{
    std::vector< Vec3d > N;
    std::vector< double > W;
    std::vector< double > d;
    
    for ( size_t i = 0; i < incident_triangles.size(); ++i )
    {
        size_t triangle_index = incident_triangles[i];
        N.push_back( triangle_normals[triangle_index] );
        W.push_back( triangle_areas[triangle_index] );
        
        assert( triangle_areas[triangle_index] > 0.0 );
    }
    
    zero( quadric_metric_tensor );
    
    // Ax = b from N^TWNx = N^TWd
    for ( size_t i = 0; i < N.size(); ++i )
    {
        quadric_metric_tensor(0,0) += N[i][0] * W[i] * N[i][0];
        quadric_metric_tensor(1,0) += N[i][1] * W[i] * N[i][0];
        quadric_metric_tensor(2,0) += N[i][2] * W[i] * N[i][0];
        
        quadric_metric_tensor(0,1) += N[i][0] * W[i] * N[i][1];
        quadric_metric_tensor(1,1) += N[i][1] * W[i] * N[i][1];
        quadric_metric_tensor(2,1) += N[i][2] * W[i] * N[i][1];
        
        quadric_metric_tensor(0,2) += N[i][0] * W[i] * N[i][2];
        quadric_metric_tensor(1,2) += N[i][1] * W[i] * N[i][2];
        quadric_metric_tensor(2,2) += N[i][2] * W[i] * N[i][2];
    }
    
    double val[3], vec[9];
    eigen_decomposition_descending( quadric_metric_tensor, val, vec );
    
    if ( ! ( val[0] >= -1e-10 && val[1] >= -1e-10 && val[2] >= -1e-10 ) )
    {
        std::cout << "quadric metric tensor is not positive semi-definite" << std::endl;
        assert(0);
    }
    
}

// ---------------------------------------------------------
///
/// Return intersection point between a set of planes in the least-squares sense
///
// ---------------------------------------------------------

void FaceOffDriver::intersection_point( const std::vector<Vec3d>& triangle_normals, 
                                       const std::vector<double>& triangle_plane_distances,
                                       const std::vector<double>& triangle_areas, 
                                       const std::vector<size_t>& incident_triangles,
                                       Vec3d& out )
{
    
    std::vector< Vec3d > N;
    std::vector< double > W;
    std::vector< double > d;
    
    for ( size_t i = 0; i < incident_triangles.size(); ++i )
    {
        size_t triangle_index = incident_triangles[i];
        N.push_back( triangle_normals[triangle_index] );
        W.push_back( triangle_areas[triangle_index] );
        d.push_back( triangle_plane_distances[triangle_index] );
    }
    
    Mat33d A(0,0,0,0,0,0,0,0,0);
    Vec3d b(0,0,0);
    
    compute_quadric_metric_tensor( triangle_normals, triangle_areas, incident_triangles, A );
    
    for ( size_t i = 0; i < N.size(); ++i )
    {
        
        b[0] += N[i][0] * W[i] * d[i];
        b[1] += N[i][1] * W[i] * d[i];
        b[2] += N[i][2] * W[i] * d[i];      
    }
    
    double eigenvalues[3];
    double work[9];
    int info = ~0, n = 3, lwork = 9;
    LAPACK::get_eigen_decomposition( &n, A.a, &n, eigenvalues, work, &lwork, &info );
    
    if ( info != 0 )
    {
        assert(0);
    }
    
    Vec3d displacement_vector(0,0,0);
    
    for ( unsigned int i = 0; i < 3; ++i )
    {
        if ( eigenvalues[i] > G_EIGENVALUE_RANK_RATIO * eigenvalues[2] )
        {
            Vec3d eigenvector( A(0,i), A(1,i), A(2,i) );
            displacement_vector += dot(eigenvector, b) * eigenvector / eigenvalues[i];
        }
    }
    
    
    
    out = displacement_vector;
}


// ---------------------------------------------------------
///
/// Assign a velocity vector to each vertex in the input mesh, and adjust the time step size
///
// ---------------------------------------------------------

void FaceOffDriver::set_predicted_vertex_positions( const SurfTrack& surf, 
                                                   std::vector<Vec3d>& new_positions, 
                                                   double current_t, 
                                                   double& adaptive_dt )
{
    const NonDestructiveTriMesh& mesh = surf.m_mesh;
    std::vector<double> triangle_areas;
    triangle_areas.reserve(mesh.num_triangles());
    std::vector<Vec3d> triangle_normals;
    triangle_normals.reserve(mesh.num_triangles());
    std::vector<Vec3d> triangle_centroids;
    triangle_centroids.reserve(mesh.num_triangles());
    std::vector<double> triangle_plane_distances;
    triangle_plane_distances.reserve(mesh.num_triangles());
    
    const std::vector<Vec3st>& tris = mesh.get_triangles();
    for ( size_t i = 0; i < tris.size(); ++i )
    {
        if ( tris[i][0] == tris[i][1] )
        {
            triangle_areas.push_back( 0 );
            triangle_normals.push_back( Vec3d(0,0,0) );
            triangle_centroids.push_back( Vec3d(0,0,0) );
        }
        else
        {
            triangle_areas.push_back( surf.get_triangle_area( i ) );
            triangle_normals.push_back( surf.get_triangle_normal( i ) );
            triangle_centroids.push_back( (surf.get_position(tris[i][0]) + surf.get_position(tris[i][1]) + surf.get_position(tris[i][2])) / 3 );
        }
        
        double switch_speed = (current_t >= 1.0) ? -speed : speed;
        triangle_plane_distances.push_back( adaptive_dt * switch_speed );
    }
    
    std::vector<Vec3d> displacements;
    displacements.resize( surf.get_num_vertices() );
    
    //
    // Null space smoothing
    //
    
    {
        for ( size_t i = 0; i < surf.get_num_vertices(); ++i )
        {
            const MeshSmoother& smoother = surf.m_smoother;
            smoother.null_space_smooth_vertex( i, triangle_areas, triangle_normals, triangle_centroids, displacements[i] );
        }
    }
    
    //
    // Primary space displacement
    //
    
    for ( size_t p = 0; p < surf.get_num_vertices(); ++p )
    {
        Vec3d normal_dispacement;
        intersection_point( triangle_normals, triangle_plane_distances, triangle_areas, mesh.m_vertex_to_triangle_map[p], normal_dispacement );
        displacements[p] += normal_dispacement;
        
        //
        // Entropy solution
        //
        
        if ( surf.m_mesh.m_vertex_to_triangle_map[p].empty() )
        {
            continue;
        }
        
        double sum_mu_l = 0, sum_mu = 0;
        
        const std::vector<size_t>& incident_triangles = mesh.m_vertex_to_triangle_map[p];
        
        for ( size_t j = 0; j < incident_triangles.size(); ++j )
        {
            size_t triangle_index = incident_triangles[j];
            
            const Vec3st& tri = surf.m_mesh.get_triangle( triangle_index );
            
            Vec3d edge_vector;
            if ( tri[0] == p )
            {
                edge_vector = surf.get_position(tri[1]) - surf.get_position(tri[2]);
            }
            else if ( tri[1] == p )
            {
                edge_vector = surf.get_position(tri[2]) - surf.get_position(tri[0]);
            }
            else
            {
                edge_vector = surf.get_position(tri[0]) - surf.get_position(tri[1]);
            }
            
            Vec3d s = cross( triangle_normals[triangle_index], edge_vector );   // orthogonal to normal and edge oposite vertex
            
            bool contracting = dot( s, displacements[p] ) >= 0.0;
            
            double cos_theta = dot( triangle_normals[triangle_index], normal_dispacement ) / mag(normal_dispacement);
            
            double mu = triangle_areas[triangle_index];
            if ( contracting )
            {
                mu *= cos_theta * cos_theta;
            }
            
            double li = fabs( triangle_plane_distances[triangle_index] ); 
            
            if ( contracting )
            {
                li /= fabs( cos_theta );
            }
            
            sum_mu_l += mu * li;
            sum_mu += mu;
        }
        
        double length = sum_mu_l / sum_mu;
        
        displacements[p] += length * normal_dispacement / mag(normal_dispacement);
        
    }
    
    double beta = MeshSmoother::compute_max_timestep_quadratic_solve( surf.m_mesh.get_triangles(), surf.get_positions(), displacements, false );
    
    adaptive_dt *= beta;
    
    for(size_t i = 0; i < surf.get_num_vertices(); i++)
    {
        new_positions[i] = surf.get_position(i) + beta * displacements[i];
    }
    
    
}


// ---------------------------------------------------------
///
/// Compute the L1 error between the current mesh state and the analytical final state
///
// ---------------------------------------------------------

double FaceOffDriver::compute_l1_error( const SurfTrack& surf )
{
    
    double total_error = 0.0;
    double total_area = 0.0;
    
    for ( size_t i = 0; i < surf.get_num_vertices(); ++i )
    {
        if ( surf.m_mesh.m_vertex_to_triangle_map[i].empty() )
        {
            // ignore deleted vertices
            continue;
        }
        
        double dist = fabs( signed_distance_entropy( surf.get_position(i), sphere_a_centre, sphere_b_centre, max_radius, interior_radius ) );
        
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
/// Compute the L_inf error between the current mesh state and the analytical final state
///
// ---------------------------------------------------------

double FaceOffDriver::compute_inf_error( const SurfTrack& surf )
{
    
    double max_error = -1.0;
    
    for ( size_t i = 0; i < surf.get_num_vertices(); ++i )
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

void FaceOffDriver::compute_error( const SurfTrack& surf, double current_t )
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

