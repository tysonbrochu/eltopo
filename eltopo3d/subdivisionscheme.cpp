// ---------------------------------------------------------
//
//  subdivisionscheme.cpp
//  Tyson Brochu 2008
//  
//  A collection of interpolation schemes for generating vertex locations.
//
// ---------------------------------------------------------

// ---------------------------------------------------------
// Includes
// ---------------------------------------------------------

#include <subdivisionscheme.h>

#include <mat.h>
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

// --------------------------------------------------------
///
/// Midpoint scheme: simply places the new vertex at the midpoint of the edge
///
// --------------------------------------------------------

void MidpointScheme::generate_new_midpoint( size_t edge_index, const SurfTrack& surface, Vec3d& new_point )
{
    const NonDestructiveTriMesh& mesh = surface.m_mesh;
    const std::vector<Vec3d>& positions = surface.get_positions();
    size_t p1_index = mesh.m_edges[edge_index][0];
	size_t p2_index = mesh.m_edges[edge_index][1];   
    
    new_point = 0.5 * ( positions[ p1_index ] + positions[ p2_index ] );
}


// --------------------------------------------------------
///
/// Butterfly scheme: uses a defined weighting of nearby vertices to determine the new vertex location
///
// --------------------------------------------------------

void ButterflyScheme::generate_new_midpoint( size_t edge_index, const SurfTrack& surface, Vec3d& new_point )
{
    const NonDestructiveTriMesh& mesh = surface.m_mesh;
    const std::vector<Vec3d>& positions = surface.get_positions();
    
    size_t p1_index = mesh.m_edges[edge_index][0];
	size_t p2_index = mesh.m_edges[edge_index][1];
	
    size_t tri0 = mesh.m_edge_to_triangle_map[edge_index][0];
    size_t tri1 = mesh.m_edge_to_triangle_map[edge_index][1];
    
	size_t p3_index = mesh.get_third_vertex( mesh.m_edges[edge_index][0], mesh.m_edges[edge_index][1], mesh.get_triangle(tri0) );
	size_t p4_index = mesh.get_third_vertex( mesh.m_edges[edge_index][0], mesh.m_edges[edge_index][1], mesh.get_triangle(tri1) );
	
	size_t adj_edges[4] = { mesh.get_edge_index( p1_index, p3_index ),
        mesh.get_edge_index( p2_index, p3_index ),
        mesh.get_edge_index( p1_index, p4_index ),
        mesh.get_edge_index( p2_index, p4_index ) };
    
	size_t q_indices[4];
	
	for ( size_t i = 0; i < 4; ++i )
	{
		const std::vector<size_t>& adj_tris = mesh.m_edge_to_triangle_map[ adj_edges[i] ];
		if ( adj_tris.size() != 2 )
		{
            // abort
			new_point = 0.5 * ( positions[ p1_index ] + positions[ p2_index ] );
            return;
		}
		
		if ( adj_tris[0] == tri0 || adj_tris[0] == tri1 )
		{
			q_indices[i] = mesh.get_third_vertex( mesh.m_edges[ adj_edges[i] ][0], mesh.m_edges[ adj_edges[i] ][1], mesh.get_triangle( adj_tris[1] ) );
		}
		else
		{
			q_indices[i] = mesh.get_third_vertex( mesh.m_edges[ adj_edges[i] ][0], mesh.m_edges[ adj_edges[i] ][1], mesh.get_triangle( adj_tris[0] ) );
		}
	}
    
	new_point =   8. * positions[ p1_index ] + 8. * positions[ p2_index ] + 2. * positions[ p3_index ] + 2. * positions[ p4_index ]
    - positions[ q_indices[0] ] - positions[ q_indices[1] ] - positions[ q_indices[2] ] - positions[ q_indices[3] ];
    
	new_point *= 0.0625;
    
}


// --------------------------------------------------------
///
/// Get the set of vertices adjacent to the given vertices in order of triangle adjacency, starting with the given first edge.
///
// --------------------------------------------------------

bool ModifiedButterflyScheme::get_adjacent_vertices_ordered( size_t vertex_index, 
                                                            size_t first_edge,
                                                            const NonDestructiveTriMesh& mesh,
                                                            std::vector<size_t>& adjacent_vertices )
{
    adjacent_vertices.clear();
    
    size_t current_edge = first_edge;
    size_t current_triangle = mesh.m_edge_to_triangle_map[current_edge][0];
    
    if ( mesh.m_edge_to_triangle_map[current_edge].size() != 2 )
    {
        // edge not incident on two triangles
        return false;
    }
    
    size_t count = 0;
    while ( count++ < mesh.m_vertex_to_edge_map[vertex_index].size() + 1 )
    {
        // find the other edge contained in the current triangle which has the current vertex as one of its endpoints
        unsigned int i = 0;
        for ( ; i < 3; ++i )
        {
            size_t edge_index = mesh.m_triangle_to_edge_map[current_triangle][i];
            if ( edge_index == current_edge ) { continue; }
            
            if ( mesh.m_edge_to_triangle_map[current_edge].size() != 2 )
            {
                // edge not incident on two triangles
                return false;
            }
            
            const Vec2st& edge_i = mesh.m_edges[edge_index];
            if ( (edge_i[0] == vertex_index || edge_i[1] == vertex_index) && edge_index != current_edge )
            {
                
                // save the vertex before switching to the next edge
                if ( edge_i[0] != vertex_index ) 
                {
                    assert( edge_i[1] == vertex_index );
                    adjacent_vertices.push_back( edge_i[0] ); 
                }
                else
                {
                    assert( edge_i[1] != vertex_index );
                    adjacent_vertices.push_back( edge_i[1] ); 
                }
                
                // now update the current edge
                current_edge = edge_index;
                
                break;
            }
        }
        assert( i < 3 );
        
        if ( current_edge == first_edge )
        {
            // we're done
            return true;
        }
        
        // update the current triangle
        
        if ( mesh.m_edge_to_triangle_map[current_edge].size() != 2 )
        {
            return false;
        }
        
        if ( mesh.m_edge_to_triangle_map[current_edge][0] == current_triangle )
        {
            current_triangle = mesh.m_edge_to_triangle_map[current_edge][1];
        }
        else
        {
            assert( mesh.m_edge_to_triangle_map[current_edge][1] == current_triangle );
            current_triangle = mesh.m_edge_to_triangle_map[current_edge][0];
        }
    }
    
    assert( !"get_adjacent_vertices_ordered: did not complete loop" );
    return false;
    
}



// --------------------------------------------------------
///
/// Get weighting coefficients for modified Butterfly scheme.
///
// --------------------------------------------------------

void ModifiedButterflyScheme::get_coefficients( size_t vertex_degree, std::vector<double>& coefficients )
{
    
    coefficients.resize( vertex_degree );
    
    if ( vertex_degree >= 5 )
    {
        for ( unsigned int i = 0; i < vertex_degree; ++i )
        {
            double K = static_cast<double>(vertex_degree);
            coefficients[i] = ( 0.25 + cos( 2.0 * M_PI * i / K ) + 0.5 * cos( 4.0 * M_PI * i / K ) ) / K;
        }
    }
    else if ( vertex_degree == 3 )
    {   
        coefficients[0] = 5.0/12.0;
        coefficients[1] = -1.0/12.0;
        coefficients[2] = -1.0/12.0;            
    }
    else if ( vertex_degree == 4 )
    {
        coefficients[0] = 3.0/8.0;
        coefficients[1] = 0.0;
        coefficients[2] = -1.0/8.0;
        coefficients[3] = 0.0;
    }
    else
    {
        assert(false);
    }
}


// --------------------------------------------------------
///
/// Modified Butterfly scheme: Handles vertices with degree != 6 better than Butterfly.  [Zorin et al. 1996]
///
// --------------------------------------------------------

void ModifiedButterflyScheme::generate_new_midpoint( size_t edge_index, const SurfTrack& surface, Vec3d& new_point )
{
    const NonDestructiveTriMesh& mesh = surface.m_mesh;
    const std::vector<Vec3d>& positions = surface.get_positions();
    size_t p1_index = mesh.m_edges[edge_index][0];
	size_t p2_index = mesh.m_edges[edge_index][1];
    
    const std::vector<size_t>& p1_edges = mesh.m_vertex_to_edge_map[p1_index];
    size_t p1_degree = p1_edges.size();
    std::vector<size_t> p1_adjacent_vertices;
    bool locally_manifold = get_adjacent_vertices_ordered( p1_index, edge_index, mesh, p1_adjacent_vertices );      
    if ( !locally_manifold )
    {
        new_point = 0.5 * ( positions[p1_index] + positions[p2_index] );
        return;
    }
    assert( p1_adjacent_vertices.size() == p1_degree );
    
    const std::vector<size_t>& p2_edges = mesh.m_vertex_to_edge_map[p2_index];
    size_t p2_degree = p2_edges.size();
    std::vector<size_t> p2_adjacent_vertices;
    locally_manifold = get_adjacent_vertices_ordered( p2_index, edge_index, mesh, p2_adjacent_vertices );      
    if ( !locally_manifold )
    {
        new_point = 0.5 * ( positions[p1_index] + positions[p2_index] );
        return;
    }   
    assert( p2_adjacent_vertices.size() == p2_degree );
    
    if ( p1_degree == 6 && p2_degree == 6 )
    {
        // use the usual 10-point stencil
        new_point = 0.5 * ( positions[p1_index] + positions[p2_index] );
        new_point += 0.125 * ( positions[ p1_adjacent_vertices[1] ] + positions[ p1_adjacent_vertices[5] ] );
        new_point += -0.0625 * ( positions[ p1_adjacent_vertices[2] ] + positions[ p1_adjacent_vertices[4] ] );
        new_point += -0.0625 * ( positions[ p2_adjacent_vertices[2] ] + positions[ p2_adjacent_vertices[4] ] );
        return;
    }      
    
    // One or both vertices are irregular
    
    std::vector<double> p1_neighbour_coefficients;
    get_coefficients( p1_degree, p1_neighbour_coefficients );      
    assert( p1_neighbour_coefficients.size() == p1_adjacent_vertices.size() );
    
    std::vector<double> p2_neighbour_coefficients;
    get_coefficients( p2_degree, p2_neighbour_coefficients );               
    assert( p2_neighbour_coefficients.size() == p2_adjacent_vertices.size() );
    
    if ( p1_degree != 6 && p2_degree == 6 )
    {      
        new_point = 0.75 * positions[p1_index];
        for ( size_t i = 0; i < p1_degree; ++i )
        {
            new_point += p1_neighbour_coefficients[i] * positions[ p1_adjacent_vertices[i] ];
        }
        return;
        
    }
    else if ( p1_degree == 6 && p2_degree != 6 )
    {
        new_point = 0.75 * positions[p2_index];
        for ( size_t i = 0; i < p2_degree; ++i )
        {
            new_point += p2_neighbour_coefficients[i] * positions[ p2_adjacent_vertices[i] ];
        }
        return;
        
    }
    else
    {
        
        Vec3d p1_contribution = 0.75 * positions[p1_index];
        for ( size_t i = 0; i < p1_degree; ++i )
        {
            p1_contribution += p1_neighbour_coefficients[i] * positions[ p1_adjacent_vertices[i] ];
        }      
        
        Vec3d p2_contribution = 0.75 * positions[p2_index];
        for ( size_t i = 0; i < p2_degree; ++i )
        {
            p2_contribution += p2_neighbour_coefficients[i] * positions[ p2_adjacent_vertices[i] ];
        }
        
        new_point = 0.5 * ( p1_contribution + p2_contribution );
        return;
    }
    
    assert( !"Should not get here" );
    
}


// --------------------------------------------------------
///
/// Quadric error minimization scheme: places the new vertex at the location that minimizes the change in the quadric metric tensor along the edge.
///
// --------------------------------------------------------

void QuadraticErrorMinScheme::generate_new_midpoint( size_t edge_index, const SurfTrack& surface, Vec3d& new_point )
{
    const NonDestructiveTriMesh& mesh = surface.m_mesh;
    const std::vector<Vec3d>& positions = surface.get_positions();
    
    size_t v0 = mesh.m_edges[edge_index][0];
    size_t v1 = mesh.m_edges[edge_index][1];
    
    Mat33d Q;
    zero(Q);
    Vec3d b;
    zero(b);
    
    std::vector<size_t> triangles_counted;
    
    Mat<1,1,double> constant_dist;
    constant_dist.a[0] = 0;
    
    for ( size_t i = 0; i < mesh.m_vertex_to_triangle_map[v0].size(); ++i )
    {
        size_t t = mesh.m_vertex_to_triangle_map[v0][i];
        const Vec3d& plane_normal = surface.get_triangle_normal( t );
        Q += outer( plane_normal, plane_normal );
        b += dot( positions[v0], plane_normal ) * plane_normal;
        constant_dist.a[0] += dot( plane_normal, positions[v0] ) * dot( plane_normal, positions[v0] );
        triangles_counted.push_back(t);
    }
    
    for ( size_t i = 0; i < mesh.m_vertex_to_triangle_map[v1].size(); ++i )
    {
        size_t t = mesh.m_vertex_to_triangle_map[v1][i];
        
        bool already_counted = false;
        for ( size_t j = 0; j < triangles_counted.size(); ++j ) 
        {
            if ( t == triangles_counted[j] )
            {
                already_counted = true;
            }
        }
        
        if ( !already_counted )
        {
            const Vec3d& plane_normal = surface.get_triangle_normal( t );
            Q += outer( plane_normal, plane_normal );
            b += dot( positions[v1], plane_normal ) * plane_normal;
            constant_dist.a[0] += dot( plane_normal, positions[v1] ) * dot( plane_normal, positions[v1] );
        }
    }
    
    // Compute normal direction
    Vec3d normal = 0.5 * (surface.get_vertex_normal(v0) + surface.get_vertex_normal(v1));
    normalize(normal);
    
    Mat<3,1,double> n;
    n(0,0) = normal[0];
    n(1,0) = normal[1];
    n(2,0) = normal[2];
    
    // Compute edge midpoint
    Vec3d midpoint = 0.5 * (positions[v0] + positions[v1]);   
    Mat<3,1,double> m;
    m(0,0) = midpoint[0];
    m(1,0) = midpoint[1];
    m(2,0) = midpoint[2]; 
    
    Mat<3,1,double> d;
    d(0,0) = b[0];
    d(1,0) = b[1];
    d(2,0) = b[2];
    
    double LHS = 2.0 * (n.transpose()*Q*n).a[0];              // result of multiplication is Mat<1,1,double>, hence the .a[0]
    double RHS = ( 2.0 * (n.transpose()*d) - (n.transpose()*Q*m) - (m.transpose()*Q*n) ).a[0];
    
    double a;
    if ( fabs(LHS) > 1e-10 )
    {
        a = RHS / LHS;
    }
    else
    {
        a = 0.0;
    }
    
    Mat<3,1,double> v = m + (a * n);
    
    double v_error = (v.transpose() * Q * v - 2.0 * (v.transpose() * d) + constant_dist).a[0];
    double m_error = (m.transpose() * Q * m - 2.0 * (m.transpose() * d) + constant_dist).a[0];
    
    //assert( v_error < m_error + 1e-8 );
    
    if ( surface.m_verbose )
    {
        std::cout << "a: " << a << std::endl;
        std::cout << "error at v: " << v_error << std::endl;
        std::cout << "error at midpoint: " << m_error << std::endl;
    }
    
    new_point = Vec3d( v.a[0], v.a[1], v.a[2] );
    
}
