// ---------------------------------------------------------
//
//  geometryinit.cpp
//  Tyson Brochu 2008
//
//  A set of geometry initialization functions.
//
// ---------------------------------------------------------

// ---------------------------------------------------------
// Includes
// ---------------------------------------------------------

#include <geometryinit.h>
#include <array2.h>
#include <fstream>
#include <marching_tiles_hires.h>
#include <surftrack.h>

// ---------------------------------------------------------
// Global externs
// ---------------------------------------------------------

// ---------------------------------------------------------
// Local constants, typedefs, macros
// ---------------------------------------------------------

// ---------------------------------------------------------
// Static and non-member function definitions
// ---------------------------------------------------------

void append_mesh( std::vector<Vec3st>& tris, 
                 std::vector<Vec3d>& verts,
                 std::vector<double>& masses,
                 const std::vector<Vec3st>& new_tris, 
                 const std::vector<Vec3d>& new_verts,
                 const std::vector<double>& new_masses )
{
    size_t old_num_verts = verts.size();
    
    for ( size_t i = 0; i < new_verts.size(); ++i )
    {
        verts.push_back( new_verts[i] );
    }
    
    for ( size_t i = 0; i < new_tris.size(); ++i )
    {
        tris.push_back( new_tris[i] + Vec3st(old_num_verts) );
    }
    
    for ( size_t i = 0; i < new_masses.size(); ++i )
    {
        masses.push_back( new_masses[i] );
    }
    
}


// ---------------------------------------------------------
///
/// Use MarchingTilesHiRes to create a triangle mesh from a signed distance field
///
// ---------------------------------------------------------

void contour_phi( const Vec3d& domain_low, double domain_dx, Array3d& phi, std::vector<Vec3st>& tris, std::vector<Vec3d>& verts )
{
    MarchingTilesHiRes tiles( domain_low, domain_dx, phi );
    
    std::cout << "Contouring..." << std::endl;
    tiles.contour();
    
    std::cout << "Improving..." << std::endl;
    tiles.improve_mesh();
    
    std::cout << "Copying..." << std::endl;
    for ( size_t i = 0; i < tiles.tri.size(); ++i )
    {
        tris.push_back( Vec3st( tiles.tri[i] ) );
    }
    
    for ( size_t i = 0; i < tiles.x.size(); ++i )
    {
        verts.push_back( tiles.x[i] );
    }   
    
    std::cout << "done" << std::endl;   
}


void create_circle( std::vector<Vec3d>& verts, 
                   std::vector<Vec3st>& tris, 
                   std::vector<double>& masses,
                   const Vec3d& centre,
                   double radius,
                   size_t nx )
{
    Vec3d low = centre - Vec3d( radius, 0.0, radius );
    
    double dx = 2.0 * radius / (double)nx;
    
    Array2<Vec3d> grid_points;
    grid_points.resize(nx, nx);
    for ( size_t i = 0; i < nx; ++i )
    {
        for ( size_t j = 0; j < nx; ++j )
        {
            grid_points(i,j) = low + dx * Vec3d( i, 0, j );
            verts.push_back( low + dx * Vec3d( i, 0, j ) );
        }
    }
    
    // cells
    for ( size_t i = 0; i < nx-1; ++i )
    {
        for ( size_t j = 0; j < nx-1; ++j )
        {
            size_t a = i   + nx*j;
            size_t b = i+1 + nx*j;
            size_t c = i+1 + nx*(j+1);
            size_t d = i   + nx*(j+1);
            
            
            if ( ( dist( verts[a], centre ) < radius ) &&
                ( dist( verts[b], centre ) < radius ) &&
                ( dist( verts[d], centre ) < radius ) )
            {
                tris.push_back( Vec3st( a, b, d ) );
            }
            
            if ( ( dist( verts[b], centre ) < radius ) &&
                ( dist( verts[c], centre ) < radius ) &&
                ( dist( verts[d], centre ) < radius ) )            
            {
                tris.push_back( Vec3st( b, c, d ) );
            }
        }
    }
    
    masses.clear();
    masses.resize( verts.size(), 1.0 );
    
    NonDestructiveTriMesh temp_tri_mesh;
    temp_tri_mesh.set_num_vertices( verts.size() );
    for ( size_t i = 0; i < tris.size(); ++i )
    {
        temp_tri_mesh.add_triangle( tris[i] );
    }
    
}


std::vector<Vec2d> rest_vertices;


// ---------------------------------------------------------
///
/// Create a vertical or horizontal, triangulated sheet.
///
// ---------------------------------------------------------

void create_sheet( std::vector<Vec3d>& verts, 
                  std::vector<Vec3st>& tris, 
                  const Vec3d& low_corner, 
                  double dx, 
                  size_t nx, 
                  size_t ny )
{
    
    for(size_t i = 0; i < ny; i++)
    {
        for(size_t j = 0; j < nx; j++)
        {
            
            // plane normal is pointing in +y direction
            Vec3d offset;
            offset[0] = dx*j;  //width*((double)j/nx); 
            offset[1] = 0.0;
            offset[2] = dx*i;  //width*((double)i/ny);              
            
            verts.push_back( low_corner + offset );
            
            rest_vertices.push_back( Vec2d(offset[0], offset[2]) );         
        }
    }
    
    for(size_t i = 0; i < ny-1; i++)
    {
        for(size_t j = 0; j < nx-1; j++)
        {
            size_t idx = i*(nx)+j;
            tris.push_back(Vec3st(idx, idx+(nx), idx+1));
            tris.push_back(Vec3st(idx+1, idx+(nx), idx+(nx)+1));
        }
    }
}


// ---------------------------------------------------------

void create_curved_sheet( std::vector<Vec3d>& verts, 
                         std::vector<Vec3st>& tris, 
                         const Vec3d& low_corner, 
                         double dx, 
                         size_t nx, 
                         size_t ny )
{
    std::cout << "sheet: " << low_corner << ", " << dx << std::endl;
    std::cout << "resolution: " << nx << "x" << ny << std::endl;
    std::cout << "dimensions: " << dx*nx << "x" << dx*ny << std::endl;
    
    unsigned int num_curves = 4;
    double amplitude = dx*nx / 20.0;
    
    for(size_t i = 0; i < ny; i++)
    {
        for(size_t j = 0; j < nx; j++)
        {
            
            double theta = 2.0 * M_PI * (double)j / (double)nx * (double)num_curves;
            
            // plane normal is pointing in +y direction
            double x = dx*j;  //width*((double)j/nx); 
            double z = dx*i;  //width*((double)i/ny);              
            double y = 1e-3*z + amplitude * sin( theta );  
            verts.push_back(Vec3d(x,y,z) + low_corner);
            
        }
    }
    
    for(size_t i = 0; i < ny-1; i++)
    {
        for(size_t j = 0; j < nx-1; j++)
        {
            size_t idx = i*(nx)+j;
            tris.push_back(Vec3st(idx, idx+(nx), idx+1));
            tris.push_back(Vec3st(idx+1, idx+(nx), idx+(nx)+1));
        }
    }
}


//void create_curved_sheet( std::vector<Vec3d>& verts, 
//                          std::vector<Vec3st>& tris, 
//                          const Vec3d& low_corner, 
//                          const Vec3d& plane_normal, 
//                          double dx, 
//                          size_t nx, 
//                          size_t ny )
//{
//   std::cout << "sheet: " << low_corner << ", " << plane_normal << ", " << dx << std::endl;
//   std::cout << "resolution: " << nx << "x" << ny << std::endl;
//   std::cout << "dimensions: " << dx*nx << "x" << dx*ny << std::endl;
//   
//
//   double radius = nx/2 * dx / (2 * M_PI);
//
//   size_t num_periods = 2;
//   
//   double total_arc_length = num_periods * 2 * M_PI * radius;
//
//   size_t nk = (nx-1) / (2*num_periods);
//   
//   std::cout << "nk: " << nk << std::endl;
//   
//   double ds = total_arc_length / (double)(nx-1);
//   
//   
//   for( size_t i = 0; i < ny; i++ )
//   {
//      
//      double z = dx*i;
//      
//      for(size_t j = 0; j < nx; j++)
//      {
//
//         double s = j * ds;
//         
//         size_t k = j / nk;
//         
//         std::cout << "current half-period, k: " << k << ", ";
//         
//         size_t jk = j - k * nk;
//         
//         std::cout << "jk: " << jk << ", ";
//         std::cout << "s: " << s << ", ";
//         
//         Vec3d centre( radius + 2*k*radius, 0.0, 0.0 );
//         //Vec3d centre( radius, 0.0, 0.0 );
//         
//         std::cout << "centre: " << centre << std::endl;
//         
//         //double t = s - 2*jk * period_arc_length;
//         //std::cout << "t: " << t << std::endl;
//         
//         double t = s / radius;
//         
//         double x = -radius * cos( t ) + centre[0];
//         double y = radius * sin( t ) + centre[1];
//         
//         if ( k % 2 == 1 )
//         {
//            x = radius * cos( t ) + centre[0];
//         }
//        
////            double x = dx*j;
////            double x_c = x - radius;
////            double x_mod = x_c; // - 2.0 * radius * std::floor(x_c / (2.0*radius));
////            assert( radius*radius - x_mod*x_mod >= 0 );
////            double y = sqrt( radius*radius - x_mod*x_mod );
//         
//         verts.push_back(Vec3d(x,y,z) + low_corner);
//         
//         
//         Vec2d rest_vert( s, z );
//         rest_vertices.push_back( rest_vert );
//         
//         
//      }
//   }
//   
//   for(size_t i = 0; i < ny-1; i++)
//   {
//      for(size_t j = 0; j < nx-1; j++)
//      {
//         size_t idx = i*(nx)+j;
//         tris.push_back(Vec3st(idx, idx+(nx), idx+1));
//         tris.push_back(Vec3st(idx+1, idx+(nx), idx+(nx)+1));
//      }
//   }
//}
//


// ---------------------------------------------------------
///
/// Project mesh vertices onto an analytic cube
///
// ---------------------------------------------------------

void project_to_exact_cube( std::vector<Vec3d>& verts, 
                           const Vec3d& cube_low, 
                           const Vec3d& cube_high )
{
    
    for ( size_t i = 0; i < verts.size(); ++i )
    {
        Vec3d& v = verts[i];
        Vec3d p(0,0,0);
        double min_dist = 1e+30;
        
        if ( fabs( v[0] - cube_low[0] ) < min_dist ) { min_dist = fabs( v[0] - cube_low[0] ); p = v; p[0] = cube_low[0]; }
        if ( fabs( v[1] - cube_low[1] ) < min_dist ) { min_dist = fabs( v[1] - cube_low[1] ); p = v; p[1] = cube_low[1]; }
        if ( fabs( v[2] - cube_low[2] ) < min_dist ) { min_dist = fabs( v[2] - cube_low[2] ); p = v; p[2] = cube_low[2]; }
        if ( fabs( v[0] - cube_high[0] ) < min_dist ) { min_dist = fabs( v[0] - cube_high[0] ); p = v; p[0] = cube_high[0]; }
        if ( fabs( v[1] - cube_high[1] ) < min_dist ) { min_dist = fabs( v[1] - cube_high[1] ); p = v; p[1] = cube_high[1]; }
        if ( fabs( v[2] - cube_high[2] ) < min_dist ) { min_dist = fabs( v[2] - cube_high[2] ); p = v; p[2] = cube_high[2]; }
        
        v = p;
    }   
}

// ---------------------------------------------------------
///
/// Project mesh vertices onto an analytic sphere
///
// ---------------------------------------------------------

void project_to_exact_sphere( std::vector<Vec3d>& verts, const Vec3d& sphere_centre, double sphere_radius )
{
    for ( size_t i = 0; i < verts.size(); ++i )
    {
        Vec3d& v = verts[i];
        
        double dist = mag( v - sphere_centre ) - sphere_radius;
        
        v -= dist * ( v - sphere_centre ) / mag( v - sphere_centre );
    }
}


// ---------------------------------------------------------
///
/// Project mesh vertices onto an analytic dumbbell
///
// ---------------------------------------------------------

void project_to_exact_dumbbell( std::vector<Vec3d>& verts, 
                               const Vec3d& sphere_a_centre, 
                               const Vec3d sphere_b_centre, 
                               double sphere_radius, 
                               double handle_width )
{
    
    for ( size_t i = 0; i < verts.size(); ++i )
    {
        Vec3d& pt = verts[i];
        double dx = 1e-6;
        
        double x_pos = signed_distance_dumbbell( pt + Vec3d(dx,0,0), sphere_a_centre, sphere_b_centre, sphere_radius, handle_width );
        double x_neg = signed_distance_dumbbell( pt - Vec3d(dx,0,0), sphere_a_centre, sphere_b_centre, sphere_radius, handle_width );
        double y_pos = signed_distance_dumbbell( pt + Vec3d(0,dx,0), sphere_a_centre, sphere_b_centre, sphere_radius, handle_width );
        double y_neg = signed_distance_dumbbell( pt - Vec3d(0,dx,0), sphere_a_centre, sphere_b_centre, sphere_radius, handle_width );
        double z_pos = signed_distance_dumbbell( pt + Vec3d(0,0,dx), sphere_a_centre, sphere_b_centre, sphere_radius, handle_width );
        double z_neg = signed_distance_dumbbell( pt - Vec3d(0,0,dx), sphere_a_centre, sphere_b_centre, sphere_radius, handle_width );
        
        Vec3d n( x_pos - x_neg, y_pos - y_neg, z_pos - z_neg );
        normalize( n );
        
        double dist = signed_distance_dumbbell( pt, sphere_a_centre, sphere_b_centre, sphere_radius, handle_width );
        
        pt -= dist * n;
        
        assert( fabs( signed_distance_dumbbell( pt, sphere_a_centre, sphere_b_centre, sphere_radius, handle_width )) < 1e-8 ); 
    }
    
}


// ---------------------------------------------------------
///
/// Read a signed distance field from an ASCII file
///
// ---------------------------------------------------------

void read_signed_distance( const char* filename, Array3d& signed_distance )
{
    FILE* file = fopen( filename, "r" );
    if ( file == NULL )
    {
        std::cout << "failed to open file: " << filename << std::endl;
        assert(0);
        return;
    }
    
    int num_dim;
    fscanf( file, "%d ", &num_dim );
    
    if ( num_dim != 3 )
    {
        std::cout << "num_dim: " << num_dim << std::endl;
        assert(0);
    }
    
    int dims[3];
    fscanf( file, "%d %d %d ", &(dims[0]), &(dims[1]), &(dims[2]) );
    
    unsigned int num_elements = dims[0]*dims[1]*dims[2];
    
    signed_distance.ni = dims[0];
    signed_distance.nj = dims[1];
    signed_distance.nk = dims[2];
    
    signed_distance.a.resize(num_elements);
    
    //fread( &buf, sizeof(double), num_elements, file );
    for ( unsigned int i = 0; i < num_elements; ++i )
    {
        float buf;
        fscanf( file, "%f ", &buf );
        signed_distance.a[i] = buf;
    }
    
    fclose(file);
    
    std::cout << "Read mat file" << std::endl;
    
}



// ---------------------------------------------------------
///
/// Create signed distance field for a capsule with the given geometry
///
// ---------------------------------------------------------

void create_capsule_signed_distance( const Vec3d& capsule_end_a, 
                                    const Vec3d& capsule_end_b, 
                                    double capsule_radius,
                                    double dx,
                                    const Vec3d& domain_low,
                                    const Vec3d& domain_high,                                     
                                    Array3d& phi )
{ 
    phi.resize( (int) ceil( (domain_high[0]-domain_low[0]) / dx), (int) ceil( (domain_high[1]-domain_low[1]) / dx), (int) ceil( (domain_high[2]-domain_low[2]) / dx) );
    
    std::cout << "Generating signed distance function.  Grid resolution: " << phi.ni << " x " << phi.nj << " x " << phi.nk << std::endl;
    
    for ( int i = 0; i < phi.ni; ++i )
    {
        for ( int j = 0; j < phi.nj; ++j )
        {
            for ( int k = 0; k < phi.nk; ++k )
            {
                Vec3d pt = domain_low + dx * Vec3d(i,j,k);
                
                double distance;
                
                Vec3d central_segment(capsule_end_b - capsule_end_a);
                double m2=mag2(central_segment);
                
                // find parameter value of closest point on infinite line
                double s = dot(capsule_end_b - pt, central_segment)/m2;
                
                if ( s < 0.0 )
                {
                    // dist = distance to the cylinder disc at end b
                    
                    distance = dist(pt, capsule_end_b);
                    distance -= capsule_radius;
                    
                    
                } 
                else if ( s > 1.0 )
                {
                    // dist = distance to the cylinder disc at end b
                    distance = dist(pt, capsule_end_a );
                    distance -= capsule_radius;
                    
                }
                else
                {
                    // dist = distance to the cylinder's central axis
                    
                    distance = dist(pt, s*capsule_end_a + (1-s)*capsule_end_b);
                    distance -= capsule_radius;
                }
                
                phi(i,j,k) = distance;
                
            }
        }
    }
}


// ---------------------------------------------------------
///
/// Create signed distance field for a cube with the given geometry
///
// ---------------------------------------------------------

void create_cube_signed_distance( const Vec3d& cube_low, 
                                 const Vec3d& cube_high, 
                                 double dx,
                                 const Vec3d& domain_low,
                                 const Vec3d& domain_high,                                     
                                 Array3d& phi )
{
    phi.resize( (int) ceil( (domain_high[0]-domain_low[0]) / dx), (int) ceil( (domain_high[1]-domain_low[1]) / dx), (int) ceil( (domain_high[2]-domain_low[2]) / dx) );
    
    std::cout << "Generating signed distance function.  Grid resolution: " << phi.ni << " x " << phi.nj << " x " << phi.nk << std::endl;
    
    for ( int i = 0; i < phi.ni; ++i )
    {
        for ( int j = 0; j < phi.nj; ++j )
        {
            for ( int k = 0; k < phi.nk; ++k )
            {
                Vec3d pt = domain_low + dx * Vec3d(i,j,k);
                
                double dist_low_x = cube_low[0] - pt[0];
                double dist_high_x = pt[0] - cube_high[0];
                
                double dist_low_y = cube_low[1] - pt[1];
                double dist_high_y = pt[1] - cube_high[1];
                
                double dist_low_z = cube_low[2] - pt[2];
                double dist_high_z = pt[2] - cube_high[2];
                
                phi(i,j,k) = max( dist_low_x, dist_high_x, dist_low_y, dist_high_y, dist_low_z, dist_high_z );
                
            }
        }
    }
}


// ---------------------------------------------------------
///
/// Create signed distance field for a sphere with the given geometry
///
// ---------------------------------------------------------

void create_sphere_signed_distance( const Vec3d& sphere_centre, 
                                   double sphere_radius, 
                                   double dx,
                                   const Vec3d& domain_low,
                                   const Vec3d& domain_high,                                     
                                   Array3d& phi )
{
    phi.resize( (int) ceil( (domain_high[0]-domain_low[0]) / dx), (int) ceil( (domain_high[1]-domain_low[1]) / dx), (int) ceil( (domain_high[2]-domain_low[2]) / dx) );
    
    std::cout << "Generating signed distance function.  Grid resolution: " << phi.ni << " x " << phi.nj << " x " << phi.nk << std::endl;
    
    for ( int i = 0; i < phi.ni; ++i )
    {
        for ( int j = 0; j < phi.nj; ++j )
        {
            for ( int k = 0; k < phi.nk; ++k )
            {
                Vec3d pt = domain_low + dx * Vec3d(i,j,k);
                double dist = mag( pt - sphere_centre );
                phi(i,j,k) = dist - sphere_radius;
            }
        }
    }
}

// ---------------------------------------------------------
///
/// Create a sphere surface
///
// ---------------------------------------------------------

void create_sphere( const Vec3d& sphere_centre,
                   double sphere_radius,
                   double dx,
                   std::vector<Vec3d>& verts, 
                   std::vector<Vec3st>& tris )
{
    
    const Vec3d domain_low = sphere_centre - Vec3d(sphere_radius + 3*dx);
    const Vec3d domain_high = sphere_centre + Vec3d(sphere_radius + 3*dx);
    Array3d phi;
    create_sphere_signed_distance( sphere_centre, sphere_radius, dx, domain_low, domain_high, phi );  
    
    MarchingTilesHiRes marching_tiles( domain_low, dx, phi );
    marching_tiles.contour();
    marching_tiles.improve_mesh();
    
    std::vector<Vec3d> new_verts( marching_tiles.x.size() );
    for ( size_t i = 0; i < new_verts.size(); ++i )
    {
        new_verts[i] = marching_tiles.x[i];
    }
    
    std::vector<Vec3st> new_tris( marching_tiles.tri.size() );
    for ( size_t i = 0; i < new_tris.size(); ++i )
    {
        new_tris[i] = Vec3st( marching_tiles.tri[i] );
    }
    
    project_to_exact_sphere( new_verts, sphere_centre, sphere_radius );
    
    Vec3st offset( verts.size() );
    for ( size_t i = 0; i < new_verts.size(); ++i )
    {
        verts.push_back( new_verts[i] );
    }
    
    for ( size_t i = 0; i < new_tris.size(); ++i )
    {
        tris.push_back( new_tris[i] + offset );
    }
    
}

// ---------------------------------------------------------

void create_icosohedron( std::vector<Vec3d>& verts, 
                        std::vector<Vec3st>& tris )
{
    
    double t = 0.5 * ( 1. - sqrt(5.) );
    double c = 1. / sqrt( 1. - t*t );
    
    Vec3st offset( verts.size() );
    
    verts.push_back( c * Vec3d(  t,  1,  0 ) );
    verts.push_back( c * Vec3d( -t,  1,  0 ) );
    verts.push_back( c * Vec3d(  t, -1,  0 ) );
    verts.push_back( c * Vec3d( -t, -1,  0 ) );
    verts.push_back( c * Vec3d(  1,  0,  t ) );
    verts.push_back( c * Vec3d(  1,  0, -t ) );   
    verts.push_back( c * Vec3d( -1,  0,  t ) );
    verts.push_back( c * Vec3d( -1,  0, -t ) );
    verts.push_back( c * Vec3d(  0,  t,  1 ) );
    verts.push_back( c * Vec3d(  0, -t,  1 ) );
    verts.push_back( c * Vec3d(  0,  t, -1 ) );
    verts.push_back( c * Vec3d(  0, -t, -1 ) );
    
    tris.push_back( Vec3st(0,8,4) + offset );
    tris.push_back( Vec3st(0,5,10) + offset );
    tris.push_back( Vec3st(2,4,9) + offset );
    tris.push_back( Vec3st(2,11,5) + offset );
    tris.push_back( Vec3st(1,6,8) + offset );
    tris.push_back( Vec3st(1,10,7) + offset );
    tris.push_back( Vec3st(3,9,6) + offset );   
    tris.push_back( Vec3st(3,7,11) + offset );   
    tris.push_back( Vec3st(0,10,8) + offset );   
    tris.push_back( Vec3st(1,8,10) + offset );   
    tris.push_back( Vec3st(2,9,11) + offset );   
    tris.push_back( Vec3st(3,9,11) + offset );      
    tris.push_back( Vec3st(4,2,0) + offset );      
    tris.push_back( Vec3st(5,0,2) + offset );      
    tris.push_back( Vec3st(6,1,3) + offset );      
    tris.push_back( Vec3st(7,3,1) + offset );      
    tris.push_back( Vec3st(8,6,4) + offset );      
    tris.push_back( Vec3st(9,4,6) + offset );      
    tris.push_back( Vec3st(10,5,7) + offset );      
    tris.push_back( Vec3st(11,7,5) + offset );         
    
}


// ---------------------------------------------------------
///
/// Analytic entropy solution to motion in the normal direction of two spheres
///
// ---------------------------------------------------------

double signed_distance_entropy(const Vec3d& x,
                               const Vec3d& a,
                               const Vec3d& b,
                               double max_radius,
                               double interior_radius)
{
    double d=dist(a, b);
    double alpha=dot(x-a, b-a)/d; // distance along vector from a to b of projection of x onto line
    
    if(alpha<=0){ // x is on the far side of a from b
        return dist(x, a) - max_radius + interior_radius; // regular sphere distance
        
    }else if(alpha>=d){ // x is on the far side of b from a
        return dist(x, b) - max_radius + interior_radius; // regular sphere distance
        
    }else if(alpha<=0.5*d){ // x is between a and b, but closer to a
        if(max_radius<=0.5*d){ // if the spheres at maximum radius never intersected, life is simple
            return dist(x, a) - max_radius + interior_radius; // regular sphere distance
        }
        double beta=std::sqrt(dist2(x, a) - sqr(alpha)); // distance between x and projection onto line through a and b
        double gamma=std::sqrt(sqr(max_radius) - sqr(0.5*d)); // radius of intersection curve between spheres at max_radiu
        if(beta/alpha>=gamma/(0.5*d)){ // if closest point is still on a's sphere (not in the intersection region)
            return dist(x, a) - max_radius + interior_radius; // regular sphere distance
        }else{
            // use the distance from closest point on the intersection circle of the two spheres
            return interior_radius - sqrt(sqr(gamma-beta) + sqr(0.5*d-alpha));
        }
        
    }else{ // x is between a and b, but closer to b
        if(max_radius<=0.5*d){ // if the spheres at maximum radius never intersected, life is simple
            return dist(x, b) - max_radius + interior_radius; // regular sphere distance
        }
        double beta=std::sqrt(dist2(x, a) - sqr(alpha)); // distance between x and projection onto line through a and b
        double gamma=std::sqrt(sqr(max_radius) - sqr(0.5*d)); // radius of intersection curve between spheres at max_radiu
        if(beta/(d-alpha)>=gamma/(0.5*d)){ // if closest point is still on b's sphere (not in the intersection region)
            return dist(x, b) - max_radius + interior_radius; // regular sphere distance
        }else{
            // use the distance from closest point on the intersection circle of the two spheres
            return interior_radius - sqrt(sqr(gamma-beta) + sqr(alpha-0.5*d));
            
        }
    }
}


// ---------------------------------------------------------
///
/// Create signed distance field for two spheres
///
// ---------------------------------------------------------

void create_two_sphere_signed_distance( const Vec3d& sphere_a_centre, 
                                       const Vec3d& sphere_b_centre, 
                                       double sphere_radius, 
                                       double dx,
                                       const Vec3d& domain_low,
                                       const Vec3d& domain_high,                                     
                                       Array3d& phi )
{  
    phi.resize( (int) ceil( (domain_high[0]-domain_low[0]) / dx), (int) ceil( (domain_high[1]-domain_low[1]) / dx), (int) ceil( (domain_high[2]-domain_low[2]) / dx) );
    
    std::cout << "Generating signed distance function.  Grid resolution: " << phi.ni << " x " << phi.nj << " x " << phi.nk << std::endl;
    
    for ( int i = 0; i < phi.ni; ++i )
    {
        for ( int j = 0; j < phi.nj; ++j )
        {
            for ( int k = 0; k < phi.nk; ++k )
            {
                Vec3d pt = domain_low + dx * Vec3d(i,j,k);
                phi(i,j,k) = signed_distance_entropy( pt, sphere_a_centre, sphere_b_centre, sphere_radius, 0.0 );
            }
        }
    }
    
}


// ---------------------------------------------------------
///
/// Analytic signed distance function for a dumbbell
///
// ---------------------------------------------------------

double signed_distance_dumbbell( const Vec3d& pt,
                                const Vec3d& sphere_a_centre, 
                                const Vec3d& sphere_b_centre, 
                                double sphere_radius, 
                                double handle_width )
{
    double dist_spheres = ( min( mag(pt - sphere_a_centre), mag(pt - sphere_b_centre) ) - sphere_radius );
    return min( dist_spheres, max( fabs(pt[0]) - sphere_b_centre[0], sqrt( pt[1]*pt[1] + pt[2]*pt[2] ) - handle_width ) );   
}


// ---------------------------------------------------------
///
/// Create signed distance field for a dumbbell
///
// ---------------------------------------------------------

void create_dumbbell_signed_distance( const Vec3d& sphere_a_centre, 
                                     const Vec3d& sphere_b_centre, 
                                     double sphere_radius, 
                                     double handle_width,
                                     double dx,
                                     const Vec3d& domain_low,
                                     const Vec3d& domain_high,                                     
                                     Array3d& phi )
{
    
    phi.resize( (int) ceil( (domain_high[0]-domain_low[0]) / dx), (int) ceil( (domain_high[1]-domain_low[1]) / dx), (int) ceil( (domain_high[2]-domain_low[2]) / dx) );
    
    std::cout << "Generating signed distance function.  Grid resolution: " << phi.ni << " x " << phi.nj << " x " << phi.nk << std::endl;
    
    for ( int i = 0; i < phi.ni; ++i )
    {
        for ( int j = 0; j < phi.nj; ++j )
        {
            for ( int k = 0; k < phi.nk; ++k )
            {
                Vec3d pt = domain_low + dx * Vec3d(i,j,k);
                phi(i,j,k) = signed_distance_dumbbell( pt, sphere_a_centre, sphere_b_centre, sphere_radius, handle_width ); 
                
            }
        }
    }
}


// ---------------------------------------------------------
// Member function definitions
// ---------------------------------------------------------


