// ---------------------------------------------------------
//
//  geometryinit.h
//  Tyson Brochu 2008
//
//  A set of geometry initialization functions.
//
// ---------------------------------------------------------

#ifndef EL_TOPO_GEOMETRYINIT_H
#define EL_TOPO_GEOMETRYINIT_H

// ---------------------------------------------------------
// Nested includes
// ---------------------------------------------------------

#include <vec.h>

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

class SurfTrack;

// ---------------------------------------------------------
//  Function declarations
// ---------------------------------------------------------

void append_mesh( std::vector<Vec3st>& tris, 
                 std::vector<Vec3d>& verts,
                 std::vector<double>& masses,
                 const std::vector<Vec3st>& new_tris, 
                 const std::vector<Vec3d>& new_verts,
                 const std::vector<double>& new_masses );

void contour_phi( const Vec3d& domain_low, double domain_dx, Array3d& phi, std::vector<Vec3st>& tris, std::vector<Vec3d>& verts );

void create_circle( std::vector<Vec3d>& verts, 
                   std::vector<Vec3st>& tris,
                   std::vector<double>& masses,
                   const Vec3d& centre,
                   double radius,
                   size_t nx );

/// Create a vertical or horizontal, triangulated sheet.
///
void create_sheet( std::vector<Vec3d>& verts, 
                  std::vector<Vec3st>& tris, 
                  const Vec3d& low_corner,  
                   double width, 
                  size_t nx,
                  size_t ny );  

/// Create a vertical or horizontal, triangulated sheet.
///
void create_curved_sheet( std::vector<Vec3d>& verts, 
                         std::vector<Vec3st>& tris, 
                         const Vec3d& low_corner, 
                         double width, 
                         size_t nx,
                         size_t ny );

/// 
///
void project_to_exact_cube( std::vector<Vec3d>& verts, 
                           const Vec3d& cube_low, 
                           const Vec3d& cube_high );

/// Project mesh vertices onto an analytic sphere
///
void project_to_exact_sphere( std::vector<Vec3d>& verts, 
                             const Vec3d& sphere_centre, 
                             double radius );


/// Project mesh vertices onto analytic dumbbell
///
void project_to_exact_dumbbell( std::vector<Vec3d>& verts, 
                               const Vec3d& sphere_a_centre, 
                               const Vec3d sphere_b_centre, 
                               double sphere_radius, 
                               double handle_width );

/// Read a signed distance field from an ASCII file
///
void read_signed_distance( const char* filename, Array3d& signed_distance );

/// Create signed distance field for a capsule with the given geometry
///
void create_capsule_signed_distance( const Vec3d& capsule_end_a, 
                                    const Vec3d& capsule_end_b, 
                                    double capsule_radius,
                                    double dx,
                                    const Vec3d& domain_low,
                                    const Vec3d& domain_high,                                     
                                    Array3d& phi );

/// Create signed distance field for a cube with the given geometry
///
void create_cube_signed_distance( const Vec3d& cube_low, 
                                 const Vec3d& cube_high, 
                                 double dx,
                                 const Vec3d& domain_low,
                                 const Vec3d& domain_high,                                     
                                 Array3d& phi );

/// Create signed distance field for a sphere with the given geometry
///
void create_sphere_signed_distance( const Vec3d& sphere_centre, 
                                   double sphere_radius, 
                                   double dx,
                                   const Vec3d& domain_low,
                                   const Vec3d& domain_high,                                     
                                   Array3d& phi );

void create_icosohedron( std::vector<Vec3d>& verts, 
                        std::vector<Vec3st>& tris );

///
///
void create_sphere( const Vec3d& sphere_centre,
                   double sphere_radius,
                   double dx,
                   std::vector<Vec3d>& verts, 
                   std::vector<Vec3st>& tris );


///  Analytic entropy solution to motion in the normal direction of two spheres
///
double signed_distance_entropy( const Vec3d& pt,
                               const Vec3d& sphere_a_centre, 
                               const Vec3d& sphere_b_centre, 
                               double sphere_max_radius,
                               double sphere_interior_radius );

/// Create signed distance field for two spheres
///
void create_two_sphere_signed_distance( const Vec3d& sphere_a_centre, 
                                       const Vec3d& sphere_b_centre, 
                                       double sphere_radius, 
                                       double dx,
                                       const Vec3d& domain_low,
                                       const Vec3d& domain_high,                                     
                                       Array3d& phi );

/// Analytic signed distance function for a dumbbell
///
double signed_distance_dumbbell( const Vec3d& pt,
                                const Vec3d& sphere_a_centre, 
                                const Vec3d& sphere_b_centre, 
                                double sphere_radius, 
                                double handle_width );

/// Create signed distance field for a dumbbell
///
void create_dumbbell_signed_distance( const Vec3d& sphere_a_centre, 
                                     const Vec3d& sphere_b_centre, 
                                     double sphere_radius, 
                                     double handle_width,
                                     double dx,
                                     const Vec3d& domain_low,
                                     const Vec3d& domain_high,                                     
                                     Array3d& phi );

#endif


