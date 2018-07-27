// ---------------------------------------------------------
//
//  iomesh.h
//
//  Non-member functions for reading and writing various mesh file formats.
//
// ---------------------------------------------------------

#ifndef EL_TOPO_IOMESH_H
#define EL_TOPO_IOMESH_H

// ---------------------------------------------------------
// Nested includes
// ---------------------------------------------------------

#include <fstream>
#include <vec.h>
#include <vector>

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

namespace ElTopo
{
class NonDestructiveTriMesh;
}

namespace Gluvi
{
    struct Target3D;
}

// ---------------------------------------------------------
//  Function declarations
// ---------------------------------------------------------

// ---------------------------------------------------------
//
// Read/write mesh in our own binary format
//
// ---------------------------------------------------------

bool write_binary_file( const ElTopo::NonDestructiveTriMesh &mesh,  const std::vector<ElTopo::Vec3d> &x, const std::vector<double> &masses, double curr_t, const char *filename_format, ...);
bool write_binary_file_with_velocities( const ElTopo::NonDestructiveTriMesh &mesh, const std::vector<ElTopo::Vec3d> &x, const std::vector<double> &masses, const std::vector<ElTopo::Vec3d> &v, double curr_t, const char *filename_format, ...);

bool write_surface_ids( const std::vector<size_t> &ids,
                       const char *filename_format, ... );

bool read_binary_file(ElTopo::NonDestructiveTriMesh &mesh, std::vector<ElTopo::Vec3d> &x, std::vector<double> &masses, double& curr_t, const char *filename_format, ...);
bool read_binary_file_with_velocities(ElTopo::NonDestructiveTriMesh &mesh, std::vector<ElTopo::Vec3d> &x, std::vector<double> &masses, std::vector<ElTopo::Vec3d> &v, double& curr_t, const char *filename_format, ...);

bool read_surface_ids( std::vector<unsigned int> &ids,
                      const char *filename_format, ... );


// ---------------------------------------------------------
//
// Read/write mesh in Wavefront OBJ format (ASCII)
//
// ---------------------------------------------------------

bool write_objfile(const ElTopo::NonDestructiveTriMesh &mesh, const std::vector<ElTopo::Vec3d> &x, const char *filename_format, ...);
bool read_objfile(ElTopo::NonDestructiveTriMesh &mesh, std::vector<ElTopo::Vec3d> &x, const char *filename_format, ...);

// ---------------------------------------------------------
//
// Write mesh in Renderman RIB format (geometry only)
//
// ---------------------------------------------------------

bool write_ribfile(const ElTopo::NonDestructiveTriMesh &mesh, const std::vector<float> &x, const char *filename_format, ...);
bool write_ribfile(const ElTopo::NonDestructiveTriMesh &mesh, const std::vector<float> &x, std::ostream &output);
bool write_ribfile(const ElTopo::NonDestructiveTriMesh &mesh, const std::vector<float> &x, FILE *output);

// ---------------------------------------------------------
///
/// Write an RIB file for the shadow map for the given light
///
// ---------------------------------------------------------

bool output_shadow_rib( Gluvi::Target3D& light, const std::vector<ElTopo::Vec3d>& positions,  const ElTopo::NonDestructiveTriMesh& mesh, const char *filename_format, ...);

// ---------------------------------------------------------
///
/// Write a render-ready RIB file.
///
// ---------------------------------------------------------

bool output_rib( const std::vector<ElTopo::Vec3d>& positions, const ElTopo::NonDestructiveTriMesh& mesh, const char *filename_format, ...);

// ---------------------------------------------------------
//
// Write mesh in PBRT format (geometry only)
//
// ---------------------------------------------------------

bool write_pbrtfile(const ElTopo::NonDestructiveTriMesh &mesh, const std::vector<float> &x, const char *filename_format, ...);
bool write_pbrtfile(const ElTopo::NonDestructiveTriMesh &mesh, const std::vector<float> &x, std::ostream &output);


// ---------------------------------------------------------
//
// Write an STL vector to an ASCII file.  Not really mesh-related, but useful.
//
// ---------------------------------------------------------

template<class T> void dump_vector_to_file( const char* filename, const std::vector<T, std::allocator<T> >& vec )
{
    std::ofstream outfile( filename, std::ios::out|std::ios::trunc );
    for ( unsigned int i = 0; i < vec.size(); ++i )
    {
        outfile << vec[i] << std::endl;
    }         
    outfile.close();
}


#endif
