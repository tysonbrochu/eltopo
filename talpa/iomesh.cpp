
// ---------------------------------------------------------
//
//  iomesh.cpp
//
//  Non-member functions for reading and writing various mesh file formats.
//
// ---------------------------------------------------------


#include <iomesh.h>

#include <bfstream.h>
#include <cstdarg>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <nondestructivetrimesh.h>

#ifndef NO_GUI
#include <gluvi.h>
#endif

#define LINESIZE 1024 // maximum line size when reading .OBJ files


// ---------------------------------------------------------
///
/// Write mesh in binary format
///
// ---------------------------------------------------------

bool write_binary_file( const NonDestructiveTriMesh &mesh, 
                       const std::vector<Vec3d> &x,
                       const std::vector<double> &masses, 
                       double curr_t, 
                       const char *filename_format, ...)
{
    va_list ap;
    va_start(ap, filename_format);   
    bofstream outfile( filename_format, ap );
    va_end(ap);
    
    outfile.write_endianity();
    
    outfile << curr_t;
    
    unsigned int nverts = static_cast<unsigned int>( x.size() );
    outfile << nverts;
    for ( unsigned int i = 0; i < x.size(); ++i )
    {
        outfile << x[i][0];
        outfile << x[i][1];
        outfile << x[i][2];
    }
    
    assert( x.size() == masses.size() );
    
    for ( unsigned int i = 0; i < masses.size(); ++i )
    {
        outfile << masses[i];
    }
    
    unsigned int ntris = static_cast<unsigned int>( mesh.num_triangles() );
    outfile << ntris;
    
    for ( unsigned int t = 0; t < mesh.num_triangles(); ++t )
    {
        const Vec3ui& tri = Vec3ui( mesh.get_triangle(t) );
        outfile << tri[0];
        outfile << tri[1];
        outfile << tri[2];      
    }
    
    outfile.close();
    
    return outfile.good();
}

// ---------------------------------------------------------
///
/// Write mesh in binary format, with per-vertex velocities
///
// ---------------------------------------------------------

bool write_binary_file_with_velocities( const NonDestructiveTriMesh &mesh, 
                                       const std::vector<Vec3d> &x,
                                       const std::vector<double> &masses,                                       
                                       const std::vector<Vec3d> &v,
                                       double curr_t, 
                                       const char *filename_format, ...)
{
    va_list ap;
    va_start(ap, filename_format);   
    bofstream outfile( filename_format, ap );
    va_end(ap);
    
    outfile.write_endianity();
    
    outfile << curr_t;
    
    unsigned int nverts = static_cast<unsigned int>( x.size() );
    outfile << nverts;
    
    for ( unsigned int i = 0; i < x.size(); ++i )
    {
        outfile << x[i][0];
        outfile << x[i][1];
        outfile << x[i][2];
    }
    
    assert( x.size() == masses.size() );
    
    for ( unsigned int i = 0; i < masses.size(); ++i )
    {
        outfile << masses[i];
    }
    
    for ( unsigned int i = 0; i < v.size(); ++i )
    {
        outfile << v[i][0];
        outfile << v[i][1];
        outfile << v[i][2];
    }
    
    unsigned int ntris = static_cast<unsigned int>( mesh.num_triangles() );
    outfile << ntris;
    
    for ( unsigned int t = 0; t < mesh.num_triangles(); ++t )
    {
        const Vec3ui& tri = Vec3ui( mesh.get_triangle(t) );
        outfile << tri[0];
        outfile << tri[1];
        outfile << tri[2];      
    }
    
    outfile.close();
    
    return outfile.good();
}

// ---------------------------------------------------------

bool write_surface_ids( const std::vector<size_t> &ids,
                       const char *filename_format, ... )
{
    va_list ap;
    va_start(ap, filename_format);   
    bofstream outfile( filename_format, ap );
    va_end(ap);
    
    outfile.write_endianity();
    
    size_t nids = ids.size(); 
    outfile << nids;
    
    for ( unsigned int i = 0; i < ids.size(); ++i )
    {
        unsigned int curr_id = static_cast<unsigned int>(ids[i]);
        outfile << curr_id;
    }
    
    outfile.close();
    
    return outfile.good();
}


// ---------------------------------------------------------
///
/// Read mesh in binary format
///
// ---------------------------------------------------------

bool read_binary_file( NonDestructiveTriMesh &mesh, 
                      std::vector<Vec3d> &x, 
                      std::vector<double> &masses, 
                      double& curr_t, 
                      const char *filename_format, ...)
{
    va_list ap;
    va_start(ap, filename_format);   
    bifstream infile( filename_format, ap );
    va_end(ap);
    
    assert( infile.good() );
    
    infile.read_endianity();
    
    infile >> curr_t;
    
    unsigned int nverts;
    infile >> nverts;
    x.resize( nverts );
    for ( unsigned int i = 0; i < nverts; ++i )
    {
        infile >> x[i][0];
        infile >> x[i][1];
        infile >> x[i][2];  
        mesh.add_vertex();
    }
    
    masses.resize( nverts );
    for ( unsigned int i = 0; i < nverts; ++i )
    {
        infile >> masses[i];
    }
    
    unsigned int ntris;
    infile >> ntris;
    
    for ( unsigned int t = 0; t < ntris; ++t )
    {
        Vec3ui tri;
        infile >> tri[0];
        infile >> tri[1];
        infile >> tri[2];
        
        mesh.add_triangle( Vec3st(tri) );
    }
    
    infile.close();
    
    return infile.good();
}


// ---------------------------------------------------------
///
/// Read mesh in binary format, with per-vertex velocities
///
// ---------------------------------------------------------

bool read_binary_file_with_velocities( NonDestructiveTriMesh &mesh, 
                                      std::vector<Vec3d> &x, 
                                      std::vector<double> &masses,
                                      std::vector<Vec3d> &v, 
                                      double& curr_t, 
                                      const char *filename_format, ...)
{
    va_list ap;
    va_start(ap, filename_format);   
    bifstream infile( filename_format, ap );
    va_end(ap);
    
    infile.read_endianity();
    
    infile >> curr_t;
    
    unsigned int nverts;
    infile >> nverts;
    
    mesh.set_num_vertices(nverts);
    
    x.resize( nverts );
    for ( unsigned int i = 0; i < nverts; ++i )
    {
        infile >> x[i][0];
        infile >> x[i][1];
        infile >> x[i][2];      
    }
    
    masses.resize( nverts );
    for ( unsigned int i = 0; i < masses.size(); ++i )
    {
        infile >> masses[i];
    }
    
    v.resize( nverts );
    for ( unsigned int i = 0; i < nverts; ++i )
    {
        infile >> v[i][0];
        infile >> v[i][1];
        infile >> v[i][2];      
    }
    
    
    unsigned int ntris;
    infile >> ntris;
    
    for ( unsigned int t = 0; t < ntris; ++t )
    {
        Vec3st tri;
        infile >> tri[0];
        infile >> tri[1];
        infile >> tri[2];
        mesh.add_triangle(tri);
    }
    
    infile.close();
    
    return infile.good();
}


// ---------------------------------------------------------

bool read_surface_ids( std::vector<unsigned int>& ids,
                      const char *filename_format, ... )
{
    
    va_list ap;
    va_start(ap, filename_format);   
    bifstream infile( filename_format, ap );
    va_end(ap);
    
    infile.read_endianity();
    
    unsigned int n;
    infile >> n;
    
    ids.resize(n);
    for ( unsigned int i = 0; i < ids.size(); ++i )
    {
        infile >> ids[i];
    }
    
    infile.close();
    
    return infile.good();
    
}


// ---------------------------------------------------------
///
/// Write mesh in Wavefront OBJ format
///
// ---------------------------------------------------------

bool write_objfile(const NonDestructiveTriMesh &mesh, const std::vector<Vec3d> &x, const char *filename_format, ...)
{

#ifdef _MSC_VER
    va_list ap;
    va_start(ap, filename_format);
    int len=_vscprintf(filename_format, ap) // _vscprintf doesn't count
    +1; // terminating '\0'
    char *filename=new char[len];
    vsprintf(filename, filename_format, ap);
    std::ofstream output(filename, std::ofstream::binary);

    delete[] filename;
    va_end(ap);
#else
    va_list ap;
    va_start(ap, filename_format);
    char *filename;
    vasprintf(&filename, filename_format, ap);

    std::ofstream output(filename, std::ofstream::binary);
    std::free(filename);
    va_end(ap);
#endif    
    if(!output.good()) return false;
    
    output<<"# generated by editmesh"<<std::endl;
    for(unsigned int i=0; i<x.size(); ++i)
        output<<"v "<<x[i]<<std::endl;
    for(unsigned int t=0; t<mesh.num_triangles(); ++t)
        output<<"f "<<mesh.get_triangle(t)[0]+1<<' '<<mesh.get_triangle(t)[1]+1<<' '<<mesh.get_triangle(t)[2]+1<<std::endl; // correct for 1-based indexing in OBJ files
    return output.good();
}

namespace {
    
    // ---------------------------------------------------------
    ///
    /// Helper for reading OBJ file
    ///
    // ---------------------------------------------------------
    
    bool read_int(const char *s, int &value, bool &leading_slash, int &position)
    {
        leading_slash=false;
        for(position=0; s[position]!=0; ++position){
            switch(s[position]){
                case '/':
                    leading_slash=true;
                    break;
                case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7': case '8': case '9':
                    goto found_int;
            }
        }
        return false;
        
    found_int:
        value=0;
        for(;; ++position){
            switch(s[position]){
                case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7': case '8': case '9':
                    value=10*value+s[position]-'0';
                    break;
                default:
                    return true;
            }
        }
        return true; // should never get here, but keeps compiler happy
    }
    
    // ---------------------------------------------------------
    ///
    /// Helper for reading OBJ file
    ///
    // ---------------------------------------------------------
    
    void read_face_list(const char *s, std::vector<int> &vertex_list)
    {
        vertex_list.clear();
        int v, skip;
        bool leading_slash;
        for(int i=0;;){
            if(read_int(s+i, v, leading_slash, skip)){
                if(!leading_slash)
                    vertex_list.push_back(v-1); // correct for 1-based index
                i+=skip;
            }else
                break;
        }
    }
    
}  // unnamed namespace

// ---------------------------------------------------------
///
/// Read mesh in Wavefront OBJ format
///
// ---------------------------------------------------------

bool read_objfile(NonDestructiveTriMesh &mesh, std::vector<Vec3d> &x, const char *filename_format, ...)
{

#ifdef _MSC_VER
    va_list ap;
    va_start(ap, filename_format);
    int len=_vscprintf(filename_format, ap) // _vscprintf doesn't count
        +1; // terminating '\0'
    char *filename=new char[len];
    vsprintf(filename, filename_format, ap);
    std::ifstream input(filename, std::ifstream::binary);
    delete[] filename;
    va_end(ap);
#else
    va_list ap;
    va_start(ap, filename_format);
    char *filename;
    vasprintf(&filename, filename_format, ap);
    std::ifstream input(filename, std::ifstream::binary);
    std::free(filename);
    va_end(ap);
#endif

    if(!input.good()) return false;
    
    x.clear();
    mesh.clear();
    
    char line[LINESIZE];
    std::vector<int> vertex_list;
    while(input.good()){
        input.getline(line, LINESIZE);
        switch(line[0]){
            case 'v': // vertex data
                if(line[1]==' '){
                    Vec3d new_vertex;
                    std::sscanf(line+2, "%lf %lf %lf", &new_vertex[0], &new_vertex[1], &new_vertex[2]);
                    x.push_back(new_vertex);
                }
                break;
            case 'f': // face data
                if(line[1]==' '){
                    read_face_list(line+2, vertex_list);
                    for(int j=0; j<(int)vertex_list.size()-2; ++j)
                        mesh.add_triangle( Vec3st(vertex_list[0], vertex_list[j+1], vertex_list[j+2]) );
                }
                break;
        }
    }
    return true;
}

// ---------------------------------------------------------
///
/// Write mesh in Renderman RIB format.
///
// ---------------------------------------------------------

bool write_ribfile(const NonDestructiveTriMesh &mesh, const std::vector<float> &x, const char *filename_format, ...)
{
#ifdef _MSC_VER
    va_list ap;
    va_start(ap, filename_format);
    int len=_vscprintf(filename_format, ap) // _vscprintf doesn't count
        +1; // terminating '\0'
    char *filename=new char[len];
    vsprintf(filename, filename_format, ap);
    std::ofstream output(filename, std::ofstream::binary);
    delete[] filename;
    va_end(ap);
#else
    va_list ap;
    va_start(ap, filename_format);
    char *filename;
    vasprintf(&filename, filename_format, ap);
    std::ofstream output(filename, std::ofstream::binary);
    std::free(filename);
    va_end(ap);
#endif
    if(!output.good()) return false;
    return write_ribfile(mesh, x, output);
}

// ---------------------------------------------------------
///
/// Write mesh in Renderman RIB format.
///
// ---------------------------------------------------------

bool write_ribfile(const NonDestructiveTriMesh &mesh, const std::vector<float> &x, std::ostream &output)
{
    output<<"# generated by editmesh"<<std::endl;
    output<<"PointsPolygons"<<std::endl;
    output<<" [ ";
    const std::vector<Vec3st>& tris = mesh.get_triangles();
    for(unsigned int i=0; i<tris.size(); ++i){
        output<<"3 ";
        if(i%38==37 && i!=tris.size()-1) output<<std::endl;
    }
    output<<"]"<<std::endl;
    output<<" [ ";
    for(unsigned int i=0; i<tris.size(); ++i){
        output<<tris[i]<<"  ";
        if(i%6==5 && i!=tris.size()-1) output<<std::endl;
    }
    output<<"]"<<std::endl;
    output<<" \"P\" [";
    for(unsigned int i=0; i<x.size(); ++i){
        output<<x[i]<<"  ";
        if(i%4==3 && i!=x.size()-1) output<<std::endl;
    }
    output<<"]"<<std::endl;
    
    return output.good();
}

// ---------------------------------------------------------
///
/// Write mesh in Renderman RIB format.
///
// ---------------------------------------------------------

bool write_ribfile(const NonDestructiveTriMesh &mesh, const std::vector<float> &x, FILE *output)
{
    fprintf( output, "# generated by editmesh\n" );
    fprintf( output, "PointsPolygons\n" );
    fprintf( output, " [ " );
    const std::vector<Vec3st>& tris = mesh.get_triangles();
    for(unsigned int i=0; i<tris.size(); ++i){
        fprintf( output, "3 " );
        if(i%38==37 && i!=tris.size()-1) fprintf( output, "\n" );
    }
    fprintf( output, "]\n" );
    fprintf( output, " [ " );
    for(unsigned int i=0; i<tris.size(); ++i){
        Vec3ui new_tri;
        fprintf( output, " %d %d %d ", (int)tris[i][0], (int)tris[i][1], (int)tris[i][2] );
        if(i%6==5 && i!=tris.size()-1) fprintf( output, "\n" ); 
    }
    fprintf( output, "]\n" );
    fprintf( output, " \"P\" [" );
    for(unsigned int i=0; i<x.size(); ++i){
        fprintf( output, " %f ", x[i] );
        if(i%4==3 && i!=x.size()-1) fprintf( output, "\n" ); 
    }
    fprintf( output, "]\n" );
    
    return true; 
}


#ifndef NO_GUI

// ---------------------------------------------------------
///
/// Write an RIB file for the shadow map for the given light
///
// ---------------------------------------------------------

bool output_shadow_rib( Gluvi::Target3D& light, const std::vector<Vec3d>& positions, const NonDestructiveTriMesh& mesh, const char *filename_format, ...)
{
#ifdef _MSC_VER
    va_list ap;
    va_start(ap, filename_format);
    int len=_vscprintf(filename_format, ap) // _vscprintf doesn't count
        +1; // terminating '\0'
    char *filename=new char[len];
    vsprintf(filename, filename_format, ap);

    std::ofstream out;
    out.open( filename );

    delete[] filename;

    if( out == NULL )
    {
        return false;
    }

    len=_vscprintf("track%04d_shadow.tiff", ap) // _vscprintf doesn't count
        +1; // terminating '\0'
    filename=new char[len];
    vsprintf( filename, "track%04d_shadow.tiff", ap );

    va_end(ap);
#else
    va_list ap;
    va_start(ap, filename_format);
    char *filename;
    vasprintf(&filename, filename_format, ap);   
    
    std::ofstream out;
    out.open( filename );
    
    delete[] filename;
    
    if( out == NULL )
    {
        return false;
    }
    
    vasprintf( &filename, "track%04d_shadow.tiff", ap );
    
    va_end(ap);
#endif
    // flatten
    std::vector<float> xs;
    for ( unsigned int i = 0; i < positions.size(); ++i )
    {
        xs.push_back( (float) positions[i][0] );
        xs.push_back( (float) positions[i][1] );
        xs.push_back( (float) positions[i][2] );
    }
    
    out << "Display \"" << filename << "\" \"file\" \"z\"" << std::endl;
    delete[] filename;
    
    // next line: image format (width and height in pixels, pixel aspect ratio)
    out << "Format " << "1024 1024 1" << std::endl;
    out << "PixelFilter \"box\" 1 1 " << std::endl;
    
    // then write out the camera specification
    light.export_rib(out);
    
    // start the scene
    out << "WorldBegin\n";
    
    out << "AttributeBegin\n";
    out << "  Color [0.6 0.6 0.2]\n";
    out << "  Opacity [1 1 1]\n";
    out << "  Surface \"matte\" \"Kd\" 1\n";
    
    const float plane_limit = 50.0f;
    char buf[256];
    sprintf( buf, "  Polygon \"P\" [-%f -%f -%f %f -%f -%f  %f %f -%f  -%f %f -%f]\n", 
            plane_limit, plane_limit, plane_limit, plane_limit, plane_limit, plane_limit, plane_limit, plane_limit, plane_limit, plane_limit, plane_limit, plane_limit );
    out << buf;
    
    out << "AttributeEnd\n";
    
    out << "Color [0.7 0.7 0.9]\n";
    out << "Surface \"matte\" \n";
    
    write_ribfile( mesh, xs, out);
    
    // finish the scene
    out << "WorldEnd\n";
    
    out.flush();
    out.close();
    
    return true;
}   


// ---------------------------------------------------------
///
/// Write a render-ready RIB file.
///
// ---------------------------------------------------------

bool output_rib( const std::vector<Vec3d>& positions, const NonDestructiveTriMesh& mesh, const char *filename_format, ...)
{
#ifdef _MSC_VER
    va_list ap;
    va_start(ap, filename_format);
    int len=_vscprintf(filename_format, ap) // _vscprintf doesn't count
        +1; // terminating '\0'
    char *filename=new char[len];
    vsprintf(filename, filename_format, ap);

    std::ofstream out;
    out.open( filename );

    delete[] filename;

    if( out == NULL )
    {
        return false;
    }

    // first line: what image file this RIB file should produce
    len=_vscprintf("track%04d.tiff", ap) // _vscprintf doesn't count
        +1; // terminating '\0'
    filename=new char[len];
    vsprintf( filename, "track%04d.tiff", ap );

    len=_vscprintf("track%04d_shadow.shad", ap) // _vscprintf doesn't count
        +1; // terminating '\0'
    char *shadow_filename=new char[len];
    vsprintf( shadow_filename, "track%04d_shadow.shad", ap );
    
    va_end(ap);
#else
    va_list ap;
    va_start(ap, filename_format);
    char *filename;
    vasprintf(&filename, filename_format, ap);   
    
    std::ofstream out;
    out.open( filename );
    
    delete[] filename;
    
    if( out == NULL )
    {
        return false;
    }
    
    // first line: what image file this RIB file should produce
    vasprintf( &filename, "track%04d.tiff", ap );
    
    char *shadow_filename;
    vasprintf( &shadow_filename, "track%04d_shadow.shad", ap );
    
    va_end(ap);
#endif
    // flatten
    std::vector<float> xs;
    for ( unsigned int i = 0; i < positions.size(); ++i )
    {
        xs.push_back( (float) positions[i][0] );
        xs.push_back( (float) positions[i][1] );
        xs.push_back( (float) positions[i][2] );
    }
    
    std::vector<Vec3f> normals;
    for ( unsigned int i = 0; i < positions.size(); ++i )
    {
        Vec3f n(0,0,0);
        for ( unsigned int j = 0; j < mesh.m_vertex_to_triangle_map[i].size(); ++j )
        {
            const Vec3st& tri = mesh.get_triangle( mesh.m_vertex_to_triangle_map[i][j] );
            Vec3d u = positions[ tri[1] ] - positions[ tri[0] ];
            Vec3d v = positions[ tri[2] ] - positions[ tri[0] ];
            Vec3d tn = normalized(cross(u, v));          
            n += Vec3f( (float)tn[0], (float)tn[1], (float)tn[2] ); 
        }
        normals.push_back( n / (float) mesh.m_vertex_to_triangle_map[i].size() );
    }
    
    
    out << "Display \"" << filename << "\" \"file\" \"rgb\"" << std::endl;
    delete[] filename;
    
    // next line: image format (width and height in pixels, pixel aspect ratio)
    out << "Format " << Gluvi::winwidth << " " << Gluvi::winheight << " 1" << std::endl;
    out << "PixelSamples 2 2" << std::endl;
    out << "Exposure 1.0 2.2" << std::endl;
    
    // then write out the camera specification
    Gluvi::camera->export_rib(out);
    
    // start the scene
    out << "WorldBegin\n";
    
    out << "LightSource \"ambientlight\" 1 \"intensity\" .3 \"lightcolor\" [1 1 1]\n";
    
    /*
     for ( unsigned int i = 0; i < lights.size(); ++i )
     {
     // compute location of light
     Vec3f light_pos( 0, 0, lights[i].dist );
     rotate( light_pos, lights[i].pitch, Vec3f(1,0,0) );
     rotate( light_pos, lights[i].heading, Vec3f(0,1,0) );
     light_pos += Vec3f( lights[i].target[0], lights[i].target[1], lights[i].target[2] );
     
     out << "LightSource \"singleshadowpoint\" 2 \"intensity\" 30 \"from\" [" << light_pos << "] \"shadowmap\" \"" << shadow_filename << "\"\n";
     }
     */
    
    //out << "LightSource \"distantlight\" 3 \"intensity\" 0.3 \"from\" [-5 -10 20] \"to\" [0 0 0]\n";
    
    out << "AttributeBegin\n";
    out << "  Color [0.6 0.6 0.2]\n";
    out << "  Opacity [1 1 1]\n";
    out << "  Surface \"matte\" \"Kd\" 1\n";
    
    const float plane_limit = 50.0f;
    const float plane_distance = 10.0f;
    char buf[256];
    sprintf( buf, "  Polygon \"P\" [-%f -%f -%f %f -%f -%f  %f %f -%f  -%f %f -%f]\n", 
            plane_limit, plane_limit, plane_distance, 
            plane_limit, plane_limit, plane_distance, 
            plane_limit, plane_limit, plane_distance, 
            plane_limit, plane_limit, plane_distance );
    out << buf;
    
    out << "AttributeEnd\n";
    
    out << "Color [0.3 0.3 0.9]\n";
    out << "Surface \"matte\" \n";
    
    write_ribfile( mesh, xs, out);
    
    out << " \"N\" [";
    for(unsigned int i=0; i<normals.size(); ++i)
    {
        out << normals[i] << "  ";
    }
    out << "]" << std::endl;
    
    // finish the scene
    out << "WorldEnd\n";
    
    out.flush();
    out.close();
    
    return true;
}

#endif

// ---------------------------------------------------------
//
// Write mesh in PBRT format
//
// ---------------------------------------------------------

bool write_pbrtfile(const NonDestructiveTriMesh &mesh, const std::vector<float> &x, const char *filename_format, ...)
{

#ifdef _MSC_VER
    va_list ap;
    va_start(ap, filename_format);
    int len=_vscprintf(filename_format, ap) // _vscprintf doesn't count
        +1; // terminating '\0'
    char *filename=new char[len];
    vsprintf(filename, filename_format, ap);
    std::ofstream output(filename, std::ofstream::binary);
    delete[] filename;
    va_end(ap);
#else
    va_list ap;
    va_start(ap, filename_format);
    char *filename;
    vasprintf(&filename, filename_format, ap);
    std::ofstream output(filename, std::ofstream::binary);
    std::free(filename);
    va_end(ap);
#endif
    if(!output.good()) return false;
    return write_ribfile(mesh, x, output);
}

// ---------------------------------------------------------
//
// Write mesh in PBRT format
//
// ---------------------------------------------------------

bool write_pbrtfile(const NonDestructiveTriMesh &mesh, const std::vector<float> &x, std::ostream &output)
{
    output<<"# generated by editmesh"<<std::endl;
    
    //output<<"\"integer nlevels\" [3]"<<std::endl;
    output<<"\"point P\" ["<<std::endl;
    for(unsigned int i=0; i<x.size(); ++i){
        output<<x[i]<<"  ";
        if(i%4==3 && i!=x.size()-1) output<<std::endl;
    }
    output<<"]"<<std::endl;
    output<<" \"integer indices\" ["<<std::endl;
    const std::vector<Vec3st>& tris = mesh.get_triangles();
    for(unsigned int i=0; i<tris.size(); ++i){
        output<<tris[i]<<"  ";
        if(i%6==5 && i!=tris.size()-1) output<<std::endl;
    }
    output<<"]"<<std::endl;
    
    return output.good();
}





