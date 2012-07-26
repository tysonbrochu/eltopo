// ---------------------------------------------------------
//
//  meshrenderer.h
//  Tyson Brochu 2011
//  
//  OpenGL rendering for a triangle mesh.
//
// ---------------------------------------------------------

#ifndef EL_TOPO_MESHRENDERER_H
#define EL_TOPO_MESHRENDERER_H

// ---------------------------------------------------------
//  Nested includes
// ---------------------------------------------------------

#include <vec.h>

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

class DynamicSurface;

// ---------------------------------------------------------
//  Class definitions
// ---------------------------------------------------------

class MeshRenderer
{
    
public:
    
    MeshRenderer() :
    render_edges( true ),
    render_fill_triangles( true ),
    render_vertex_rank( false ),
    smooth_shading( true ),
    two_sided( true )
    {}
    
    bool render_edges;
    bool render_fill_triangles;
    bool render_vertex_rank;
    bool smooth_shading;
    bool two_sided;
    
    /// Display the surface in OpenGL using the current options settings
    ///
    void render( const DynamicSurface& surface );
    
    void render( const std::vector<Vec3d>& xs,
                const std::vector<Vec3d>& normals,
                const std::vector<Vec3st>& triangles,
                const std::vector<Vec2st>& edges );
    
    
};


#endif
