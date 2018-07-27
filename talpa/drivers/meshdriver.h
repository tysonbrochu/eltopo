// ---------------------------------------------------------
//
//  meshdriver.h
//  Tyson Brochu 2008
//
//  Interface for a class that specifies vertex velocities over a mesh.
//
// ---------------------------------------------------------

#ifndef EL_TOPO_MESHDRIVER_H
#define EL_TOPO_MESHDRIVER_H

// ---------------------------------------------------------
// Nested includes
// ---------------------------------------------------------

#include <surftrack.h>
#include <vec.h>
#include <vector>

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

using ElTopo::SurfTrack;

// ---------------------------------------------------------
//  Class definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Base class for anything that moves a dynamic mesh
///
// ---------------------------------------------------------


class MeshDriver
{
public:
    
    class MeshDriverRenderer
    {
    public:
        virtual ~MeshDriverRenderer() {}
        virtual void render() const {}
    };
    
    virtual const MeshDriverRenderer* get_renderer()
    {
        return NULL;
    }
    
    
    /// Silence compiler warning
    ///
    virtual ~MeshDriver() {}
    
    /// Initialize the driver with the surface at t = 0.  Default does nothing
    ///
    virtual void initialize( const ElTopo::SurfTrack& ) {}
    
    /// Draw something with OpenGL.  Default does nothing.
    ///
    virtual void display() {}
    
    /// Set velocities on each mesh vertex
    ///
    virtual void set_predicted_vertex_positions(
        const ElTopo::SurfTrack &surf,
        std::vector<ElTopo::Vec3d> &predicted_positions, double current_t,
        double &adaptive_dt) = 0;

    /// Compute and output error.  Default does nothing.
    ///
    virtual void compute_error( const ElTopo::SurfTrack&, double /*current_t*/ ) {}
    
    ///
    ///
    virtual void notify_done_integration( const std::vector<ElTopo::Vec3d>& /*initial_positions*/,
                                         const std::vector<ElTopo::Vec3d>& /*final_positions*/,
                                         double /*actual_dt*/ ) 
    {}
    
    ///
    ///
    virtual void update_simulation_elements(ElTopo::SurfTrack& surf )
    {
        std::vector<ElTopo::VertexUpdateEvent>& vertex_history = surf.m_vertex_change_history;
        vertex_history.clear();
        
        std::vector<ElTopo::TriangleUpdateEvent>& triangle_history = surf.m_triangle_change_history;
        triangle_history.clear();        
    }
    
    virtual void write_to_disk( char* /*output_path*/, unsigned int /*frame_num*/ ) {}
    
};



#endif

