// ---------------------------------------------------------
//
//  edgeflipper.h
//  Tyson Brochu 2011
//  
//  Functions supporting the "edge flip" operation: replacing non-delaunay edges with their dual edge.
//
// ---------------------------------------------------------

#ifndef EL_TOPO_EDGEFLIPPER_H
#define EL_TOPO_EDGEFLIPPER_H

// ---------------------------------------------------------
//  Nested includes
// ---------------------------------------------------------

#include <cstddef>
#include <meshoperator.h>
#include <vector>

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

class SurfTrack;
template<unsigned int N, class T> struct Vec;
typedef Vec<3,double> Vec3d;
typedef Vec<2,size_t> Vec2st;
typedef Vec<3,size_t> Vec3st;
class EdgeFlipObserver;

// ---------------------------------------------------------
//  Class definitions
// ---------------------------------------------------------

class EdgeFlipper : public MeshOperator
{
    
public:
    
    // ---------------------------------------------------------
    ///
    /// Information about the edge flip operation before it occurs 
    ///
    // ---------------------------------------------------------
    
    class PreEdgeFlipInfo
    {
    public:
        PreEdgeFlipInfo( size_t edge, size_t tri0, size_t tri1 );        
        size_t m_edge;    // index of edge to be flipped, before flipping
        size_t m_tri0;
        size_t m_tri1;
    };
    
    // ---------------------------------------------------------
    ///
    /// Information about the edge split operation after it has occurred
    ///
    // ---------------------------------------------------------
    
    class PostEdgeFlipInfo
    {
    public:
        PostEdgeFlipInfo(const PreEdgeFlipInfo& pre_split_info, size_t new_triangle_index_0, size_t new_triangle_index_1 ); 
        PreEdgeFlipInfo m_pre_flip_info;
        size_t m_new_triangle_index_0;
        size_t m_new_triangle_index_1;
    };
    

    EdgeFlipper( SurfTrack& surf, double edge_flip_min_length_change ) :
        MeshOperator( surf ),
        m_edge_flip_min_length_change( edge_flip_min_length_change )
    {}
    
    /// Flip all non-delaunay edges
    ///
    void process_mesh();
    
    void add_observer( EdgeFlipObserver* observer );
    
private:
    
    /// Minimum edge length improvement in order to flip an edge
    double m_edge_flip_min_length_change;
    
    ///
    ///
    bool flip_introduces_collision( size_t edge_index, 
                                   const Vec2st& new_edge, 
                                   const Vec3st& new_triangle_a, 
                                   const Vec3st& new_triangle_b );
    
    /// Flip an edge: remove the edge and its incident triangles, then add a new edge and two new triangles
    ///
    bool flip_edge(size_t edge, size_t tri0, size_t tri1, size_t third_vertex_0, size_t third_vertex_1 );
    
    size_t vertex_valence( size_t vertex_index );
    
    int total_mesh_regularity();
    
    bool regularity_improves( size_t edge, 
                             size_t tri0, 
                             size_t tri1, 
                             size_t third_vertex_0, 
                             size_t third_vertex_1 );
    
    ///
    ///
    std::vector<EdgeFlipObserver*> m_observers;
    
};

// ---------------------------------------------------------
///
/// 
///
// ---------------------------------------------------------

class EdgeFlipObserver
{
public:
    
    virtual ~EdgeFlipObserver() {}
    
    virtual bool operationOK( const SurfTrack& /*surf*/, const EdgeFlipper::PreEdgeFlipInfo& /*info*/ )
    {
        return true;
    }
    
    virtual void operationOccurred( const SurfTrack& /*surf*/, const EdgeFlipper::PostEdgeFlipInfo& /*info*/ )
    {
    }
    
};



#endif
