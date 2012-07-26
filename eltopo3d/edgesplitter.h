// ---------------------------------------------------------
//
//  edgesplitter.h
//  Tyson Brochu 2011
//  
//  Functions supporting the "edge split" operation: subdividing an edge into two shorter edges.
//
// ---------------------------------------------------------

#ifndef EL_TOPO_EDGESPLITTER_H
#define EL_TOPO_EDGESPLITTER_H

// ---------------------------------------------------------
//  Nested includes
// ---------------------------------------------------------

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

class EdgeSplitObserver;

// ---------------------------------------------------------
//  Class definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Splits long edges
///
// ---------------------------------------------------------

class EdgeSplitter : public MeshOperator
{
    
public:
    
    // ---------------------------------------------------------
    ///
    /// Information about the edge split operation before it occurs 
    ///
    // ---------------------------------------------------------
    
    class PreEdgeSplitInfo
    {
    public:
        PreEdgeSplitInfo( size_t edge, 
                          size_t tri0, 
                          size_t tri1,
                          size_t vertex_a, 
                          size_t vertex_b, 
                          size_t vertex_c, 
                          size_t vertex_d );
        
        size_t m_edge;    // index of edge to be split, before splitting
        size_t m_tri0;    // triangles incident on the edge to be split
        size_t m_tri1;
        size_t m_vertex_a;
        size_t m_vertex_b;
        size_t m_vertex_c;
        size_t m_vertex_d;
    };

    // ---------------------------------------------------------
    ///
    /// Information about the edge split operation after it has occurred
    ///
    // ---------------------------------------------------------

    class PostEdgeSplitInfo
    {
    public:
        PostEdgeSplitInfo(const PreEdgeSplitInfo& pre_split_info, 
                          size_t new_vertex, 
                          size_t new_tri0, 
                          size_t new_tri1, 
                          size_t new_tri2, 
                          size_t new_tri3 );
        
        PreEdgeSplitInfo m_pre_split_info;
        size_t m_new_vertex;
        size_t m_new_tri0, m_new_tri1, m_new_tri2, m_new_tri3;
    };

    
    EdgeSplitter( SurfTrack& surf, bool use_curvature, double max_curvature_multiplier );
        
        
    /// Split all long edges
    ///
    void process_mesh();
        
    ///
    ///
    void add_observer( EdgeSplitObserver* observer );

    bool m_use_curvature;
    
    ///
    double m_max_curvature_multiplier;
    
    /// Maximum edge length.  Edges longer than this will be subdivided.
    double m_max_edge_length;   
    
private:
    
    bool split_edge_edge_collision(size_t neighbour_index, 
                                   const Vec3d& new_vertex_position, 
                                   const Vec3d& new_vertex_smooth_position, 
                                   const Vec2st& edge );
    
    bool split_triangle_vertex_collision(const Vec3st& triangle_indices, 
                                         const Vec3d& new_vertex_position, 
                                         const Vec3d& new_vertex_smooth_position, 
                                         size_t overlapping_vert_index, 
                                         const Vec3d& vert );
    
    // helper
    bool split_edge_pseudo_motion_introduces_intersection(const Vec3d& new_vertex_position, 
                                                          const Vec3d& new_vertex_smooth_position, 
                                                          size_t edge,
                                                          size_t tri0, 
                                                          size_t tri1,
                                                          size_t vertex_a,
                                                          size_t vertex_b,
                                                          size_t vertex_c,
                                                          size_t vertex_d );
    
    bool edge_is_splittable( size_t edge_index );
    
    /// Split an edge, using subdivision_scheme to determine the new vertex location, if safe to do so.
    ///
    bool split_edge( size_t edge );
    
    bool large_angle_split_pass();

    ///
    ///
    std::vector<EdgeSplitObserver*> m_observers;
        
};

// ---------------------------------------------------------
///
/// 
///
// ---------------------------------------------------------

class EdgeSplitObserver
{
public:
    
    virtual ~EdgeSplitObserver() {}
    
    virtual bool operationOK( const SurfTrack& /*surf*/, const EdgeSplitter::PreEdgeSplitInfo& /*info*/ )
    {
        return true;
    }
    
    virtual void operationOccurred( const SurfTrack& /*surf*/, const EdgeSplitter::PostEdgeSplitInfo& /*info*/ )
    {
    }
    
};



#endif
