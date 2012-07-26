// ---------------------------------------------------------
//
//  meshmerger.h
//  Tyson Brochu 2011
//  
//  Search for mesh edges which are near to each other, zipper their neighbouring triangles together.
//
// ---------------------------------------------------------

#ifndef EL_TOPO_MESHMERGER_H
#define EL_TOPO_MESHMERGER_H

// ---------------------------------------------------------
//  Nested includes
// ---------------------------------------------------------

#include <meshoperator.h>
#include <vector>
#include <vec.h>

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

class SurfTrack;
class MeshMergeObserver;

// ---------------------------------------------------------
//  Class definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Pair of proximal edges, sortable by distance.  Used to build a list of edge pairs in ascending order of proximity, so we can 
/// handle them from nearest to farthest.
///
// ---------------------------------------------------------

struct SortableEdgeEdgeProximity
{
    SortableEdgeEdgeProximity( size_t a, size_t b, double d ) :
    edge_a( a ),
    edge_b( b ),
    distance( d )
    {}
    
    size_t edge_a;
    size_t edge_b;
    double distance;
    
    bool operator<( const SortableEdgeEdgeProximity& other ) const
    {
        return distance < other.distance;
    }
};

// ---------------------------------------------------------
///
///
///
// ---------------------------------------------------------

class MeshMerger : public MeshOperator
{
    
public:
    
    class PreMeshMergeInfo
    {
    public:
        PreMeshMergeInfo(size_t edge_index_a, 
                         size_t edge_index_b, 
                         const std::vector<size_t>& triangles_to_delete, 
                         const std::vector<Vec3st>& triangles_to_add );
        
        size_t m_edge_index_a;
        size_t m_edge_index_b;
        std::vector<size_t> m_triangles_to_delete;
        std::vector<Vec3st> m_triangles_to_add;
    };

    class PostMeshMergeInfo
    {
    public:
        PostMeshMergeInfo( const PreMeshMergeInfo& pre_merge_info );
        
        PreMeshMergeInfo m_pre_merge_info;
    };

    
    MeshMerger( SurfTrack& surf ) :
      MeshOperator( surf )
    {}
    
    /// 
    ///
    void process_mesh();
    
    void add_observer( MeshMergeObserver* observer );
    
private:
    
    // zipper
    void twist_vertices( size_t *zipper_vertices );
    bool get_zipper_triangles( size_t edge_index_0, size_t edge_index_1, std::vector<Vec3st>& output_triangles );
    bool zippering_introduces_collision( const std::vector<Vec3st>& new_triangles, const std::vector<size_t>& deleted_triangles );
    bool get_vertex_triangle_zipper_triangles( size_t v, size_t t, std::vector<Vec3st>& new_triangles );   
    bool stitch_triangle_triangle( size_t ta, size_t tb, std::vector<Vec3st>& new_tris );
    
    bool zipper_edges( size_t edge_index_a, size_t edge_index_b );
    
    std::vector<MeshMergeObserver*> m_observers;
    
};


// ---------------------------------------------------------
///
/// 
///
// ---------------------------------------------------------

class MeshMergeObserver
{
public:
    
    virtual ~MeshMergeObserver() {}
    
    virtual bool operationOK( const SurfTrack& /*surf*/, const MeshMerger::PreMeshMergeInfo& /*info*/ )
    {
        return true;
    }
    
    virtual void operationOccurred( const SurfTrack& /*surf*/, const MeshMerger::PostMeshMergeInfo& /*info*/ )
    {
    }
    
};



#endif
