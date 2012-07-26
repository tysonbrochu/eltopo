// ---------------------------------------------------------
//
//  meshpincher.h
//  Tyson Brochu 2011
//  
//  Look for a loop of three vertices that can be "pinched" to separate a volume into two.
//
// ---------------------------------------------------------

#ifndef EL_TOPO_MESHPINCHER_H
#define EL_TOPO_MESHPINCHER_H

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
class MeshPinchObserver;

// ---------------------------------------------------------
//  Class definitions
// ---------------------------------------------------------

class MeshPincher : public MeshOperator
{
    
public:
    
    typedef std::vector<size_t> TriangleSet;
    
    class PrePinchInfo
    {
    public:
        PrePinchInfo(size_t vertex_index, 
                     const std::vector<TriangleSet>& connected_components,
                     const TriangleSet& triangles_to_delete,
                     const std::vector<Vec3st>& triangles_to_add,
                     const std::vector<size_t>& vertices_added );
        
        size_t m_vertex_index;        
        const std::vector<TriangleSet>& m_connected_components;
        const TriangleSet& m_triangles_to_delete;
        const std::vector<Vec3st>& m_triangles_to_add;
        const std::vector<size_t>& m_vertices_added;
    };
    
    class PostPinchInfo
    {
    public:
        PostPinchInfo( PrePinchInfo& pre_pinch_info );
        PrePinchInfo m_pre_pinch_info;
    };
    
    MeshPincher( SurfTrack& surf ) :
        MeshOperator( surf )
    {}
    
    /// Find vertices which connect two otherwise disjoint surface patches and pull apart the mesh at those vertices.
    ///
    void process_mesh();
    
    void add_observer( MeshPinchObserver* observer );
        
private:
        
    void partition_vertex_neighbourhood( size_t vertex_index, std::vector< TriangleSet >& connected_components );
    
    bool pull_apart_vertex( size_t vertex_index, const std::vector< TriangleSet >& connected_components );
    
    std::vector<MeshPinchObserver*> m_observers;
        
};


// ---------------------------------------------------------
///
/// 
///
// ---------------------------------------------------------

class MeshPinchObserver
{
public:
    
    virtual ~MeshPinchObserver() {}
    
    virtual bool operationOK( const SurfTrack& /*surf*/, const MeshPincher::PrePinchInfo& /*info*/ )
    {
        return true;
    }
    
    virtual void operationOccurred( const SurfTrack& /*surf*/, const MeshPincher::PostPinchInfo& /*info*/ )
    {
    }
    
};


#endif
