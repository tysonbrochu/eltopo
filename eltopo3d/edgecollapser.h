// ---------------------------------------------------------
//
//  edgecollapser.h
//  Tyson Brochu 2011
//  
//  Functions supporting the "edge collapse" operation: removing short edges from the mesh.
//
// ---------------------------------------------------------

#ifndef EL_TOPO_EDGECOLLAPSER_H
#define EL_TOPO_EDGECOLLAPSER_H

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
class EdgeCollapseObserver;

// ---------------------------------------------------------
//  Class definitions
// ---------------------------------------------------------

class EdgeCollapser : public MeshOperator
{
    
public:
    
    class PreEdgeCollapseInfo
    {
    public:
        PreEdgeCollapseInfo( size_t edge, size_t vertex_to_keep, size_t vertex_to_delete, const Vec3d& vertex_new_position );
        size_t m_edge;
        size_t m_vertex_to_keep;
        size_t m_vertex_to_delete;
        Vec3d m_vertex_new_position;
    };

    class PostEdgeCollapseInfo
    {
    public:
        PostEdgeCollapseInfo( const PreEdgeCollapseInfo& pre_collapse_info,
                              const std::vector<size_t>& triangles_incident_to_edge,
                              const std::vector<size_t>& triangles_incident_to_deleted_vertex,
                              const std::vector<size_t>& new_triangles,
                              const std::vector<size_t>& edges_incident_to_deleted_vertex );

        PreEdgeCollapseInfo m_pre_collapse_info;
        std::vector< size_t > m_triangles_incident_to_edge;
        std::vector< size_t > m_triangles_incident_to_deleted_vertex;
        std::vector< size_t > m_new_triangles;
        std::vector< size_t > m_edges_incident_to_deleted_vertex;
        
    };

    
    EdgeCollapser( SurfTrack& surf, bool use_curvature, double min_curvature_multiplier );
    
    bool small_angle_collapse_pass();
    
    void process_mesh();
    
    void add_observer( EdgeCollapseObserver* observer );
    
    /// Mimimum edge length.  Edges shorter than this will be collapsed.
    double m_min_edge_length;   
    
    bool m_use_curvature;
    double m_min_curvature_multiplier;
    
    
private:
    
    friend class SurfTrack;
    
    bool edge_is_collapsible( size_t edge_index );
    
    void get_moving_triangles( size_t source_vertex, 
                              size_t destination_vertex, 
                              size_t edge_index,
                              std::vector<size_t>& moving_triangles );
    
    void get_moving_edges( size_t source_vertex, 
                          size_t destination_vertex, 
                          size_t edge_index,
                          std::vector<size_t>& moving_edges );
    
    bool collapse_edge_pseudo_motion_introduces_collision( size_t source_vertex, 
                                                          size_t destination_vertex, 
                                                          size_t edge_index, 
                                                          const Vec3d& vertex_new_position );
    
    bool collapse_edge_introduces_normal_inversion( size_t source_vertex, 
                                                   size_t destination_vertex, 
                                                   size_t edge_index, 
                                                   const Vec3d& vertex_new_position );
    
    bool collapse_edge_introduces_volume_change( size_t source_vertex, 
                                                size_t edge_index, 
                                                const Vec3d& vertex_new_position );   
    
    bool collapse_edge_introduces_bad_angle( size_t source_vertex, 
                                            size_t destination_vertex, 
                                            size_t edge_index, 
                                            const Vec3d& vertex_new_position );
    
    bool collapse_edge( size_t edge );
    
    
    ///
    ///
    std::vector<EdgeCollapseObserver*> m_observers;

};


// ---------------------------------------------------------
///
/// 
///
// ---------------------------------------------------------

class EdgeCollapseObserver
{
public:
    
    virtual ~EdgeCollapseObserver() {}
    
    virtual bool operationOK( const SurfTrack& /*surf*/, const EdgeCollapser::PreEdgeCollapseInfo& /*info*/ )
    {
        return true;
    }
    
    virtual void operationOccurred( const SurfTrack& /*surf*/, const EdgeCollapser::PostEdgeCollapseInfo& /*info*/ )
    {
    }
    
};

#endif


