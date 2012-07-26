
// ---------------------------------------------------------
//
//  surftrack.h
//  Tyson Brochu 2008
//  
//  The SurfTrack class: a dynamic mesh with topological changes and mesh maintenance operations.
//
// ---------------------------------------------------------

#ifndef EL_TOPO_SURFTRACK_H
#define EL_TOPO_SURFTRACK_H

#include <dynamicsurface.h>
#include <edgecollapser.h>
#include <edgeflipper.h>
#include <edgesplitter.h>
#include <meshmerger.h>
#include <meshpincher.h>
#include <meshsmoother.h>

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

class SubdivisionScheme;
typedef std::vector<size_t> TriangleSet;
class DefragObserver;

// ---------------------------------------------------------
//  Class definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Structure for setting up a SurfTrack object with some initial parameters.  This is passed to the SurfTrack constructor.
///
// ---------------------------------------------------------

struct SurfTrackInitializationParameters
{
    
    /// Set default values for parameters which are not likely to be specified
    SurfTrackInitializationParameters();
    
    /// Elements closer than this are considered "near" (or proximate)
    double m_proximity_epsilon;
    
    double m_friction_coefficient;
    
    double m_min_triangle_area;
    
    /// Collision epsilon to use during mesh improvment operations (i.e. if any mesh elements are closer than this, the operation is 
    /// aborted).  NOTE: This should be greater than collision_epsilon, to prevent improvement operations from moving elements into 
    /// a collision configuration.
    double m_improve_collision_epsilon;
    
    /// Whether to set the min and max edge lengths as fractions of the initial average edge length
    bool m_use_fraction;
    
    /// If use_fraction is true, these are taken to be fractions of the average edge length of the new surface.
    /// If use_fraction is false, these are absolute.
    double m_min_edge_length;
    double m_max_edge_length; 
    double m_max_volume_change;
    
    // In-triangle angles to enforce
    double m_min_triangle_angle;
    double m_max_triangle_angle;   
    
    
    bool m_use_curvature_when_splitting;
    bool m_use_curvature_when_collapsing;
    
    // Clamp curvature scaling to these values
    double m_min_curvature_multiplier;
    double m_max_curvature_multiplier;
    
    bool m_allow_vertex_movement;
    
    /// Minimum edge length improvement in order to flip an edge
    double m_edge_flip_min_length_change;
    
    /// Elements within this distance will trigger a merge attempt   
    double m_merge_proximity_epsilon;
    
    /// Type of subdivision to use when collapsing or splitting (butterfly, quadric error minimization, etc.)
    SubdivisionScheme *m_subdivision_scheme;   
    
    /// Whether to enforce collision-free surfaces (including during mesh maintenance operations)
    bool m_collision_safety;
    
    /// Whether to allow changes in topology
    bool m_allow_topology_changes;
    
    /// Wether to allow non-manifold (edges incident on more than two triangles)
    bool m_allow_non_manifold;
    
    /// Whether to allow mesh improvement
    bool m_perform_improvement;
    
};

// ---------------------------------------------------------
///
/// Used to build a list of edges sorted in order of increasing length.
/// 
// ---------------------------------------------------------

struct SortableEdge
{
    size_t edge_index;
    double edge_length;
    
    SortableEdge( size_t ei, double el ) : edge_index(ei), edge_length(el) {}
    
    bool operator<( const SortableEdge& other ) const
    {
        return (this->edge_length < other.edge_length);
    }
};


// ---------------------------------------------------------
///
/// Keeps track of a vertex removal or addition.  If it's an addition, it also points to the edge that was split to create it.
///
// ---------------------------------------------------------

struct VertexUpdateEvent
{
    VertexUpdateEvent( bool _is_remove = false, size_t _vertex_index = (size_t)~0, const Vec2st& _split_edge = Vec2st((size_t)~0) ) :
    is_remove( _is_remove ),
    vertex_index( _vertex_index ),
    split_edge( _split_edge )
    {}
    
    static const bool VERTEX_REMOVE = true;
    static const bool VERTEX_ADD = false;
    
    bool is_remove;
    size_t vertex_index;   
    Vec2st split_edge;
    
};


// ---------------------------------------------------------
///
/// Keeps track of a triangle removal or addition. If addition, contains the three vertices that form the new triangle.
///
// ---------------------------------------------------------

struct TriangleUpdateEvent
{
    TriangleUpdateEvent( bool _is_remove = false, size_t _triangle_index = (size_t)~0, const Vec3st& _triangle = Vec3st((size_t)~0) ) :
    is_remove( _is_remove ),
    triangle_index( _triangle_index ),
    tri( _triangle )
    {}
    
    static const bool TRIANGLE_REMOVE = true;
    static const bool TRIANGLE_ADD = false;
    
    bool is_remove;
    size_t triangle_index;   
    Vec3st tri;
    
};


// ---------------------------------------------------------
///
/// A DynamicSurface with topological and mesh maintenance operations.
///
// ---------------------------------------------------------

class SurfTrack : public DynamicSurface
{
    
public:
    
    class PreDefragInfo
    {
        
    };
    
    class PostDefragInfo
    {
    public:
        std::vector<Vec2st> m_defragged_triangle_map;
        std::vector<Vec2st> m_defragged_vertex_map;
        std::vector<size_t> m_defragged_edge_map;
    };
    
    /// Create a SurfTrack object from a set of vertices and triangles using the specified paramaters
    ///
    SurfTrack( const std::vector<Vec3d>& vs, 
              const std::vector<Vec3st>& ts, 
              const std::vector<double>& masses,
              const SurfTrackInitializationParameters& initial_parameters );
    
    
    ~SurfTrack();
    
private:
    
    // Disallow copying and assignment by declaring private
    //
    SurfTrack( const SurfTrack& );
    SurfTrack& operator=( const SurfTrack& );
    
    
public:
    
    // ---------------------------------------------------------
    // Mesh bookkeeping
    
    size_t add_triangle(const Vec3st& t);
    
    void remove_triangle(size_t t);  
    
    size_t add_vertex( const Vec3d& new_vertex_position, double new_vertex_mass );
    
    void remove_vertex(size_t v);
    
    void defrag_mesh();
    
    
    // ---------------------------------------------------------
    // Main operations
    
    /// run mesh maintenance operations
    ///
    void improve_mesh( );
    
    /// run edge-edge merging
    ///
    void topology_changes( );
    
    
    // ---------------------------------------------------------
    // mesh cleanup
    // ---------------------------------------------------------
    
    /// Check for and delete flaps and zero-area triangles among the given triangle indices, then separate singular vertices.
    ///
    void trim_non_manifold( std::vector<size_t>& triangle_indices );
        
    /// Fire an assert if any degenerate triangles or tets (flaps) are found.
    /// 
    void assert_no_degenerate_triangles();
    
    
    void add_observer( DefragObserver* observer )
    {
        m_observers.push_back(observer);
    }
    
    
    // ---------------------------------------------------------
    // Member variables
    // ---------------------------------------------------------
    
    EdgeCollapser m_collapser;
    EdgeSplitter m_splitter;
    EdgeFlipper m_flipper;
    MeshSmoother m_smoother;
    
    MeshMerger m_merger;
    MeshPincher m_pincher;
    
    /// Collision epsilon to use during mesh improvment operations
    double m_improve_collision_epsilon;
    
    /// Minimum edge length improvement in order to flip an edge
    double m_edge_flip_min_length_change;
    
    /// Maximum volume change allowed when flipping or collapsing an edge
    double m_max_volume_change;
    
    /// Mimimum edge length.  Edges shorter than this will be collapsed.
    double m_min_edge_length;   
    
    /// Maximum edge length.  Edges longer than this will be subdivided.
    double m_max_edge_length;   
    
    /// Elements within this distance will trigger a merge attempt
    double m_merge_proximity_epsilon;
    
    /// Try to prevent triangles with area less than this
    double m_min_triangle_area;
    
    /// Don't create triangles with angles less than this.  If angles less than this do exist, try to remove them.
    double m_min_triangle_angle;
    
    /// Don't create triangles with angles greater than this.  If angles greater than this do exist, try to remove them.
    double m_max_triangle_angle;
    
    /// Interpolation scheme, determines edge midpoint location
    SubdivisionScheme *m_subdivision_scheme;
    
    /// If we allocate our own SubdivisionScheme object, we must delete it in this object's deconstructor.
    bool should_delete_subdivision_scheme_object;
    
    /// Triangles which are involved in connectivity changes which may introduce degeneracies
    std::vector<size_t> m_dirty_triangles;
    
    /// Whether to allow merging and separation
    bool m_allow_topology_changes;
    
    /// Allow more than two triangles incident on an edge
    bool m_allow_non_manifold;
    
    /// Whether to perform adaptivity operations
    bool m_perform_improvement;
    
    /// When doing mesh optimization, whether to allow the vertices to move.  If set to false, we allow edge flipping, edge 
    /// splitting, and edge collapsing (where the edge is collapsed down to one of its endpoints).  If true, we do mesh smoothing,
    /// as well as allowing a collapsed edge to collapse down to some point other than an endpoint.
    bool m_allow_vertex_movement;
    
    std::vector<VertexUpdateEvent> m_vertex_change_history;
    std::vector<TriangleUpdateEvent> m_triangle_change_history;
        
    std::vector<DefragObserver*> m_observers;
    
};

class DefragObserver
{
public:
    
    virtual ~DefragObserver() {}
    
    virtual bool operationOK( const SurfTrack& /*surf*/, const SurfTrack::PreDefragInfo& /*info*/ )
    {
        return true;
    }
    
    virtual void operationOccurred( const SurfTrack& /*surf*/, const SurfTrack::PostDefragInfo& /*info*/ )
    {
    }
    
};



#endif

