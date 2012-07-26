
// ---------------------------------------------------------
//
//  eltopo.h
//  Tyson Brochu 2009
//
//  C-callable API for El Topo
//
// ---------------------------------------------------------

#ifndef ELTOPO_H
#define ELTOPO_H

#ifdef __cplusplus
extern "C" {
#endif
    
    
    // =========================================================
    //  STRUCTURES FOR SPECIFYING OPTIONS
    // =========================================================   
    
    // ---------------------------------------------------------
    ///
    /// Options common to static operations and position integration.
    ///
    // ---------------------------------------------------------
    
    struct ElTopoGeneralOptions
    {
        int m_verbose;               // whether to output a lot of information to 
        // the console
        int m_collision_safety;      // whether to enforce an intersection-free mesh
        // at all times
        
        double m_proximity_epsilon;
        
    };
    
    // ---------------------------------------------------------
    ///
    /// Options for static operations (mesh maintenance and topology change)
    ///
    // ---------------------------------------------------------
    
    struct ElTopoStaticOperationsOptions
    {
        // whether to perform mesh maintenance
        int m_perform_improvement;            
        
        // whether to allow merging and separation
        int m_allow_topology_changes;     
        
        // maximum allowable change in volume when performing mesh maintenance
        double m_max_volume_change;       
        
        // edges shorter than this length will be collapsed
        double m_min_edge_length;              
        
        // edges longer then this length will be subdivided
        double m_max_edge_length;              
        
        double m_min_triangle_area;
        
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
        void *m_subdivision_scheme;   
        
        /// Whether to enforce collision-free surfaces (including during mesh maintenance operations)
        bool m_collision_safety;
        
        /// Wether to allow non-manifold (edges incident on more than two triangles)
        bool m_allow_non_manifold;
        
    };
    
    // ---------------------------------------------------------
    ///
    /// Options for position integration
    ///
    // ---------------------------------------------------------
    
    struct ElTopoIntegrationOptions
    {
        double m_friction_coefficient;
        
        // integration timestep size
        double m_dt;                     
    };
    
    
    // ---------------------------------------------------------
    ///
    /// 
    ///
    // ---------------------------------------------------------
    
    struct ElTopoMesh
    {
        
        // Number of vertices
        int num_vertices;
        
        // Vertex positions in 3D space.  This array should be shaped like: [x0 y0 z0 x1 y1 z1 ... ]
        double* vertex_locations;
        
        // Number of traingles
        int num_triangles;
        
        // Triangles, indexing into the vertex array. [a0 b0 c0 a1 b1 c1 ... ], where triangle i has vertices (ai,bi,ci).
        int* triangles;
        
        // For each vertex i, array element i is 1 if the vertex is solid, 0 otherwise
        double* vertex_masses;
        
    };
    
    
    // ---------------------------------------------------------
    ///
    /// 
    ///
    // ---------------------------------------------------------
    
    struct ElTopoDefragInformation
    {
        int num_vertex_changes;    // = N
        int* vertex_is_remove;     // boolean, size N
        int* vertex_index;         // size N
        int* split_edge;           // size 2*N
        
        int num_triangle_changes;  // = N
        int* triangle_is_remove;   // boolean, size N
        int* triangle_index;       // N
        int* new_tri;              // 3*N
        
        int defragged_triangle_map_size;    // = N
        int* defragged_triangle_map;        // 2*N
        
        int defragged_vertex_map_size;      // = N
        int* defragged_vertex_map;          // 2*N
        
    };
    
    
    // =========================================================
    //  API FUNCTIONS
    // =========================================================   
    
    
    // ---------------------------------------------------------
    ///
    /// Static operations: edge collapse, edge split, edge flip, null-space 
    /// smoothing, and topological changes
    ///
    /// Parameters:
    ///   in_num_vertices:           (Input) Number of input vertices.
    ///   in_vertex_locations:       (Input) Flattened array of vertex coordinates.
    ///                                      Elements 0, 1, and 2 are x, y, and z
    ///                                      coordinates of the first vertex, etc.
    ///   in_num_triangles           (Input) Number of input triangles.
    ///   in_triangles               (Input) Triples of vertex indices.  Elements 
    ///                                      0, 1, and 2 are indices of vertices
    ///                                      comprising the first triangle, etc.
    ///   in_masses                  (Input) Input vertex masses.
    ///   general_otions             (Input) Structure specifying options common to
    ///                                      static operations and integration.
    ///   options                    (Input) Structure specifying options specific 
    ///                                      to static operations.
    ///   out_num_vertices           (Output) Number of vertices after static 
    ///                                       operations.
    ///   out_vertex_locations       (Output, allocated by El Topo) Vertex 
    ///                                       locations after static operations.
    ///   out_num_triangles          (Output) Number of triangles after static 
    ///                                       operations.
    ///   out_triangles              (Output, allocated by El Topo) Triangles after
    ///                                       static operations.
    ///   out_masses                 (Output, allocated by El Topo) Vertex masses 
    ///                                       after static operations.
    ///
    // ---------------------------------------------------------
    
    void el_topo_static_operations( const struct ElTopoMesh* inputs,
                                   const struct ElTopoGeneralOptions* general_options,
                                   const struct ElTopoStaticOperationsOptions* options, 
                                   struct ElTopoDefragInformation* defrag_info, 
                                   struct ElTopoMesh* outputs );
    
    // ---------------------------------------------------------
    ///
    /// Free memory allocated by static operations.
    ///
    // ---------------------------------------------------------
    
    void el_topo_free_static_operations_results( struct ElTopoMesh* outputs, struct ElTopoDefragInformation* defrag_info  );
    
    
    // ---------------------------------------------------------
    ///
    /// Surface vertex position integration.
    ///
    /// Parameters:
    ///   num_vertices:              (Input) Number of input vertices (not changed 
    ///                                      by integration).
    ///   in_vertex_locations:       (Input) Flattened array of vertex coordinates.
    ///                                      Elements 0, 1, and 2 are x, y, and z
    ///                                      coordinates of the first vertex, etc.
    ///   in_new_vertex_locations:   (Input) Predicted vertex coordinates.
    ///   num_triangles              (Input) Number of input triangles (not changed
    ///                                      by integration).
    ///   triangles                  (Input) Triples of vertex indices (not changed
    ///                                      by integration).  Elements 0, 1, and 
    ///                                      2 are indices of vertices comprising 
    ///                                      the first triangle, and so on.
    ///   masses                     (Input) Input vertex masses (not changed by 
    ///                                      integration).
    ///   general_otions             (Input) Structure specifying options common to
    ///                                      static operations and integration.
    ///   options                    (Input) Structure specifying options specific 
    ///                                      to integration.
    ///   out_vertex_locations       (Output, allocated by El Topo) Vertex 
    ///                                      locations after integration.
    ///   out_dt                     (Output) Actual timestep used during integration
    ///
    // ---------------------------------------------------------
    
    void el_topo_integrate(const ElTopoMesh* inputs,
                           const double* in_vertex_new_locations,
                           const struct ElTopoGeneralOptions* general_options,
                           const struct ElTopoIntegrationOptions* options,
                           double **out_vertex_locations,
                           double *out_dt );
    
    
    // ---------------------------------------------------------
    ///
    /// Free memory allocated by integration.
    ///
    // ---------------------------------------------------------
    
    void el_topo_free_integrate_results( double* out_vertex_locations );
    
    
#ifdef __cplusplus
}
#endif

#endif
