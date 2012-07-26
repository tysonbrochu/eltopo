// ---------------------------------------------------------
//
//  scriptinit.h
//  Tyson Brochu 2011
//
//  Parse a script text file to initialize the simulation and mesh obects.
//
// ---------------------------------------------------------

#ifndef EL_TOPO_SCRIPTINIT_H
#define EL_TOPO_SCRIPTINIT_H

#include <newparser.h>
#include <surftrack.h>
#include <vec.h>
#include <vector>

class MeshDriver;

class ScriptInit
{
    
public:
    
    ScriptInit() :
    output_path(),
    output_path_is_relative(false),
    frame_dt( UNINITIALIZED_DOUBLE ),
    sim_dt( UNINITIALIZED_DOUBLE ),
    end_sim_t( UNINITIALIZED_DOUBLE ),
    curr_t_specified(false),
    curr_t(0.0),
    vertices(),
    triangles(),
    masses(),
    surf_track_params(),
    driver( NULL ),
    camera_target( UNINITIALIZED_DOUBLE ),
    camera_distance( UNINITIALIZED_DOUBLE ),
    camera_heading( UNINITIALIZED_DOUBLE ),
    camera_pitch( UNINITIALIZED_DOUBLE )
    {}
    
    void parse_script( const char* filename );
    
private:
    
    ScriptInit( const ScriptInit& other );
    const ScriptInit& operator=( const ScriptInit& other );
    
    void parse_surftrack_parameters( const ParseTree& surftrack_branch );
    
    // Mesh drivers
    // --------
    
    void parse_faceoff( const ParseTree& faceoeff_sim_branch );   
    void parse_normal( const ParseTree& normal_sim_branch );   
    void parse_mean_curvature( const ParseTree& mean_curvature_sim_branch );
    void parse_sisc_curl_noise( const ParseTree& sisc_curl_noise_sim_branch );
    void parse_enright( const ParseTree& enright_sim_branch );      

    void parse_camera( const ParseTree& camera_branch );
    
    // Geometry
    // --------
    
    void parse_sheet( const ParseTree& sheet_branch );
    void parse_curved_sheet( const ParseTree& curved_sheet_branch );
    void parse_sphere( const ParseTree& sphere_branch );
    void parse_dumbbell( const ParseTree& dumbbell_branch );
    
public:
    
    
    // Simulation settings
    // --------
    
    std::string output_path;
    bool output_path_is_relative;
    
    double frame_dt;
    double sim_dt;
    double end_sim_t;
    
    bool curr_t_specified;
    double curr_t;
        
    // Surface geometry
    // --------
    
    std::vector<Vec3d> vertices;
    std::vector<Vec3st> triangles;
    std::vector<double> masses;
    

    // SurfTrack
    // --------
    
    SurfTrackInitializationParameters surf_track_params;
    
    // MeshDriver
    // --------
    
    MeshDriver* driver;
    
    // GUI
    // --------
    
    Vec3d camera_target;
    double camera_distance, camera_heading, camera_pitch;
    
};

#endif


