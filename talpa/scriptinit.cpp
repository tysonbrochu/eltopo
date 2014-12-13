// ---------------------------------------------------------
//
//  ScriptInit.cpp
//  Tyson Brochu 2011
//
//  Parse a script text file to initialize the simulation and mesh obects.
//
// ---------------------------------------------------------


#include <scriptinit.h>

#include <enrightdriver.h>
#include <faceoff.h>
#include <fstream>
#include <geometryinit.h>
#include <iomesh.h>
#include <meancurvature.h>
#include <normaldriver.h>
#include <sisccurlnoisedriver.h>
#include <subdivisionscheme.h>

#ifdef _MSC_VER
//this is kind of hacky, but seems to do the trick for now (on Windows, the name is prefaced by an underscore)
//so we'll just rename is here
#define snprintf _snprintf 
#endif
// ---------------------------------------------------------

void ScriptInit::parse_surftrack_parameters( const ParseTree& surftrack_branch )
{
    
    int use_fraction, perform_improvement, topology_changes, collision_safety, non_manifold;
    std::string subdivision_scheme;
    
    surftrack_branch.get_int( "use_fraction", use_fraction );
    surftrack_branch.get_number( "min_edge_length", surf_track_params.m_min_edge_length );
    surftrack_branch.get_number( "max_edge_length", surf_track_params.m_max_edge_length  );
    surftrack_branch.get_number( "max_volume_change", surf_track_params.m_max_volume_change );
    surftrack_branch.get_number( "min_triangle_angle", surf_track_params.m_min_triangle_angle );   
    surftrack_branch.get_number( "max_triangle_angle", surf_track_params.m_max_triangle_angle );      
    surftrack_branch.get_number( "min_triangle_area", surf_track_params.m_min_triangle_area );
    
    int use_curvature_when_splitting;
    if ( surftrack_branch.get_int( "use_curvature_when_splitting", use_curvature_when_splitting ) )
    {
        surf_track_params.m_use_curvature_when_splitting = ( use_curvature_when_splitting != 0 );
    }
    
    int use_curvature_when_collapsing;
    if ( surftrack_branch.get_int( "use_curvature_when_collapsing", use_curvature_when_collapsing ) )
    {
        surf_track_params.m_use_curvature_when_collapsing = ( use_curvature_when_collapsing != 0 );
    }
    
    surftrack_branch.get_number( "min_curvature_multiplier", surf_track_params.m_min_curvature_multiplier );
    surftrack_branch.get_number( "max_curvature_multiplier", surf_track_params.m_max_curvature_multiplier );
    surftrack_branch.get_number( "merge_proximity", surf_track_params.m_merge_proximity_epsilon );
    surftrack_branch.get_number( "repulsion_proximity", surf_track_params.m_proximity_epsilon );
    surftrack_branch.get_number( "friction_coefficient", surf_track_params.m_friction_coefficient );
    surftrack_branch.get_int( "perform_improvement", perform_improvement );
    surftrack_branch.get_int( "allow_topology_changes", topology_changes );   
    surftrack_branch.get_int( "collision_safety", collision_safety );   
    surftrack_branch.get_string( "subdivision_scheme", subdivision_scheme );
    
    surf_track_params.m_use_fraction = ( use_fraction != 0 );   
    surf_track_params.m_perform_improvement = (perform_improvement != 0);
    surf_track_params.m_allow_topology_changes = (topology_changes != 0);
    surf_track_params.m_collision_safety = (collision_safety != 0);   
    
    if ( surftrack_branch.get_int( "allow_non_manifold", non_manifold ) )
    {
      surf_track_params.m_allow_non_manifold = (non_manifold != 0);
    }
    
    if ( strcmp( subdivision_scheme.c_str(), "butterfly" ) == 0 )
    {
        surf_track_params.m_subdivision_scheme = new ButterflyScheme();
    }
    else
    {
        surf_track_params.m_subdivision_scheme = new MidpointScheme();
    }
    
    int allow_vertex_movement;
    if ( surftrack_branch.get_int( "allow_vertex_movement", allow_vertex_movement ) )
    {
        surf_track_params.m_allow_vertex_movement = ( allow_vertex_movement != 0 );
    }
    
}

// ---------------------------------------------------------

void ScriptInit::parse_faceoff( const ParseTree& faceoeff_sim_branch )
{
    double speed;
    faceoeff_sim_branch.get_number( "speed", speed );
    driver = new FaceOffDriver( speed, Vec3d( -0.25, 0.0, 0.0 ), Vec3d( 0.25, 0.0, 0.0 ), 0.4, 0.2 );
}

// ---------------------------------------------------------

void ScriptInit::parse_normal( const ParseTree& normal_sim_branch )
{
    double speed;
    normal_sim_branch.get_number( "speed", speed );
    driver = new NormalDriver( speed, Vec3d( -0.25, 0.0, 0.0 ), Vec3d( 0.25, 0.0, 0.0 ), 0.4, 0.2 );
}

// ---------------------------------------------------------

void ScriptInit::parse_mean_curvature( const ParseTree& mean_curvature_sim_branch )
{
    double speed;
    mean_curvature_sim_branch.get_number( "speed", speed );
    
    std::string ground_truth_file;
    mean_curvature_sim_branch.get_string( "ground_truth_file", ground_truth_file );
    
    Array3d sethian_final;
    read_signed_distance( ground_truth_file.c_str(), sethian_final );       
    
    Vec3d phi_domain_low;   
    mean_curvature_sim_branch.get_vec3d( "phi_domain_low", phi_domain_low );
    
    double phi_domain_dx;
    mean_curvature_sim_branch.get_number( "phi_domain_dx", phi_domain_dx );
    
    driver = new MeanCurvatureDriver( speed, sethian_final, phi_domain_low, phi_domain_dx );
}

// ---------------------------------------------------------

void ScriptInit::parse_sisc_curl_noise( const ParseTree&  )
{
    driver = new SISCCurlNoiseDriver( );
}

// ---------------------------------------------------------
void ScriptInit::parse_enright( const ParseTree&  )
{
    driver = new EnrightDriver( );
}

// ---------------------------------------------------------

void ScriptInit::parse_camera( const ParseTree& camera_branch )
{
    camera_branch.get_vec3d( "target", camera_target );
    camera_branch.get_number( "distance", camera_distance );
    camera_branch.get_number( "heading", camera_heading );   
    camera_branch.get_number( "pitch", camera_pitch );
}

// ---------------------------------------------------------

void ScriptInit::parse_sheet( const ParseTree& sheet_branch )
{
    Vec3d lower_corner;
    sheet_branch.get_vec3d( "corner", lower_corner );
    
    double sheet_dx;
    sheet_branch.get_number( "dx", sheet_dx );
    
    int sheet_ni, sheet_nj;
    sheet_branch.get_int( "ni", sheet_ni );
    sheet_branch.get_int( "nj", sheet_nj );
    
    std::vector<Vec3d> sheet_verts;
    std::vector<Vec3st> sheet_tris;
    
    create_sheet( sheet_verts, sheet_tris, lower_corner, sheet_dx, sheet_ni, sheet_nj );      
    
    
    Vec3d rotate_axis;
    double rotate_radians;
    if ( sheet_branch.get_vec3d( "rotate_axis", rotate_axis ) )
    {
        if ( sheet_branch.get_number( "rotate_radians", rotate_radians ) )
        {
            Vec3d centre(0,0,0);
            for ( unsigned int i = 0; i < sheet_verts.size(); ++i )
            {
                centre += sheet_verts[i];
            }
            centre /= static_cast<double>( sheet_verts.size() );
            
            for ( unsigned int i = 0; i < sheet_verts.size(); ++i )
            {
                sheet_verts[i] -= centre;
            }
            
            for ( unsigned int i = 0; i < sheet_verts.size(); ++i )
            {
                sheet_verts[i] = rotate( sheet_verts[i], rotate_radians, rotate_axis );
                sheet_verts[i] += centre;
            }
            
        }
    }
    
    
    int is_solid = 0;
    sheet_branch.get_int( "is_solid", is_solid );  
    
    std::vector<double> sheet_masses;
    
    if ( is_solid )
    {
        sheet_masses.resize( sheet_verts.size(), std::numeric_limits<double>::infinity() );      
    }
    else
    {
        sheet_masses.resize( sheet_verts.size(), 1.0 );
    }
    
    append_mesh( triangles, vertices, masses, sheet_tris, sheet_verts, sheet_masses );
    
}


// ---------------------------------------------------------

void ScriptInit::parse_curved_sheet( const ParseTree& sheet_branch )
{
    Vec3d lower_corner;
    sheet_branch.get_vec3d( "corner", lower_corner );
    
    double sheet_dx;
    sheet_branch.get_number( "dx", sheet_dx );
    
    int sheet_ni, sheet_nj;
    sheet_branch.get_int( "ni", sheet_ni );
    sheet_branch.get_int( "nj", sheet_nj );
    
    std::vector<Vec3d> sheet_verts;
    std::vector<Vec3st> sheet_tris;
    
    create_curved_sheet( sheet_verts, sheet_tris, lower_corner, sheet_dx, sheet_ni, sheet_nj );      
    
    Vec3d rotate_axis;
    double rotate_radians;
    if ( sheet_branch.get_vec3d( "rotate_axis", rotate_axis ) )
    {
        if ( sheet_branch.get_number( "rotate_radians", rotate_radians ) )
        {
            Vec3d centre(0,0,0);
            for ( unsigned int i = 0; i < sheet_verts.size(); ++i )
            {
                centre += sheet_verts[i];
            }
            centre /= static_cast<double>( sheet_verts.size() );
            
            for ( unsigned int i = 0; i < sheet_verts.size(); ++i )
            {
                sheet_verts[i] -= centre;
            }
            
            for ( unsigned int i = 0; i < sheet_verts.size(); ++i )
            {
                sheet_verts[i] = rotate( sheet_verts[i], rotate_radians, rotate_axis );
                sheet_verts[i] += centre;
            }
            
        }
    }
    
    
    std::vector<double> sheet_masses( sheet_verts.size(), 1.0 );
    
    append_mesh( triangles, vertices, masses, sheet_tris, sheet_verts, sheet_masses );
    
}

// ---------------------------------------------------------

void ScriptInit::parse_sphere( const ParseTree& sphere_branch )
{
    Vec3d sphere_center;
    sphere_branch.get_vec3d( "sphere_center", sphere_center );
    double sphere_radius;
    sphere_branch.get_number( "sphere_radius", sphere_radius );
    
    double dx;
    sphere_branch.get_number( "sphere_dx", dx );
    
    int is_solid = 0;
    sphere_branch.get_int( "is_solid", is_solid );
    
    std::vector<Vec3d> sphere_vertices;
    std::vector<Vec3st> sphere_triangles;
    
    create_sphere( sphere_center, sphere_radius, dx, sphere_vertices, sphere_triangles );
    
    std::vector<Vec3d> sphere_velocities( sphere_vertices.size(), Vec3d(0) );
    
    std::vector<double> sphere_masses;
    if ( is_solid == 0 )
    {
        sphere_masses.resize( sphere_vertices.size(), 1.0 );
    }
    else
    {
        sphere_masses.resize( sphere_vertices.size(), std::numeric_limits<double>::infinity() );         
    }
    
    append_mesh( triangles, vertices, masses, sphere_triangles, sphere_vertices, sphere_masses );
}

// ---------------------------------------------------------


void ScriptInit::parse_dumbbell( const ParseTree& dumbbell_branch )
{
    double domain_dx;
    dumbbell_branch.get_number( "domain_dx", domain_dx );
    
    Vec3d centre_a;
    dumbbell_branch.get_vec3d( "centre_a", centre_a );
    
    Vec3d centre_b;
    dumbbell_branch.get_vec3d( "centre_b", centre_b );
    
    double sphere_radius;
    dumbbell_branch.get_number( "sphere_radius", sphere_radius );
    
    double handle_width;
    dumbbell_branch.get_number( "handle_width", handle_width );   
    
    Vec3d domain_low = min_union( centre_a, centre_b ) - 2.0 * Vec3d( sphere_radius );
    Vec3d domain_high = max_union( centre_a, centre_b ) + 2.0 * Vec3d( sphere_radius );
    
    Array3d phi;   
    create_dumbbell_signed_distance( centre_a, centre_b, sphere_radius, handle_width, domain_dx, domain_low, domain_high, phi );      
    
    std::vector<Vec3st> new_tris;
    std::vector<Vec3d> new_verts;
    contour_phi( domain_low, domain_dx, phi, new_tris, new_verts );
    project_to_exact_dumbbell( new_verts, centre_a, centre_b, sphere_radius, handle_width );
    
    std::vector<double> new_masses( new_verts.size(), 1.0 );
    std::vector<Vec3d> new_velocities( new_verts.size(), Vec3d(0,0,0) );
    append_mesh( triangles, vertices, masses, new_tris, new_verts, new_masses );
}


// ---------------------------------------------------------

void ScriptInit::parse_script( const char* filename )
{
    
    std::ifstream filestream( filename );
    if ( !filestream.good() )
    {
        std::cerr << "Could not open script file" << std::endl;
        exit(1);
    }
    
    std::cout << "script file: " << filename << std::endl;
    
    ParseTree tree;
    parse_stream( filestream, tree );
    
    
    //
    // Frame stepper
    //
    
    bool ok = tree.get_number( "frame_dt", frame_dt );
    assert( ok );
    
    int num_substeps;
    bool substeps_specified = tree.get_int( "num_substeps", num_substeps );
    if ( substeps_specified )
    {
        sim_dt = frame_dt / (double) num_substeps;
    }
    
    double read_sim_dt;
    if ( tree.get_number( "sim_dt", read_sim_dt ) )
    {
        if ( substeps_specified )
        {
            std::cerr << "WARNING: Both sim_dt and num_substeps specified in config script.  Going with sim_dt." << std::endl;
        }
        
        sim_dt = read_sim_dt;
    }
    
    tree.get_number( "end_sim_t", end_sim_t );
    
    curr_t_specified = tree.get_number( "curr_t", curr_t );
    
    
    //
    // File output
    //
    
    if ( tree.get_string( "output_path", output_path ) )
    {
        output_path_is_relative = false;
    }
    else if ( tree.get_string( "relative_output_path", output_path ) )
    {
        output_path_is_relative = true;
    }         
    else
    {
        // no path specified
        output_path_is_relative = false;
        output_path = std::string( "./" );
    }
    
    //
    // OpenGL camera
    //
    
    const ParseTree* camera_branch = tree.get_branch( "camera" );
    parse_camera( *camera_branch );
    
    
    //
    // Surface geometry
    //
    
    {
        unsigned int curved_sheet_n = 0;
        char curved_sheet_name[256];
        snprintf( curved_sheet_name, 256, "curved_sheet%d", curved_sheet_n );
        const ParseTree* curved_sheet_branch = tree.get_branch( curved_sheet_name );
        
        while ( curved_sheet_branch != NULL )
        {
            parse_curved_sheet( *curved_sheet_branch );
            curved_sheet_branch = NULL;
            ++curved_sheet_n;
            snprintf( curved_sheet_name, 256, "curved_sheet%d", curved_sheet_n );
            curved_sheet_branch = tree.get_branch( curved_sheet_name );
        }
    }
    
    {
        unsigned int sheet_n = 0;   
        char sheet_name[256];
        snprintf( sheet_name, 256, "sheet%d", sheet_n );
        const ParseTree* sheet_branch = tree.get_branch( sheet_name );
        
        while ( sheet_branch != NULL )
        {
            parse_sheet( *sheet_branch );      
            sheet_branch = NULL;      
            ++sheet_n;
            snprintf( sheet_name, 256, "sheet%d", sheet_n );
            sheet_branch = tree.get_branch( sheet_name );
        }
    }
    
    {
        unsigned int sphere_n = 0;   
        char sphere_name[256];
        snprintf( sphere_name, 256, "sphere%d", sphere_n );
        const ParseTree* sphere_branch = tree.get_branch( sphere_name );
        
        while ( sphere_branch != NULL )
        {
            parse_sphere( *sphere_branch );
            sphere_branch = NULL;
            ++sphere_n;
            snprintf( sphere_name, 256, "sphere%d", sphere_n );
            sphere_branch = tree.get_branch( sphere_name );
        }
    }
    
    
    const ParseTree* sphere_branch = tree.get_branch( "sphere" );
    if ( sphere_branch != NULL )
    {
        parse_sphere( *sphere_branch );
    }      
    
    const ParseTree* dumbbell_branch = tree.get_branch( "dumbbell" );
    if ( dumbbell_branch != NULL )
    {
        parse_dumbbell( *dumbbell_branch );
    }      
    
    const ParseTree* trimesh_branch = tree.get_branch( "trimesh" );
    if ( trimesh_branch != NULL )
    {
        printf("Found trimesh branch\n");
        
        std::string meshpath;
        trimesh_branch->get_string("filepath", meshpath);
        printf("Got path: %s\n", meshpath.c_str());
        
        NonDestructiveTriMesh trimesh;
        
        printf("Reading file\n");
        
        std::vector<Vec3d> input_vertices;
        std::vector<double> in_masses;
        read_binary_file( trimesh, input_vertices, in_masses, curr_t, meshpath.c_str() );
        curr_t_specified = true;
        
        int is_solid = 0;
        if ( trimesh_branch->get_int( "is_solid", is_solid ) )
        {
            in_masses.clear();
            if ( is_solid )
            {
                in_masses.resize( input_vertices.size(), std::numeric_limits<double>::infinity() );
            }
            else
            {
                in_masses.resize( input_vertices.size(), 1.0 );
            }
        }
        
        Vec3d translate;
        if ( trimesh_branch->get_vec3d("translate", translate) )
        {
            for ( size_t i = 0; i < input_vertices.size(); ++i )
            {
                input_vertices[i] += translate;
            }
        }
                
        append_mesh( triangles, vertices, masses, trimesh.get_triangles(), input_vertices, in_masses );        
        
        printf("loaded file %s", meshpath.c_str());
    }
    
    
    const ParseTree* obj_branch = tree.get_branch( "objfile" );
    if ( obj_branch != NULL )
    {
        printf("Found obj branch\n");
        
        std::string meshpath;
        obj_branch->get_string("filepath", meshpath);
        printf("Got path: %s\n", meshpath.c_str());
        
        NonDestructiveTriMesh trimesh;
        
        printf("Reading file\n");
        
        std::vector<Vec3d> obj_vertices;
        read_objfile( trimesh, obj_vertices, meshpath.c_str() );
        
        std::vector<Vec3st> obj_triangles = trimesh.get_triangles();
        
        Vec3d translate;
        if ( obj_branch->get_vec3d("translate", translate) )
        {
            for ( size_t i = 0; i < obj_vertices.size(); ++i )
            {
                obj_vertices[i] += translate;
            }
        }
        
        std::vector<double> obj_masses(0);
        int is_solid = 0;
        obj_branch->get_int( "is_solid", is_solid );
        
        if ( is_solid )
        {
            obj_masses.resize( obj_vertices.size(), std::numeric_limits<double>::infinity() );
        }
        else
        {
            obj_masses.resize( obj_vertices.size(), 1.0 );
        }
                
        append_mesh( triangles, vertices, masses, obj_triangles, obj_vertices, obj_masses );
    }
    
    //
    // SurfTrack parameters
    //
    
    const ParseTree* surftrack_branch = tree.get_branch( "surftrack_parameters" );
    parse_surftrack_parameters( *surftrack_branch );
    
    
    //
    // Mesh drivers
    //
    
    const ParseTree* faceoeff_sim_branch = tree.get_branch( "faceoff_simulation" );
    if ( faceoeff_sim_branch != NULL )
    {
        parse_faceoff( *faceoeff_sim_branch );
    }
    
    const ParseTree* normal_sim_branch = tree.get_branch( "normal_simulation" );
    if ( normal_sim_branch != NULL )
    {
        parse_normal( *normal_sim_branch );
    }
    
    const ParseTree* mean_curvature_sim_branch = tree.get_branch( "mean_curvature_simulation" );
    if ( mean_curvature_sim_branch != NULL )
    {
        parse_mean_curvature( *mean_curvature_sim_branch );
    }
    
    const ParseTree* sisc_curl_noise_sim_branch = tree.get_branch( "sisc_curl_noise_simulation" );
    if ( sisc_curl_noise_sim_branch != NULL )
    {
        parse_sisc_curl_noise( *sisc_curl_noise_sim_branch );
    }
    
    const ParseTree* enright_branch = tree.get_branch( "enright_simulation" );
    if ( enright_branch != NULL )
    {
        parse_enright( *enright_branch );
    }
        
}

