#ifndef CURLNOISE_H
#define CURLNOISE_H

// base classes for 3d examples

#include <vec.h>

struct CurlNoise3
{
    double t;          // time
    double delta_x;    // used for finite difference approximations of curl
    
    CurlNoise3(void)
    : t(0), delta_x(1e-4)
    {}
    
    virtual ~CurlNoise3(void)
    {}
    
    virtual bool seed_particles( std::vector<ElTopo::Vec3d> &x, double dt) const;
    
    virtual ElTopo::Vec3d potential(double x, double y, double z) const;
    
    virtual void advance_time(double dt)
    {
        t+=dt;
    }
    
    void get_velocity(const ElTopo::Vec3d &x, ElTopo::Vec3d &v) const
    {
        v[0]=( (potential(x[0], x[1]+delta_x, x[2])[2] - potential(x[0], x[1]-delta_x, x[2])[2])
              -(potential(x[0], x[1], x[2]+delta_x)[1] - potential(x[0], x[1], x[2]-delta_x)[1]) ) / (2*delta_x);
        v[1]=( (potential(x[0], x[1], x[2]+delta_x)[0] - potential(x[0], x[1], x[2]-delta_x)[0])
              -(potential(x[0]+delta_x, x[1], x[2])[2] - potential(x[0]-delta_x, x[1], x[2])[2]) ) / (2*delta_x);
        v[2]=( (potential(x[0]+delta_x, x[1], x[2])[1] - potential(x[0]-delta_x, x[1], x[2])[1])
              -(potential(x[0], x[1]+delta_x, x[2])[0] - potential(x[0], x[1]-delta_x, x[2])[0]) ) / (2*delta_x);
    }
};

inline bool CurlNoise3::seed_particles( std::vector<ElTopo::Vec3d>&, double ) const
{ 
    return false; 
}

inline ElTopo::Vec3d CurlNoise3::potential( double, double, double ) const
{ 
    return ElTopo::Vec3d(0.0, 0.0, 0.0);
}

#endif
