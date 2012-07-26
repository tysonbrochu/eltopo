#ifndef NOISE_H
#define NOISE_H

#include <vec.h>

struct Noise2
{
    Noise2(unsigned int seed=171717);   
    virtual ~Noise2() {};
    
    void reinitialize(unsigned int seed);
    double operator()(double x, double y) const;
    double operator()(const Vec2d &x) const { return (*this)(x[0], x[1]); }
    
protected:
    static const unsigned int n=128;
    Vec2d basis[n];
    int perm[n];
    
    unsigned int hash_index(int i, int j) const
    { return perm[(perm[i%n]+j)%n]; }
};

struct Noise3
{
    Noise3(unsigned int seed=171717);
    virtual ~Noise3() {};
    
    void reinitialize(unsigned int seed);
    double operator()(double x, double y, double z) const;
    double operator()(const Vec3d &x) const { return (*this)(x[0], x[1], x[2]); }
    
protected:
    static const unsigned int n=128;
    Vec3d basis[n];
    int perm[n];
    
    unsigned int hash_index(int i, int j, int k) const
    { return perm[(perm[(perm[i%n]+j)%n]+k)%n]; }
};

struct Noise4
{
    Noise4(unsigned int seed=171717);
    void reinitialize(unsigned int seed);
    double operator()(double x, double y, double z, double t) const;
    double operator()(const Vec4d &x) const { return (*this)(x[0], x[1], x[2], x[3]); }
    
protected:
    static const unsigned int n=128;
    Vec4d basis[n];
    int perm[n];
    
    unsigned int hash_index(int i, int j, int k, int l) const
    { return perm[(perm[(perm[(perm[i%n]+j)%n]+k)%n]+l)%n]; }
};

// FlowNoise classes - time varying versions of some of the above

struct FlowNoise2: public Noise2
{
    FlowNoise2(unsigned int seed=171717, double spin_variation=0.2);
    void set_time(double t); // period of repetition is approximately 1
    
protected:
    Vec2d original_basis[n];
    double spin_rate[n];
};

struct FlowNoise3: public Noise3
{
    FlowNoise3(unsigned int seed=171717, double spin_variation=0.2);
    void set_time(double t); // period of repetition is approximately 1
    
protected:
    Vec3d original_basis[n];
    double spin_rate[n];
    Vec3d spin_axis[n];
};

#endif
