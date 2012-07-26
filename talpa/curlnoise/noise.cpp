#include <noise.h>

namespace {
    
    template<unsigned int N> 
    Vec<N,double> sample_sphere(unsigned int &seed)
    {
        Vec<N,double> v;
        double m2;
        do{
            for(unsigned int i=0; i<N; ++i){
                v[i]=randhashf(seed++,-1,1);
            }
            m2=mag2(v);
        }while(m2>1 || m2==0);
        return v/std::sqrt(m2);
    }
    
} // namespace

//============================================================================

Noise2::
Noise2(unsigned int seed)
{
    for(unsigned int i=0; i<n; ++i){
        double theta=(double)(i*2*M_PI)/n;
        basis[i][0]=std::cos(theta);
        basis[i][1]=std::sin(theta);
        perm[i]=i;
    }
    reinitialize(seed);
}

void Noise2::
reinitialize(unsigned int seed)
{
    for(unsigned int i=1; i<n; ++i){
        int j=randhash(seed++)%(i+1);
        swap(perm[i], perm[j]);
    }
}

double Noise2::
operator()(double x, double y) const
{
    double floorx=std::floor(x), floory=std::floor(y);
    int i=(int)floorx, j=(int)floory;
    const Vec2d &n00=basis[hash_index(i,j)];
    const Vec2d &n10=basis[hash_index(i+1,j)];
    const Vec2d &n01=basis[hash_index(i,j+1)];
    const Vec2d &n11=basis[hash_index(i+1,j+1)];
    double fx=x-floorx, fy=y-floory;
    double sx=fx*fx*fx*(10-fx*(15-fx*6)),
    sy=fy*fy*fy*(10-fy*(15-fy*6));
    return bilerp(    fx*n00[0] +     fy*n00[1],
                  (fx-1)*n10[0] +     fy*n10[1],
                  fx*n01[0] + (fy-1)*n01[1],
                  (fx-1)*n11[0] + (fy-1)*n11[1],
                  sx, sy);
}

//============================================================================

Noise3::
Noise3(unsigned int seed)
{
    for(unsigned int i=0; i<n; ++i){
        basis[i]=sample_sphere<3>(seed);
        perm[i]=i;
    }
    reinitialize(seed);
}

void Noise3::
reinitialize(unsigned int seed)
{
    for(unsigned int i=1; i<n; ++i){
        int j=randhash(seed++)%(i+1);
        swap(perm[i], perm[j]);
    }
}

double Noise3::
operator()(double x, double y, double z) const
{
    double floorx=std::floor(x), floory=std::floor(y), floorz=std::floor(z);
    int i=(int)floorx, j=(int)floory, k=(int)floorz;
    const Vec3d &n000=basis[hash_index(i,j,k)];
    const Vec3d &n100=basis[hash_index(i+1,j,k)];
    const Vec3d &n010=basis[hash_index(i,j+1,k)];
    const Vec3d &n110=basis[hash_index(i+1,j+1,k)];
    const Vec3d &n001=basis[hash_index(i,j,k+1)];
    const Vec3d &n101=basis[hash_index(i+1,j,k+1)];
    const Vec3d &n011=basis[hash_index(i,j+1,k+1)];
    const Vec3d &n111=basis[hash_index(i+1,j+1,k+1)];
    double fx=x-floorx, fy=y-floory, fz=z-floorz;
    double sx=fx*fx*fx*(10-fx*(15-fx*6)),
    sy=fy*fy*fy*(10-fy*(15-fy*6)),
    sz=fz*fz*fz*(10-fz*(15-fz*6));
    return trilerp(    fx*n000[0] +     fy*n000[1] +     fz*n000[2],
                   (fx-1)*n100[0] +     fy*n100[1] +     fz*n100[2],
                   fx*n010[0] + (fy-1)*n010[1] +     fz*n010[2],
                   (fx-1)*n110[0] + (fy-1)*n110[1] +     fz*n110[2],
                   fx*n001[0] +     fy*n001[1] + (fz-1)*n001[2],
                   (fx-1)*n101[0] +     fy*n101[1] + (fz-1)*n101[2],
                   fx*n011[0] + (fy-1)*n011[1] + (fz-1)*n011[2],
                   (fx-1)*n111[0] + (fy-1)*n111[1] + (fz-1)*n111[2],
                   sx, sy, sz);
}


//============================================================================

Noise4::
Noise4(unsigned int seed)
{
    for(unsigned int i=0; i<n; ++i){
        basis[i]=sample_sphere<4>(seed);
        perm[i]=i;
    }
    reinitialize(seed);
}

void Noise4::
reinitialize(unsigned int seed)
{
    for(unsigned int i=1; i<n; ++i){
        int j=randhash(seed++)%(i+1);
        swap(perm[i], perm[j]);
    }
}

double Noise4::
operator()(double x, double y, double z, double t) const
{
    double floorx=std::floor(x), floory=std::floor(y), floorz=std::floor(z), floort=std::floor(t);
    int i=(int)floorx, j=(int)floory, k=(int)floorz, l=(int)floort;
    const Vec4d &n0000=basis[hash_index(i,j,k,l)];
    const Vec4d &n1000=basis[hash_index(i+1,j,k,l)];
    const Vec4d &n0100=basis[hash_index(i,j+1,k,l)];
    const Vec4d &n1100=basis[hash_index(i+1,j+1,k,l)];
    const Vec4d &n0010=basis[hash_index(i,j,k+1,l)];
    const Vec4d &n1010=basis[hash_index(i+1,j,k+1,l)];
    const Vec4d &n0110=basis[hash_index(i,j+1,k+1,l)];
    const Vec4d &n1110=basis[hash_index(i+1,j+1,k+1,l)];
    const Vec4d &n0001=basis[hash_index(i,j,k,l+1)];
    const Vec4d &n1001=basis[hash_index(i+1,j,k,l+1)];
    const Vec4d &n0101=basis[hash_index(i,j+1,k,l+1)];
    const Vec4d &n1101=basis[hash_index(i+1,j+1,k,l+1)];
    const Vec4d &n0011=basis[hash_index(i,j,k+1,l+1)];
    const Vec4d &n1011=basis[hash_index(i+1,j,k+1,l+1)];
    const Vec4d &n0111=basis[hash_index(i,j+1,k+1,l+1)];
    const Vec4d &n1111=basis[hash_index(i+1,j+1,k+1,l+1)];
    double fx=x-floorx, fy=y-floory, fz=z-floorz, ft=t-floort;
    double sx=fx*fx*fx*(10-fx*(15-fx*6)),
    sy=fy*fy*fy*(10-fy*(15-fy*6)),
    sz=fz*fz*fz*(10-fz*(15-fz*6)),
    st=ft*ft*ft*(10-ft*(15-ft*6));
    return quadlerp(    fx*n0000[0] +     fy*n0000[1] +     fz*n0000[2] +     ft*n0000[3],
                    (fx-1)*n1000[0] +     fy*n1000[1] +     fz*n1000[2] +     ft*n1000[3],
                    fx*n0100[0] + (fy-1)*n0100[1] +     fz*n0100[2] +     ft*n0100[3],
                    (fx-1)*n1100[0] + (fy-1)*n1100[1] +     fz*n1100[2] +     ft*n1100[3],
                    fx*n0010[0] +     fy*n0010[1] + (fz-1)*n0010[2] +     ft*n0010[3],
                    (fx-1)*n1010[0] +     fy*n1010[1] + (fz-1)*n1010[2] +     ft*n1010[3],
                    fx*n0110[0] + (fy-1)*n0110[1] + (fz-1)*n0110[2] +     ft*n0110[3],
                    (fx-1)*n1110[0] + (fy-1)*n1110[1] + (fz-1)*n1110[2] +     ft*n1110[3],
                    fx*n0001[0] +     fy*n0001[1] +     fz*n0001[2] + (1-ft)*n0001[3],
                    (fx-1)*n1001[0] +     fy*n1001[1] +     fz*n1001[2] + (1-ft)*n1001[3],
                    fx*n0101[0] + (fy-1)*n0101[1] +     fz*n0101[2] + (1-ft)*n0101[3],
                    (fx-1)*n1101[0] + (fy-1)*n1101[1] +     fz*n1101[2] + (1-ft)*n1101[3],
                    fx*n0011[0] +     fy*n0011[1] + (fz-1)*n0011[2] + (1-ft)*n0011[3],
                    (fx-1)*n1011[0] +     fy*n1011[1] + (fz-1)*n1011[2] + (1-ft)*n1011[3],
                    fx*n0111[0] + (fy-1)*n0111[1] + (fz-1)*n0111[2] + (1-ft)*n0111[3],
                    (fx-1)*n1111[0] + (fy-1)*n1111[1] + (fz-1)*n1111[2] + (1-ft)*n1111[3],
                    sx, sy, sz, st);
}

//============================================================================

FlowNoise2::
FlowNoise2(unsigned int seed, double spin_variation)
: Noise2(seed)
{
    seed+=n;
    for(unsigned int i=0; i<n; ++i){
        original_basis[i]=basis[i];
        spin_rate[i]=2.0*(double)M_PI*randhashd(seed++, 1.0 - 0.5*spin_variation, 1.0 + 0.5*spin_variation);
    }
}

void FlowNoise2::
set_time(double t)
{
    for(unsigned int i=0; i<n; ++i){
        double theta=spin_rate[i]*t;
        double c=std::cos(theta), s=std::sin(theta);
        basis[i][0]= c*original_basis[i][0]+s*original_basis[i][1];
        basis[i][1]=-s*original_basis[i][0]+c*original_basis[i][1];
    }
}

//============================================================================

FlowNoise3::
FlowNoise3(unsigned int seed, double spin_variation)
: Noise3(seed)
{
    seed+=8*n; // probably avoids overlap with sequence used in initializing superclass Noise3
    for(unsigned int i=0; i<n; ++i){
        original_basis[i]=basis[i];
        spin_axis[i]=sample_sphere<3>(seed);
        spin_rate[i]=2.0*M_PI*randhashd(seed++, 0.1 - 0.5*spin_variation, 0.1 + 0.5*spin_variation);
    }
}

void FlowNoise3::
set_time(double t)
{
    for(unsigned int i=0; i<n; ++i){
        double theta=spin_rate[i]*t;
        double c=std::cos(theta), s=std::sin(theta);
        // form rotation matrix
        double R00=c+(1-c)*sqr(spin_axis[i][0]),
        R01=(1-c)*spin_axis[i][0]*spin_axis[i][1]-s*spin_axis[i][2],
        R02=(1-c)*spin_axis[i][0]*spin_axis[i][2]+s*spin_axis[i][1];
        double R10=(1-c)*spin_axis[i][0]*spin_axis[i][1]+s*spin_axis[i][2],
        R11=c+(1-c)*sqr(spin_axis[i][1]),
        R12=(1-c)*spin_axis[i][1]*spin_axis[i][2]-s*spin_axis[i][0];
        double R20=(1-c)*spin_axis[i][0]*spin_axis[i][2]-s*spin_axis[i][1],
        R21=(1-c)*spin_axis[i][1]*spin_axis[i][2]+s*spin_axis[i][0],
        R22=c+(1-c)*sqr(spin_axis[i][2]);
        basis[i][0]=R00*original_basis[i][0] + R01*original_basis[i][1] + R02*original_basis[i][2];
        basis[i][1]=R10*original_basis[i][0] + R11*original_basis[i][1] + R12*original_basis[i][2];
        basis[i][2]=R20*original_basis[i][0] + R21*original_basis[i][1] + R22*original_basis[i][2];
    }
}

