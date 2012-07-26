// ---------------------------------------------------------
//
//  meshoperator.h
//  Tyson Brochu 2012
//  
//  Base class for mesh processing operations.
//
// ---------------------------------------------------------

#ifndef EL_TOPO_MESHOPERATOR_H
#define EL_TOPO_MESHOPERATOR_H

#include <commonoptions.h>
#include <vector>

class SurfTrack;

class MeshOperator
{
    
protected:
    
    SurfTrack& m_surf;
    
public:
    
    MeshOperator( SurfTrack& surf ) : 
        m_surf(surf) 
    {}
    
    virtual ~MeshOperator() {}
    
    virtual void process_mesh() = 0;
        
};

#endif
