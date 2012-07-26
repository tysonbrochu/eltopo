#include <array2_utils.h>
#include <marching_triangles.h>

void MarchingTriangles::
contour_grid(void)
{
   edge.resize(0);
   x.resize(0);
   edge_cross.clear();
   for(int j=0; j+1<(int)phi.nj; ++j) for(int i=0; i+1<(int)phi.ni; ++i)
      contour_square(i,j);
   
   //improve_mesh();
   
}

void MarchingTriangles::
contour_square(int i, int j)
{
   contour_triangle(Vec2i(i,j), Vec2i(i+1,j), Vec2i(i+1,j+1), phi(i,j), phi(i+1,j), phi(i+1,j+1));
   contour_triangle(Vec2i(i,j), Vec2i(i+1,j+1), Vec2i(i,j+1), phi(i,j), phi(i+1,j+1), phi(i,j+1));
}

// assumes triangle is oriented counterclockwise
void MarchingTriangles::
contour_triangle(const Vec2i& x0, const Vec2i& x1, const Vec2i& x2, float p0, float p1, float p2)
{
   // guard against topological degeneracies
   if(p0==0) p0=1e-30f;
   if(p1==0) p1=1e-30f;
   if(p2==0) p2=1e-30f;

   if(p0<0){
      if(p1<0){
         if(p2<0) return; // no contour here
         else /* p2>0 */ edge.push_back(Vec2i(find_edge_cross(x1,x2,p1,p2), find_edge_cross(x0,x2,p0,p2)));
      }else{ // p1>0
         if(p2<0)        edge.push_back(Vec2i(find_edge_cross(x0,x1,p0,p1), find_edge_cross(x1,x2,p1,p2)));
         else /* p2>0 */ edge.push_back(Vec2i(find_edge_cross(x0,x1,p0,p1), find_edge_cross(x0,x2,p0,p2)));
      }
   }else{ // p0>0
      if(p1<0){
         if(p2<0)        edge.push_back(Vec2i(find_edge_cross(x0,x2,p0,p2), find_edge_cross(x0,x1,p0,p1)));
         else /* p2>0 */ edge.push_back(Vec2i(find_edge_cross(x1,x2,p1,p2), find_edge_cross(x0,x1,p0,p1)));
      }else{ // p1>0
         if(p2<0)        edge.push_back(Vec2i(find_edge_cross(x0,x2,p0,p2), find_edge_cross(x1,x2,p1,p2)));
         else /* p2>0 */ return; // no contour here
      }
   }
}

// return the vertex of the edge crossing (create it if necessary) between given grid points and function values
int MarchingTriangles::
find_edge_cross(const Vec2i& x0, const Vec2i& x1, float p0, float p1)
{
   unsigned int vertex_index;
   if(edge_cross.get_entry(Vec4i(x0.v[0], x0.v[1], x1.v[0], x1.v[1]), vertex_index)){
      return vertex_index;
   }else if(edge_cross.get_entry(Vec4i(x1.v[0], x1.v[1], x0.v[0], x0.v[1]), vertex_index)){
      return vertex_index;
   }else{
      float a=p1/(p1-p0), b=1-a;
      vertex_index=(int)x.size();
      x.push_back(Vec2f(origin[0]+dx*(a*x0[0]+b*x1[0]),
                        origin[1]+dx*(a*x0[1]+b*x1[1])));
      edge_cross.add(Vec4i(x0.v[0], x0.v[1], x1.v[0], x1.v[1]), vertex_index);
      return vertex_index;
   }
}


float MarchingTriangles::eval( float i, float j )
{
//   int p, q;
//   double f, g;
//
//   get_barycentric(i+0.5f, p, f, 1, phi.ni);
//   get_barycentric(j+0.5f, q, g, 1, phi.nj);
//   
//   double wx0, wx1, wx2, wy0, wy1, wy2;
//   quadratic_bspline_weights(f, wx0, wx1, wx2);
//   quadratic_bspline_weights(g, wy0, wy1, wy2);
//
//   return wx0*( wy0*phi(p-1,q-1) + wy1*phi(p-1,q) + wy2*phi(p-1,q+1) )
//         +wx1*( wy0*phi(p,q-1) + wy1*phi(p,q) + wy2*phi(p,q+1) )
//         +wx2*( wy0*phi(p+1,q-1) + wy1*phi(p+1,q) + wy2*phi(p+1,q+1) );

   return interpolate_value( Vec2f(i,j), phi );
}


void MarchingTriangles::eval_gradient( float y, float z, Vec2f& g )
{
   interpolate_gradient( g, Vec2f(y,z), phi );
}



void MarchingTriangles::improve_mesh()
{
   // first get adjacency information
   std::vector<Array1ui> nbr(x.size());
   for(unsigned int e=0; e<edge.size(); ++e)
   {
      unsigned int p = edge[e][0];
      unsigned int q = edge[e][1];
      nbr[p].add_unique(q);
      nbr[q].add_unique(p);
   }
   
   // then sweep through the mesh a few times incrementally improving positions
   for(unsigned int sweep=0; sweep<3; ++sweep)
   {
      for(unsigned int p=0; p<x.size(); ++p)
      {
         // get a weighted average of neighbourhood positions
         Vec2f target = x[p];
         target += x[ nbr[p][0] ];
         target += x[ nbr[p][1] ];
         target /= 3.0f;
         
         // project onto level set surface with Newton
         for(int projection_step=0; projection_step<5; ++projection_step)
         {
            float i=(target[0]-origin[0])/dx, j=(target[1]-origin[1])/dx;
            float f=eval(i,j);
            Vec2f g; 
            eval_gradient(i,j,g);
            
            float m2=mag2(g), m=std::sqrt(m2);
            float alpha = clamp( -f/(m2+1e-10f), -0.25f*m, 0.25f*m ); // clamp to avoid stepping more than a fraction of a grid cell
            
            // do line search to make sure we actually are getting closer to the zero level set
            bool line_search_success=false;
            for(int line_search_step=0; line_search_step<10; ++line_search_step)
            {
               double fnew=eval( i+alpha*g[0], j+alpha*g[1] );
               if(std::fabs(fnew) <= std::fabs(f))
               {
                  target += Vec2f( (alpha*dx)*g[0], (alpha*dx)*g[0] );
                  target += Vec2f( (alpha*dx)*g[1], (alpha*dx)*g[1] );
                  line_search_success=true;
                  break;
               }else
                  alpha*=0.5f;
            }
            if(!line_search_success)   // if we stalled trying to find the zero isocontour...
            { 
               // weight the target closer to the original x[p]
               std::cout<<"line search failed (p="<<p<<" project="<<projection_step<<" sweep="<<sweep<<")"<<std::endl;
               target= 0.5f * (x[p]+target);
            }
         }
         x[p]=target;
      }
   }
   
   
   
}






