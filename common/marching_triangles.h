#ifndef MARCHING_TRIANGLES_H
#define MARCHING_TRIANGLES_H

#include <array2.h>
#include <hashtable.h>
#include <vec.h>

struct MarchingTriangles
{
   std::vector<Vec2i> edge;
   std::vector<Vec2f> x;
   Vec2f origin;
   float dx;
   Array2f phi;

   HashTable<Vec4i,unsigned int> edge_cross; // stores vertices that have been created already at given edge crossings

   explicit MarchingTriangles(float dx_=1)
      : edge(), x(), origin(0), dx(dx_), phi(), edge_cross()

   {}

   explicit MarchingTriangles(const Vec2f &origin_, float dx_=1)
      : edge(), x(), origin(origin_), dx(dx_), phi(), edge_cross()
   {}

   void contour_grid(void);

   private:
   
   float eval( float x, float y );
   void eval_gradient( float x, float y, Vec2f& g );
   
   void contour_square(int i, int j); // add triangles for contour in the given grid square (i,j) to (i+1,j+1)
   void contour_triangle(const Vec2i& x0, const Vec2i& x1, const Vec2i& x2, float p0, float p1, float p2);
   int find_edge_cross(const Vec2i& x0, const Vec2i& x1, float p0, float p1);
   void improve_mesh();
};

#endif
