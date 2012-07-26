#ifndef COLLISIONQUERIES_H
#define COLLISIONQUERIES_H

#include <vec.h>

// 2D ====================================================================================================

/// Compute distance between a 2D point and the closest point on a line segment.
///
void check_point_edge_proximity(bool update, const Vec2d &x0, const Vec2d &x1, const Vec2d &x2,
                                double &distance);

/// Compute distance between a 2D point and the closest point on a line segment.  Also return the barycentric coordinate along the
/// segment and normal.  Normal is from points 1-2 towards 0, unless normal_multiplier < 0.
///
void check_point_edge_proximity(bool update, const Vec2d &x0, const Vec2d &x1, const Vec2d &x2,
                                double &distance, double &s, Vec2d &normal, double normal_multiplier);

// 3D ====================================================================================================


/// Compute distance between a 3D point and the closest point on a line segment.
///
void check_point_edge_proximity(bool update, const Vec3d &x0, const Vec3d &x1, const Vec3d &x2,
                                double &distance);

/// Compute distance between a 3D point and the closest point on a line segment.  Also return the barycentric coordinate along the
/// segment and normal.  Normal is from points 1-2 towards 0, unless normal_multiplier < 0.
///
void check_point_edge_proximity(bool update, const Vec3d &x0, const Vec3d &x1, const Vec3d &x2,
                                double &distance, double &s, Vec3d &normal, double normal_multiplier);

/// Compute distance between the closest points between two line segments in 3D.
///
void check_edge_edge_proximity(const Vec3d &x0, const Vec3d &x1, const Vec3d &x2, const Vec3d &x3,
                               double &distance);

/// Compute distance between the closest points between two line segments in 3D, with barycentric coordinates for closest points, 
/// and a normal that points from 0-1 towards 2-3 (unreliable if distance==0 or very small).
///
void check_edge_edge_proximity(const Vec3d &x0, const Vec3d &x1, const Vec3d &x2, const Vec3d &x3,
                               double &distance, double &s0, double &s2, Vec3d &normal);

/// Compute distance between a point (x0) and the closest point on a triangle (x1-x2-x3) in 3D.
///
void check_point_triangle_proximity(const Vec3d &x0, const Vec3d &x1, const Vec3d &x2, const Vec3d &x3,
                                    double &distance);

/// Find distance between 0 and 1-2-3, with barycentric coordinates for closest point, and a normal that points from 1-2-3 towards 
/// 0 (unreliable if distance==0 or very small).
///
void check_point_triangle_proximity(const Vec3d &x0, const Vec3d &x1, const Vec3d &x2, const Vec3d &x3,
                                    double &distance, double &s1, double &s2, double &s3, Vec3d &normal);

/// Compute the signed volume of a tetrahedron.
///
double signed_volume(const Vec3d &x0, const Vec3d &x1, const Vec3d &x2, const Vec3d &x3);


#endif
