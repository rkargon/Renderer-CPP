//
//  camera.h
//  Renderer
//
//  Represents a camera in 3D space, through which one views and renders a 3D
//  scene.
//
//  Created by Raphael Kargon on 6/8/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#ifndef __Renderer__camera__
#define __Renderer__camera__

#include "geom.h"

class camera {
public:
  vertex center; // location of the camera
  vertex focus; // What the camera is focused on, can be used to center the view
                // on an object
  vertex normal; // The normal vector pointing straight out of the camera
  vertex vert;   // The vector pointing 'up' from the camera, used along with
                 // 'normal' to determine orientation
  double fov;    // The field of view of the camera
  double mindist, maxdist; // Clipping distances
  bool ortho; // Whether the projection is orthographic or perspective
  double dof_focus_distance; // Distance at which objects are in focus for the
                             // camera
  double aperture_size;      // Width of aperture for DOF effects. 0 means no
                             // blurring due to depth of field.

  /* camera(const vertex &center = {0, -5, 0}, const vertex &focus = {}, */
  /*        const vertex &normal = {0, 1, 0}, const vertex &vert = {0, 0, 1}, */
  /*        const double fov = 0.75, const double mindist = 0.01, */
  /*        const double maxdist = 100, const bool ortho = false, */
  /*        double dof_focus_distance = 0, double aperture_size = 0); */

  camera();

  /* projections */

  ray cast_ray(const double px, const double py, const double w,
               const double h) const;
  point_2d<double> project_vertex(const vertex &v, double w, double h) const;
  double vertex_depth(const vertex &v) const;
  vertex
  view_vector(const vertex &v) const; // get vector from camera to given vertex
  double face_depth(const face &f) const;

  /* manipulation of camera */

  // Moves the camera so that is is facing the focus point
  void center_focus();
  void shift_focus(const double dx, const double dy);
  void zoom(const double zoomfactor);
  vertex image_plane_vector(const double dx, const double dy);

  /* rotation */
  void rotate_axis(const vertex &axis,
                   const double dtheta);    // rotate along arbitrary axis
  void rotate_local_x(const double dtheta); // rotate along local X axis
  void rotate_local_y(const double dtheta); // rotate along local y axis
  void rotate_local_z(const double dtheta); // rotate along local z axis
  void set_global_rotation(const double theta, const double rho,
                           const double psi); // set global rotation angles

  // recalculates image plane vectors
  void calc_image_vectors();

private:
  // keeps track of image plane vectors
  vertex _cx, _cy;
};

#endif /* defined(__Renderer__camera__) */
