//
//  kdtree.h
//  Renderer
//
//  Created by Raphael Kargon on 6/13/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#ifndef __Renderer__kdtree__
#define __Renderer__kdtree__

#include "geom.h"
#include <fstream>
#include <iomanip>

#define KDTfMIN_FACES 1
#define KDTfMAX_DEPTH 20
// used for SAH calculations
#define KDT_CTRAVERSAL 1
#define KDT_CINTERSECT 1.5

/**
 *  KDTree for hierarchically grouping the faces of a mesh in order to speed up
 * ray-face intersection.
 *  Each level splits faces into two groups along a different axis. Uses
 * "Surface Area Heuristic
 *
 *  Source: "On building fast kd-Trees for Ray Tracing, and on doing that in O(N
 * log N)" (Ingo Wald, Vlastimil Havran 2006)
 *
 *  TODO: Recursion ends up to deep, subdivides faces too much
 */
class kdtree {
public:
  std::unique_ptr<kdtree> lower;
  std::unique_ptr<kdtree> upper;
  std::vector<face *> faces;
  double pos;
  int axis = -1;
  bool planarside; // whether or not faces on the split plane go to the right
                   // (upper) child node
  bounds boundingbox;
  static double completion;

  std::vector<std::pair<vertex, vertex>> wireframe();
  std::vector<std::pair<vertex, vertex>> wireframe(bool isRoot, int depth);
  // TODO put all stats in a struct
  void calc_stats(int &nodes, int &leaf_nodes, int &non_empty_leaf_nodes,
                  int &triangles_in_leaves, double &est_traversals,
                  double &est_leaves_visited, double &est_tris_intersected,
                  double &estcost, double root_area = -1);
  void print_stats();
  static std::unique_ptr<kdtree> build_tree(const std::vector<face *> &faces);
  static face *ray_tree_intersect(const kdtree *kdt, const ray &r, bool lazy,
                                  vertex *tuv, bool isRoot = true);

private:
  typedef struct plane_data {
    double plane_pos;
    int plane_axis;
    double split_cost;
    bool planar_side; // true if planar faces go to the right
  } plane_data;

  class face_wrapper {
  public:
    static int i; // keeps track of how many objects are created, for debugging
                  // purposes.
    face *f;
    // last two bits of int correspond to left, right subset during
    // classification.
    // 00 = 0 = none
    // 01 = 1 = right
    // 10 = 2 = left
    // 11 = 3 = both
    int classification;

    face_wrapper(face *f, int classification);
    face_wrapper(face *f);
    ~face_wrapper();
    static std::vector<face *>
    to_face_list(const std::vector<face_wrapper *> &wrapped_faces);
    static std::vector<face_wrapper *>
    to_wrapper_list(const std::vector<face *> &faces);
  };

  class planar_event {
  public:
    static int i; // keeps track of how many objects are created, for debugging
                  // purposes.
    face_wrapper *fw;
    double pos;

    // should probably use enums for these
    // 0 - x
    // 1 - y
    // 2 - z
    int axis;

    // 0 - f ends at p
    // 1 - f is planar on p
    // 2 - f starts at p
    int type;

    planar_event(face_wrapper *fw, double p, int k, int type);
    ~planar_event();

    class planar_event_comparator {
    public:
      bool operator()(const planar_event *pe1, const planar_event *pe2);
    };
  };

  kdtree(std::vector<planar_event *> &events, const bounds &facebounds,
         int depth);
  static std::vector<planar_event *>
  build_event_list(const std::vector<face_wrapper *> &faces);
  static plane_data
  find_optimal_plane(int nfaces, const std::vector<planar_event *> &events,
                     bounds &boundingbox, double area);
  static inline double lambda(double lower_area, double upper_area, double Nl,
                              double Nr) {
    return ((Nl == 0 || Nr == 0) && !(lower_area == 1 || upper_area == 1)) ? 0.8
                                                                           : 1;
  }
  static std::pair<double, bool> SAH(const bounds &boundingbox, double position,
                                     int axis, double area, int Nl, int Np,
                                     int Nr);
  static void
  generate_clipped_events(face_wrapper *fw, double pos, int axis,
                          const bounds &boundingbox,
                          std::vector<planar_event *> &newevents_left,
                          std::vector<planar_event *> &newevents_right);
};
#endif /* defined(__Renderer__kdtree__) */
