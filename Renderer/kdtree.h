//
//  kdtree.h
//  Renderer
//
//  Created by Raphael Kargon on 6/13/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#ifndef __Renderer__kdtree__
#define __Renderer__kdtree__

#include <fstream>
#include <iomanip>
#include "geom.h"

#define KDTfMIN_FACES 1
#define KDTfMAX_DEPTH 20
//used for SAH calculations
#define KDT_CTRAVERSAL 1
#define KDT_CINTERSECT 1.5

class kdtree {
public:
    kdtree *lower = nullptr;
    kdtree *upper = nullptr;
    std::vector<face *> faces;
    double pos;
    int axis = -1;
    bool planarside; //whether or not faces on the split plane go to the right (upper) child node
    bounds boundingbox;
    static double completion;
    
    std::vector<edge*> wireframe();
    std::vector<edge*> wireframe(bool isRoot, int depth);
    void calcstats(int& nodes, int& leaf_nodes, int& non_empty_leaf_nodes, int& triangles_in_leaves, double& est_traversals, double& est_leaves_visited, double& est_tris_intersected, double& estcost, double rootArea=-1);
    void printstats();
    static kdtree *buildTree(const std::vector<face*>& faces);
    static face *rayTreeIntersect(kdtree *kdt, const ray& r, bool lazy, vertex *tuv, bool isRoot = true);
    
private:

    typedef struct planedata{
        double plane_pos;
        int plane_axis;
        double splitcost;
        bool planar_side; //true if planar faces go to the right
    } planedata;
    
    class faceWrapper{
    public:
        static int i; //keeps track of how many objects are created, for debugging purposes.
        face *f;
        //last two bits of int correspond to left, right subset during classification.
		//00 = 0 = none
		//01 = 1 = right
		//10 = 2 = left
		//11 = 3 = both
        int classification;
        
        faceWrapper(face *f, int classification);
        faceWrapper(face *f);
        ~faceWrapper();
        static std::vector<face*> toFaceList(const std::vector<faceWrapper*>& wrappedfaces);
        static std::vector<faceWrapper*> toWrapperList(const std::vector<face*>& faces);
    };
    
    class planarEvent{
    public:
        static int i; //keeps track of how many objects are created, for debugging purposes.
        faceWrapper * fw;
        double pos;
        
        //should probably use enums for these
        // 0 - x
		// 1 - y
		// 2 - z
        int axis;

		// 0 - f ends at p
		// 1 - f is planar on p
		// 2 - f starts at p
        int type;
        
        planarEvent(faceWrapper *fw, double p, int k, int type);
        ~planarEvent();
        
        class planarEventComparator{
        public:
            bool operator()(const planarEvent* pe1, const planarEvent* pe2);
        };
    };
    
    kdtree(std::vector<planarEvent*>& events, const bounds& facebounds, int depth);
    static std::vector<planarEvent*> buildEventList(const std::vector<faceWrapper*>& faces);
    static planedata findOptimalPlane(int nfaces, const std::vector<planarEvent*>& events, bounds& boundingbox, double area);
    static inline double lambda(double lowerarea, double upperarea, double Nl, double Nr){
        return ((Nl==0||Nr==0) && !(lowerarea==1 || upperarea==1)) ? 0.8 : 1;
    }
    static std::pair<double, bool> SAH(const bounds& boundingbox, double position, int axis, double area, int Nl, int Np, int Nr);
    static void generateClippedEvents(faceWrapper *fw, double pos, int axis, const bounds& boundingbox, std::vector<planarEvent*>& newevents_left, std::vector<planarEvent*>& newevents_right);
};
#endif /* defined(__Renderer__kdtree__) */
