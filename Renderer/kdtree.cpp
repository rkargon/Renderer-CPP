//
//  kdtree.cpp
//  Renderer
//
//  Created by Raphael Kargon on 6/13/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#include "kdtree.h"

real kdtree::completion=0;

kdtree *kdtree::buildTree(const std::vector<face*>& faces){
    if(faces.empty()) return nullptr;
    bounds facebounds = calcBoundingBox(faces);
    std::vector<faceWrapper*> faceswrapper = faceWrapper::toWrapperList(faces);
    std::vector<planarEvent*> events = buildEventList(faceswrapper);
    return new kdtree(events, facebounds, 0);
}

kdtree::kdtree(std::vector<planarEvent*>& events, const bounds& facebounds, int depth){
    if(depth==0) completion=0;
    this->boundingbox = facebounds;
    int nfaces=0;
    
    //face classifications are initially not 0 (see faceWrapper::toWrapperList), and are not 0 when passed from parent nodes. So an initial pass to reset classification values is not necessary
    for(planarEvent *e : events){
        if(e->fw->classification){
            nfaces++;
            e->fw->classification=0;
        }
    }

#define QUICKBUILD
#ifdef QUICKBUILD
    
    //exit if max recursion depth is reached or too few faces remaining
    if(depth > KDT_MAX_DEPTH || nfaces <= KDT_MIN_FACES || facebounds.volume()<EPSILON){
        //generate face list
        for(planarEvent *e: events){
            if(e->fw->classification==0){
                this->faces.push_back(e->fw->f);
                e->fw->classification=1;
            }
        }
        completion += 1.0/(1<<depth);
        return;
    }
#endif
    
    real cost = KDT_CTRAVERSAL + KDT_CINTERSECT * nfaces;
    real area = boundingbox.area();
    planedata optimalplane{};
    optimalplane = findOptimalPlane(nfaces, events, boundingbox, area);
    
    //only split if SAH finds it would reduce cost
    if(cost<=optimalplane.split_cost){
        //generate face list
        for(planarEvent *e: events){
            if(e->fw->classification==0){
                this->faces.push_back(e->fw->f);
                e->fw->classification=1;
            }
        }
        completion += 1.0/(1<<depth);
        return;
    }
    
    this->pos = optimalplane.plane_pos;
    this->axis = optimalplane.plane_axis;
    this->planarside = optimalplane.planar_side;
    bounds lowerbounds = boundingbox;
    bounds upperbounds = boundingbox;
    lowerbounds.max.vec[axis] = pos;
    upperbounds.min.vec[axis] = pos;
    
    std::vector<planarEvent*> totalevents_left, totalevents_right, oldevents_left, oldevents_right, newevents_left, newevents_right;
    std::vector<faceWrapper*> faces_bothsides;
    for(planarEvent *e: events) e->fw->classification = 3;
    for(planarEvent *e: events){
        if(e->axis==axis){
            if(e->type==0 && e->pos <= pos) e->fw->classification=2; //left only
            else if(e->type==2 && e->pos >= pos) e->fw->classification=1; //right only
            else if (e->type == 1){
                if(e->pos == pos) e->fw->classification = planarside ? 1 : 2;
                else if (e->pos > pos) e->fw->classification = 1;
                else e->fw->classification = 2;
            }
        }
    }
    for(planarEvent *e: events){
        if(e->fw->classification==1) oldevents_right.push_back(e);
        else if (e->fw->classification==2) oldevents_left.push_back(e);
        else if (e->fw->classification==3){
            faces_bothsides.push_back(e->fw);
            e->fw->classification=4;
        }
    }
    
    for(faceWrapper *fw : faces_bothsides) generateClippedEvents(fw, pos, axis, boundingbox, newevents_left, newevents_right);
    std::sort(newevents_left.begin(), newevents_left.end(), kdtree::planarEvent::planarEventComparator());
    std::sort(newevents_right.begin(), newevents_right.end(), kdtree::planarEvent::planarEventComparator());
    std::merge(oldevents_left.begin(), oldevents_left.end(), newevents_left.begin(), newevents_left.end(), std::back_inserter(totalevents_left), kdtree::planarEvent::planarEventComparator());
    std::merge(oldevents_right.begin(), oldevents_right.end(), newevents_right.begin(), newevents_right.end(), std::back_inserter(totalevents_right), kdtree::planarEvent::planarEventComparator());
    
    real lowerboundsarea = lowerbounds.area();
    real upperboundsarea = upperbounds.area();
    
    if(!totalevents_left.empty()) this->lower = new kdtree(totalevents_left, lowerbounds, depth+1);
    if(!totalevents_right.empty()) this->upper = new kdtree(totalevents_right, upperbounds, depth+1);
    
    //clean up events and face wrappers
    if(depth==0){
        std::vector<faceWrapper*> allfaces;
        for(planarEvent *e: events){
            if(e->fw->classification!=5) {
                e->fw->classification=5;
                allfaces.push_back(e->fw);
            }
            delete e;
        }
        for(faceWrapper *fw : allfaces) delete fw;
    }
    for(planarEvent* e : newevents_left) delete e;
    for(planarEvent* e : newevents_right) delete e;
}

/**
 Generates splitting plane events at the 6 planes of the bounding box of each face. These represent candidates for a splitting plane of the kdtere. The SAH algorithm finds the optimal plane out of these to split the tree
 */
std::vector<kdtree::planarEvent*> kdtree::buildEventList(const std::vector<faceWrapper*>& faces){
    std::vector<planarEvent*> events;
    real minx, maxx, miny, maxy, minz , maxz;
    for(faceWrapper *fw : faces){
        minx = fw->f->minCoord(0);
        maxx = fw->f->maxCoord(0);
        miny = fw->f->minCoord(1);
        maxy = fw->f->maxCoord(1);
        minz = fw->f->minCoord(2);
        maxz = fw->f->maxCoord(2);
        
        if(minx==maxx) events.push_back(new planarEvent(fw, minx, 0, 1));
        else{
            events.push_back(new planarEvent(fw, minx, 0, 2));
            events.push_back(new planarEvent(fw, maxx, 0, 0));
        }
        if(miny==maxy) events.push_back(new planarEvent(fw, miny, 1, 1));
        else{
            events.push_back(new planarEvent(fw, miny, 1, 2));
            events.push_back(new planarEvent(fw, maxy, 1, 0));
        }
        if(minz==maxz) events.push_back(new planarEvent(fw, minz, 2, 1));
        else{
            events.push_back(new planarEvent(fw, minz, 2, 2));
            events.push_back(new planarEvent(fw, maxz, 2, 0));
        }
    }
    std::sort(events.begin(), events.end(), planarEvent::planarEventComparator());
    return events;
}
/**
 * Finds optimal splitting plane in O(N) time.
 * From: http://www.eng.utah.edu/~cs6965/papers/kdtree.pdf
 *
 * @param nfaces
 *            The number of faces in the current bounding box
 * @param events
 *            A list of splitting plane candidates
 * @param boundingbox
 *            The current bounding box
 * @return A planedata struct containing data for the optimal splitting plane
 */

kdtree::planedata kdtree::findOptimalPlane(int nfaces, const std::vector<planarEvent*>& events, bounds& boundingbox, real area){
    
    int Nleft[3] = {0}, Nplanar[3] = {0}, Nright[3] = {nfaces, nfaces, nfaces};
    int axis_best = 0;
    real pos_best=0, cost_best = _INFINITY;
    bool p_side=0;
    
    //temp variables
    int plane_axis, start_count, planar_count, end_count;
    real plane_pos;
    planarEvent *e;
    std::pair<real, bool> planecostdata;
    
    for(int i=0; i<events.size();){
        e = events[i];
        plane_pos = e->pos;
        plane_axis = e->axis;
        start_count = planar_count = end_count = 0;
        while(i<events.size() && e->axis == plane_axis && e->pos == plane_pos && e->type==0){
            end_count++;
            i++;
            e = events[i];
        }
        while(i<events.size() && e->axis == plane_axis && e->pos == plane_pos && e->type==1){
            planar_count++;
            i++;
            e = events[i];
        }
        while(i<events.size() && e->axis == plane_axis && e->pos == plane_pos && e->type==2){
            start_count++;
            i++;
            e = events[i];
        }
        
        Nplanar[plane_axis] = planar_count;
        Nright[plane_axis] -= planar_count;
        Nright[plane_axis] -= end_count;
        
        planecostdata = SAH(boundingbox, plane_pos, plane_axis, area, Nleft[plane_axis], Nplanar[plane_axis], Nright[plane_axis]);
        if(planecostdata.first < cost_best){
            cost_best = planecostdata.first;
            pos_best = plane_pos;
            axis_best = plane_axis;
            p_side = planecostdata.second;
        }
        Nleft[plane_axis] += start_count;
        Nleft[plane_axis] += planar_count;
        Nplanar[plane_axis] = 0;
    }
    return planedata{pos_best, axis_best, cost_best, p_side};
}

/**
 * Computes approximate cost of traversing a ray through a voxel with a
 * given split plane, using the Surface Area Heuristic
 *
 * @param bounds
 *            The bounds of the whole voxel
 * @param position
 *            The position of the splitting plane
 * @param axis
 *            The axis of the splitting plane
 * @param area
 *            The area of the voxel
 * @param Nl
 *            The number of faces to the left (below) the splitting plane
 * @param Np
 *            The number of faces exactly (flat) on the splitting plane
 * @param Nr
 *            The number of faces to the right (above) the splitting plane
 * @return a double array {cost, pside}. pside = 0 or 1, depending if planar
 *         faces should be grouped in the left or right voxel, respectively
 */

std::pair<real, bool> kdtree::SAH(const bounds& boundingbox, real position, int axis, real area, int Nl, int Np, int Nr){
    bounds upperbounds = boundingbox;
    bounds lowerbounds = boundingbox;
    lowerbounds.max.vec[axis] = position;
    upperbounds.min.vec[axis] = position;
    
     //probability of a ray hitting each subvoxel, given that boundingbox is already intersected.
    real lowerarea = lowerbounds.area() / area;
    real upperarea = upperbounds.area() / area;
    //if bounding box has no area, or split cuts off just a plane,
    if(area==0 || boundingbox.d(axis)==0 || lowerarea==0 || upperarea==0) return std::pair<real, bool>(_INFINITY, false);
    
    real costleftbias = (lowerarea*(Nl+Np) + upperarea*Nr) * lambda(lowerarea, upperarea, Nl+Np, Nr);
    real costrightbias = lowerarea*Nl + upperarea*(Np+Nr) * lambda(lowerarea, upperarea, Nl, Np+Nr);
    //planar faces go into left or right depending on which results in lowest cost
    return std::pair<real, bool>(_min(costleftbias, costrightbias)*KDT_CINTERSECT+KDT_CTRAVERSAL, costrightbias<costleftbias);
}

#define ROUGH_FACE_SPLIT
#ifdef ROUGH_FACE_SPLIT
void kdtree::generateClippedEvents(kdtree::faceWrapper *fw, real pos, int axis, const bounds &boundingbox, std::vector<planarEvent *> &newevents_left, std::vector<planarEvent *> &newevents_right){
    bounds upperbounds = boundingbox;
    bounds lowerbounds = boundingbox;
    bounds facebound = fw->f->boundingBox();
    lowerbounds.max.vec[axis] = pos;
    upperbounds.min.vec[axis] = pos;
    intersectBoundingBoxes(lowerbounds, facebound, lowerbounds);
    intersectBoundingBoxes(upperbounds, facebound, upperbounds);
    
    for(int k=0; k<3; k++){
        if(lowerbounds.min.vec[k] == lowerbounds.max.vec[k]){
            newevents_left.push_back(new planarEvent(fw, lowerbounds.min.vec[k], k, 1));
        }
        else{
            newevents_left.push_back(new planarEvent(fw, lowerbounds.min.vec[k], k, 2));
            newevents_left.push_back(new planarEvent(fw, lowerbounds.max.vec[k], k, 0));
        }
        
        if(upperbounds.min.vec[k] == upperbounds.max.vec[k]){
            newevents_right.push_back(new planarEvent(fw, upperbounds.min.vec[k], k, 1));
        }
        else{
            newevents_right.push_back(new planarEvent(fw, upperbounds.min.vec[k], k, 2));
            newevents_right.push_back(new planarEvent(fw, upperbounds.max.vec[k], k, 0));
        }
    }
}
#else
void kdtree::generateClippedEvents(kdtree::faceWrapper *fw, real pos, int axis, const bounds &boundingbox, std::vector<planarEvent *> &newevents_left, std::vector<planarEvent *> &newevents_right){
    int sideflags = (fw->f->vertices[0]->vec[axis] < pos ? 0 : 1)
    | (fw->f->vertices[1]->vec[axis] < pos ? 0 : 2)
    | (fw->f->vertices[2]->vec[axis] < pos ? 0 : 4);
    
    //checks if sideflags is a power of two, ie there is only one vertex on the right
    bool oneUpperSide = (sideflags & (sideflags-1)) == 0;

    //whichever side has only one vertex, assign 'lonely' vertex to v1
    vertex *v1=nullptr, *v2=nullptr, *v3=nullptr;
    switch (oneUpperSide ? sideflags : (~sideflags) & 7) {
        case 1:
            v1 = fw->f->vertices[0];
            v2 = fw->f->vertices[1];
            v3 = fw->f->vertices[2];
            break;
        case 2:
            v1 = fw->f->vertices[1];
            v2 = fw->f->vertices[0];
            v3 = fw->f->vertices[2];
            break;
        case 4:
            v1 = fw->f->vertices[2];
            v2 = fw->f->vertices[0];
            v3 = fw->f->vertices[1];
            break;
    }
    real r12 = (pos-v1->vec[axis])/(v2->vec[axis]-v1->vec[axis]);
    real r13 = (pos-v1->vec[axis])/(v3->vec[axis]-v1->vec[axis]);
    //points on split plane along triangle edges
    vertex v12cut = vertex::lerp(*v1, *v2, r12);
    vertex v13cut = vertex::lerp(*v1, *v3, r13);

    
    int ax1 = (axis + 2) % 3, ax2 = (axis + 1) % 3; //other axes
    //clamp values of intersecting points to the bounding box
    clamp(&v12cut.vec[ax1], boundingbox.min.vec[ax1], boundingbox.max.vec[ax1]);
    clamp(&v12cut.vec[ax2], boundingbox.min.vec[ax2], boundingbox.max.vec[ax2]);
    clamp(&v13cut.vec[ax1], boundingbox.min.vec[ax1], boundingbox.max.vec[ax1]);
    clamp(&v13cut.vec[ax2], boundingbox.min.vec[ax2], boundingbox.max.vec[ax2]);
  
    bounds lowerbounds, upperbounds;
    if(oneUpperSide){
        vertex::listBounds(upperbounds, 3, v1, &v12cut, &v13cut);
        vertex::listBounds(lowerbounds, 4, v2, v3, &v12cut, &v13cut);
    }
    else {
        vertex::listBounds(lowerbounds, 3, v1, &v12cut, &v13cut);
        vertex::listBounds(upperbounds, 4, v2, v3, &v12cut, &v13cut);
    }
    //intersecting triangle bounding boxes with bounding box of volume
    vertex::intersectBoundingBoxes(lowerbounds, boundingbox, lowerbounds);
    vertex::intersectBoundingBoxes(upperbounds, boundingbox, upperbounds);
    
    for(int k=0; k<3; k++){
        if(lowerbounds.min.vec[k] == lowerbounds.max.vec[k]){
            newevents_left.push_back(new planarEvent(fw, lowerbounds.min.vec[k], k, 1));
        }
        else{
            newevents_left.push_back(new planarEvent(fw, lowerbounds.min.vec[k], k, 2));
            newevents_left.push_back(new planarEvent(fw, lowerbounds.max.vec[k], k, 0));
        }
        if(upperbounds.min.vec[k] == upperbounds.max.vec[k]){
            newevents_right.push_back(new planarEvent(fw, upperbounds.min.vec[k], k, 1));
        }
        else{
            newevents_right.push_back(new planarEvent(fw, upperbounds.min.vec[k], k, 2));
            newevents_right.push_back(new planarEvent(fw, upperbounds.max.vec[k], k, 0));
        }
    }
}
#endif

int kdtree::faceWrapper::i = 0;
kdtree::faceWrapper::faceWrapper(face *f, int classification) :f(f), classification(classification){i++;}
kdtree::faceWrapper::faceWrapper(face *f) :f(f){i++;}
kdtree::faceWrapper::~faceWrapper() {i--;}

std::vector<face*> kdtree::faceWrapper::toFaceList(const std::vector<faceWrapper*>& wrappedfaces){
    std::vector<face*> faces;
    for(faceWrapper* fw: wrappedfaces) faces.push_back(fw->f);
    return faces;
}
std::vector<kdtree::faceWrapper*> kdtree::faceWrapper::toWrapperList(const std::vector<face*>& faces){
    std::vector<faceWrapper*> wrappedfaces;
    for(face* f: faces) wrappedfaces.push_back(new faceWrapper(f, 1));
    return wrappedfaces;
}

int kdtree::planarEvent::i = 0;
kdtree::planarEvent::planarEvent(faceWrapper *fw, real p, int k, int type) :fw(fw), pos(p), axis(k), type(type){i++;}
kdtree::planarEvent::~planarEvent(){i--;}

bool kdtree::planarEvent::planarEventComparator::operator()(const kdtree::planarEvent* pe1, const kdtree::planarEvent* pe2){
    if(pe1->pos==pe2->pos){
        if(pe1->axis==pe2->axis) return pe1->type < pe2->type;
        return pe1->axis < pe2->axis;
    }
    return pe1->pos < pe2->pos;
}

face *kdtree::rayTreeIntersect(kdtree *kdt, const ray &r, bool lazy, vertex *tuv, bool isRoot){
    if (kdt == nullptr) return nullptr;
    vertex tuv1, tuv2;
    bool tuvvalid = (tuv!=nullptr); //in case no variable to store intersection coordinates is given
    face *f1, *f2;
    
    //does ray intersect this bounding box?
    if(rayAABBIntersect(kdt->boundingbox, r)){
        //leaf node or nah?
        if(kdt->lower!=nullptr || kdt->upper!=nullptr){
            f1 = rayTreeIntersect(kdt->lower, r, lazy, &tuv1, true);
            if(lazy && f1!=nullptr){
                if(tuvvalid) *tuv = tuv1;
                return f1;
            }
            f2 = rayTreeIntersect(kdt->upper, r, lazy, &tuv2, true);
            if(lazy && f2!=nullptr){
                if(tuvvalid) *tuv = tuv2;
                return f2;
            }
            
            //non-lazy compare return face from each child
            if(f1==nullptr){
                if(tuvvalid) *tuv = tuv2;
                return f2;
            }
            else if (f2==nullptr){
                if(tuvvalid) *tuv = tuv1;
                return f1;
            }
            else{
                if(tuv1.t < tuv2.t){
                    if(tuvvalid) *tuv = tuv1;
                    return f1;
                }
                else{
                    if(tuvvalid) *tuv = tuv2;
                    return f2;
                }
            }
        }
        else{
            //std::cout << "intersect " << kdt->faces.size() << " faces" << std::endl;
            f1 = rayFacesIntersect(kdt->faces, r, lazy, &tuv1);
            if(tuvvalid) *tuv = tuv1;
            return f1;
        }
    }
    return nullptr;
}

std::vector<edge*> kdtree::wireframe(){
    return wireframe(true, 0);
}

std::vector<edge*> kdtree::wireframe(bool isRoot, int depth){
    std::vector<edge*> edges;
    if(isRoot){
        
        //create corner vertices (capital letters indicate greater magnitude along given axis)
        vertex *xyz = new vertex(boundingbox.min);
        vertex *XYZ = new vertex(boundingbox.max);
        vertex *xyZ = new vertex(xyz->x, xyz->y, XYZ->z);
        vertex *xYz = new vertex(xyz->x, XYZ->y, xyz->z);
        vertex *xYZ = new vertex(xyz->x, XYZ->y, XYZ->z);
        vertex *Xyz = new vertex(XYZ->x, xyz->y, xyz->z);
        vertex *XyZ = new vertex(XYZ->x, xyz->y, XYZ->z);
        vertex *XYz = new vertex(XYZ->x, XYZ->y, xyz->z);
    
        //from bottom corner
        edges.push_back(new edge(xyz, xyZ));
        edges.push_back(new edge(xyz, xYz));
        edges.push_back(new edge(xyz, Xyz));
        
        //from top corner
        edges.push_back(new edge(XYZ, XYz));
        edges.push_back(new edge(XYZ, XyZ));
        edges.push_back(new edge(XYZ, xYZ));
        
        //connect remaining edges
        edges.push_back(new edge(xyZ, xYZ));
        edges.push_back(new edge(xyZ, XyZ));
        edges.push_back(new edge(xYz, XYz));
        edges.push_back(new edge(xYz, xYZ));
        edges.push_back(new edge(Xyz, XYz));
        edges.push_back(new edge(Xyz, XyZ));
    }
    
    if(faces.empty()){
        switch (axis) {
            case 0:
                edges.push_back(new edge(new vertex(pos, boundingbox.min.y, boundingbox.min.z), new vertex(pos, boundingbox.min.y, boundingbox.max.z)));
                edges.push_back(new edge(new vertex(pos, boundingbox.min.y, boundingbox.min.z), new vertex(pos, boundingbox.max.y, boundingbox.min.z)));
                edges.push_back(new edge(new vertex(pos, boundingbox.min.y, boundingbox.max.z), new vertex(pos, boundingbox.max.y, boundingbox.max.z)));
                edges.push_back(new edge(new vertex(pos, boundingbox.max.y, boundingbox.min.z), new vertex(pos, boundingbox.max.y, boundingbox.max.z)));
                break;
            case 1:
                edges.push_back(new edge(new vertex(boundingbox.min.x, pos, boundingbox.min.z), new vertex(boundingbox.min.x, pos, boundingbox.max.z)));
                edges.push_back(new edge(new vertex(boundingbox.min.x, pos, boundingbox.min.z), new vertex(boundingbox.max.x, pos, boundingbox.min.z)));
                edges.push_back(new edge(new vertex(boundingbox.min.x, pos, boundingbox.max.z), new vertex(boundingbox.max.x, pos, boundingbox.max.z)));
                edges.push_back(new edge(new vertex(boundingbox.max.x, pos, boundingbox.min.z), new vertex(boundingbox.max.x, pos, boundingbox.max.z)));
                break;
            case 2:
                edges.push_back(new edge(new vertex(boundingbox.min.x, boundingbox.min.y, pos), new vertex(boundingbox.min.x, boundingbox.max.y, pos)));
                edges.push_back(new edge(new vertex(boundingbox.min.x, boundingbox.min.y, pos), new vertex(boundingbox.max.x, boundingbox.min.y, pos)));
                edges.push_back(new edge(new vertex(boundingbox.max.x, boundingbox.min.y, pos), new vertex(boundingbox.max.x, boundingbox.max.y, pos)));
                edges.push_back(new edge(new vertex(boundingbox.min.x, boundingbox.max.y, pos), new vertex(boundingbox.max.x, boundingbox.max.y, pos)));
                break;
        }
        if(lower!=nullptr){
            std::vector<edge*> loweredges = lower->wireframe(false, depth+1);
            edges.insert(edges.end(), loweredges.begin(), loweredges.end());
        }
        if(upper!=nullptr){
            std::vector<edge*> upperedges = upper->wireframe(false, depth+1);
            edges.insert(edges.end(), upperedges.begin(), upperedges.end());
        }
    }
    return edges;
}

void kdtree::calcstats(int& nodes, int& leaf_nodes, int& non_empty_leaf_nodes, int& triangles_in_leaves, real& est_traversals, real& est_leaves_visited, real& est_tris_intersected, real& est_cost, real rootArea){
    if(rootArea<0) rootArea = this->boundingbox.area();
    real area_ratio = this->boundingbox.area()/rootArea;
    
    nodes++;
    est_traversals += area_ratio;
    est_cost += area_ratio * 15;
    
    if(!this->faces.empty()){
        //every leaf node is non empty;
        leaf_nodes++;
        non_empty_leaf_nodes++;
        triangles_in_leaves+=this->faces.size();
        est_leaves_visited += area_ratio;
        est_tris_intersected += area_ratio * this->faces.size();
        est_cost += area_ratio*20;
        //constants used for est_cost (15, 20) are used to match up to Ingo Wald 2006 paper, to compare for debugging purposes
    }
    
    if(this->lower!=nullptr) this->lower->calcstats(nodes, leaf_nodes, non_empty_leaf_nodes, triangles_in_leaves, est_traversals, est_leaves_visited, est_tris_intersected, est_cost, rootArea);
    if(this->upper!=nullptr) this->upper->calcstats(nodes, leaf_nodes, non_empty_leaf_nodes, triangles_in_leaves, est_traversals, est_leaves_visited, est_tris_intersected, est_cost, rootArea);
}

void kdtree::printstats(){
    clock_t begin = clock();
    for(int i=1; i<=100000; i++){
        rayTreeIntersect(this, ray(vertex(), vertex((real)rand()/RAND_MAX, (real)rand()/RAND_MAX, (real)rand()/RAND_MAX)), false, nullptr);
    }
    clock_t end = clock();
    std::cout << "100,000 kdt intersects: " << real(end - begin) / CLOCKS_PER_SEC << std::endl;

    int n=0, ln=0, ln_ne=0, t_tot=0;
    real est_t=0, est_l=0, est_tr=0, est_c=0;
    calcstats(n, ln, ln_ne, t_tot, est_t, est_l, est_tr, est_c);
    std::cout << std::endl << "total nodes: " << n << std::endl;
    std::cout << "leaf nodes : " << ln << std::endl;
    std::cout << "non-empty leaf nodes: " << ln_ne << std::endl;
    std::cout << "avg. faces per leaf node: " << real(t_tot)/ln_ne << std::endl;
    std::cout << std::endl << "Average, for each ray:" << std::endl;
    std::cout << "est. num of traversals: " << est_t << std::endl;
    std::cout << "est. num of leaf nodes visited: " << est_l << std::endl;
    std::cout << "est. num of triangles intersected: " << est_tr << std::endl;
    std::cout << "est. cost (using SAH): " << est_c << std::endl;
}
