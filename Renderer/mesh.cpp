//
//  mesh.cpp
//  Renderer
//
//  Created by Raphael Kargon on 6/5/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#include "mesh.h"
using namespace std;

mesh::mesh(ifstream& infile, string objname){
    smooth = false; //not smooth by default
    name = objname;
    cout << "Loading STL object \"" << name << "\"" << endl;
    if(infile.fail()){
        cout << "Reading of " << name << " failed, creating an empty object." << endl;
        return;
    }
    
    //read header
    char *header = new char[80];
    infile.read(header, 80);
    cout << header << endl;
    
    char *facebuffer = new char[50];
    infile.read(facebuffer, 4);
    
    //read face data
    int i=0;
    float xtmp, ytmp, ztmp;
    vertex norm, vtmp;
    meshvertex *v1, *v2, *v3;
    face *facetmp;
    unordered_map<vertex, meshvertex*, vertexHasher> vertices_hash;
    unordered_map<vertex, meshvertex*, vertexHasher>::const_iterator v_iter;
    unordered_set<edge, edgeHasher> edges_hash;
    unordered_set<edge, edgeHasher>::const_iterator e_iter;
    while (!infile.read(facebuffer, 50).eof()){
        //read face normal
        xtmp = ((float *)facebuffer)[0];
        ytmp = ((float *)facebuffer)[1];
        ztmp = ((float *)facebuffer)[2];
        norm.set(xtmp, ytmp, ztmp);
        
        //read first vertex
        xtmp = ((float *)facebuffer)[3];
        ytmp = ((float *)facebuffer)[4];
        ztmp = ((float *)facebuffer)[5];
        vtmp.set(xtmp, ytmp, ztmp);
        //if vertex is not unique, then just point to the existing vertex.
        //Otherwise, store new vertex
        if((v_iter = vertices_hash.find(vtmp)) == vertices_hash.end()){
            v1 = new meshvertex(xtmp, ytmp, ztmp);
            vertices_hash.insert({vtmp, v1});
        }
        else v1 = v_iter->second;
        
        //second vertex
        xtmp = ((float *)facebuffer)[6];
        ytmp = ((float *)facebuffer)[7];
        ztmp = ((float *)facebuffer)[8];
        vtmp.set(xtmp, ytmp, ztmp);
        if((v_iter = vertices_hash.find(vtmp)) == vertices_hash.end()){
            v2 = new meshvertex(xtmp, ytmp, ztmp);
            vertices_hash.insert({vtmp, v2});
        }
        else v2 = v_iter->second;
        
        //third vertex
        xtmp = ((float *)facebuffer)[9];
        ytmp = ((float *)facebuffer)[10];
        ztmp = ((float *)facebuffer)[11];
        vtmp.set(xtmp, ytmp, ztmp);
        if((v_iter = vertices_hash.find(vtmp)) == vertices_hash.end()){
            v3 = new meshvertex(xtmp, ytmp, ztmp);
            vertices_hash.insert({vtmp, v3});
        }
        else v3 = v_iter->second;
        
        //Store edges of triangle
        edges_hash.emplace(v1, v2);
        edges_hash.emplace(v2, v3);
        edges_hash.emplace(v3, v1);
        
        //Load face with normal, vertices, and link to this object
        facetmp = new face(norm, v1, v2, v3, this);
        faces.push_back(facetmp);
        v1->faces.push_back(facetmp);
        v2->faces.push_back(facetmp);
        v3->faces.push_back(facetmp);
        
        i++;
    }
    
    //Load unique vertices into this->vertices
    for(auto it = vertices_hash.begin(); it != vertices_hash.end(); ++it){
        vertices.push_back(it->second);
    }
    //load unique edges into this->edges
    for(auto it = edges_hash.begin(); it != edges_hash.end(); ++it){
        edges.push_back(new edge(it->v1, it->v2));
    }
    
    //generate object origin
    origin = centroid();
    
    cout << vertices.size() << " vertices." << endl;
    cout << edges.size() << " edges." << endl;
    cout << faces.size() << " faces." << endl;
    
    delete header;
    delete facebuffer;
}

void mesh::project_texture(tex_projection_t proj){
    vertex vn;
    for(meshvertex *v : vertices){
        vn = v->unitvect();
        switch (proj) {
            case TEX_PROJ_SPHERICAL:
                v->tex_u = _atan2(vn.x,vn.y);
                if(v->tex_u<0) v->tex_u += M_PI*2;
                v->tex_u /= (2*M_PI);
                v->tex_v = _acos(vn.z)/M_PI;
                break;
                
            default:
                break;
        }
    }
}

void mesh::move(const vertex& dv){
    for(vertex *v: vertices){
        (*v) += dv;
    }
    origin += dv;
}

void mesh::scale(const vertex& ds, const vertex& scale_center){
    vertex dv;
    for(vertex *v: vertices){
        dv = *v - scale_center; //get vertex relative to center
        dv *= ds; //scale vertex relative to center
        *v = scale_center + dv; //return new vertex
    }
    dv = origin - scale_center;
    dv *= ds;
    origin = scale_center + dv;
}
void mesh::scale_centered(const vertex& ds){
    scale(ds, origin);
}

vertex mesh::centroid(){
    vertex mean;
    for(vertex* v: vertices){
        mean += v;
    }
    return mean * (1.0/vertices.size());
}