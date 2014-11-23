//
//  mesh.h
//  Renderer
//
//  Created by Raphael Kargon on 6/5/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#ifndef __Renderer__mesh__
#define __Renderer__mesh__

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include "geom.h"
#include "BSDF.h"

class material;

/* Texture projections:
 Determine how a texture is mapped onto an object's faces. 
 
TEX_PROJ_SPHERICAL
    Projects vertices onto a sphere and uses spherical coordinates as texture coordinates.
    Each vertex is converted from cartesian to spherical coordinates which are normalized to [0, 1]
TEX_PROJ_CUBIC
TEXT_PROJ_CYLINDRICAL
 */
enum tex_projection_t {TEX_PROJ_SPHERICAL, TEX_PROJ_CUBIC, TEXT_PROJ_CYLINDRICAL};

class mesh
{
public:
    std::vector<meshvertex*> vertices; //List of *unique* vertices in mesh
    std::vector<edge*> edges;//List of *unique* edges in mesh
    std::vector<face*> faces;//List of faces in mesh
    vertex origin;
    std::string name;
    material *mat;
    BSDF *bsdf;
    bool smooth;
    
    /* mesh(std::ifstream& infile, std::string objname)
     
     This constructor reads a binary STL file and creates a mesh based off of the vertices and faces in the file. 
     If infile is invalid, no vertices, faces, or edges will be read.
     By default, this>smooth is set to false, and this->mat is uninitialized. 
     this->name is set to 'objname'.
     */
    mesh(std::ifstream& infile, std::string objname);
    
    
    /* void project_texture(tex_projection_t proj)
     Assigns texture coordinates to each vertex in this->vertices based on a given projection specifed 
      by proj.
     */
    void project_texture(tex_projection_t proj);
    
    /* Mesh manipulation functions 
        NOTE: These functions INVALIDATE kd-trees that have beeen constructed.
     */
    
    // moves an object (ie moves all vertices in an object)
    void move(const vertex& dv);
    
    // scales an object about a given center, by a given amount.
    // the scaling amount is given as a vector, where each component represents the scaling amount along
    // a corresponding axis
    void scale(const vertex& ds, const vertex& center);
    
    //scales an object relative to the geometric mean of its vertices.
    //This function calls scale(const vertex&, const vertex&) to do the actual scaling.
    void scale_centered(const vertex& ds);

    //calculates the centroid, or arithmetic mean, of all the vertices in an object. 
    vertex centroid();
};

#endif
