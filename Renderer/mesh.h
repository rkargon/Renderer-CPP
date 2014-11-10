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
enum tex_projection_t {TEX_PROJ_SPHERICAL, TEX_PROJ_CUBIC, TEXT_PROJ_CYLINDRICAL};

class mesh
{
public:
    std::vector<meshvertex*> vertices; //List of *unique* vertices in mesh
    std::vector<edge*> edges;//List of *unique* edges in mesh
    std::vector<face*> faces;//List of faces in mesh
    std::string name;
    material *mat;
    BSDF *bsdf;
    bool smooth;
    
    /* mesh(std::ifstream& infile, std::string objname)
     
     This constructor reads a binary STL file and creates a mesh based off of the vertices and faces in the file. 
     If infile is invalid, no vertices, faces, or edges will be read.
     By default, this>smooth is set to false, and this->mat is uninitialized. 
     this->name is set to 'objname'.
     
     Parameters:
     std::ifstrem& infile : A binary STIL file from which geomtry is loaded.
     std::string objname : An identifier for the object, stored in this->name
     */
    mesh(std::ifstream& infile, std::string objname);
    
    
    /* void project_texture(tex_projection_t proj)
     Assigns texture coordinates to each vertex in this->vertices based on a given projection specifed 
      by proj.
     
     Parameters:
     proj - Specifies which type of projection to use:
        TEX_PROJ_SPHERICAL
            Projects vertices onto a sphere and uses spherical coordinates as texture coordinates.
            Each vertex is converted from cartesian to spherical coordinates which are normalized to [0, 1]
        TEX_PROJ_CUBIC
        TEXT_PROJ_CYLINDRICAL
     */
    void project_texture(tex_projection_t proj);
};

#endif
