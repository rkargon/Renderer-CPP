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
#include "scene.h"

class mesh
{
public:
    std::vector<meshvertex*> vertices;
    std::vector<edge*> edges;
    std::vector<face*> faces;
    std::string name;
    material *mat;
    bool smooth;
    mesh(std::ifstream& infile, std::string objname);
    void project_texture(tex_projection_t proj);
};

#endif
