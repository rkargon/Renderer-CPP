//
//  renderthreads.h
//  Renderer
//
//  Created by Raphael Kargon on 5/16/15.
//  Copyright (c) 2015 Raphael Kargon. All rights reserved.
//

#ifndef __Renderer__renderthreads__
#define __Renderer__renderthreads__

#include <queue>
#include <stdio.h>
#include <thread>
#include <vector>
#include "raster.h"
#include "scene.h"

class thread_manager;
class worker_thread;

typedef struct tile {
    int x, y;
    int w, h;
} tile;

using render_method_func=std::function<color (const double, const double, const int, const int, scene*)>;

class thread_manager {
public:
    static const int tilesize;
    
    thread_manager() = delete;
    thread_manager(const int nthreads, raster **raster_ref, scene *sc, const int antialiasing_per_side = 1);
    void start();
    void stop();
    void fill_tile_queue();
    void set_render_method(render_method_func render_method);
    
    friend class worker_thread;
    
private:
    std::vector<worker_thread> thread_pool;
    std::queue<tile> tiles;
    std::mutex tile_queue_lock;
    
    //a double pointer, in case actual pointer changes when image is resized
    int antialiasing_per_side;
    raster **raster_ref;
    scene *sc;
    render_method_func render_method;
};

class worker_thread {
public:
    worker_thread(thread_manager *manager);
    void start();
    void stop();
    void render_tiles();

    
private:
    bool is_running;
    std::thread render_thread;
    thread_manager *manager;
};

#endif /* defined(__Renderer__renderthreads__) */
