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

class thread_manager {
public:
    static const int tilesize;
    
    thread_manager() = delete;
    thread_manager(const int nthreads, raster **raster_ref, scene *sc);
    void start();
    void stop();
    void fill_tile_queue();
    void set_render_method(color (*render_method)(const int, const int, const int, const int, scene*));
    
    friend class worker_thread;
    
private:
    std::vector<worker_thread> thread_pool;
    std::queue<tile> tiles;
    std::mutex tile_queue_lock;
    
    //a double pointer, in case actual pointer changes when image is resized
    raster **raster_ref;
    scene *sc;
    color (*render_method)(const int, const int, const int, const int, scene*);
};

class worker_thread {
public:
    worker_thread(thread_manager *manager);
    void start();
    void stop();
    
private:
    bool is_running;
    std::thread render_thread;
    thread_manager *manager;
};

void render_tiles(std::queue<tile> &tiles, std::mutex &tile_queue_lock, bool &is_running, raster **raster_ref, scene *sc, color (*render_method)(const int, const int, const int, const int, scene*));

#endif /* defined(__Renderer__renderthreads__) */
