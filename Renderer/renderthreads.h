//
//  renderthreads.h
//  Renderer
//
//  Created by Raphael Kargon on 5/16/15.
//  Copyright (c) 2015 Raphael Kargon. All rights reserved.
//

#ifndef __Renderer__renderthreads__
#define __Renderer__renderthreads__

#include <algorithm>
#include <condition_variable>
#include <queue>
#include <stdio.h>
#include <thread>
#include <vector>
#include "raster.h"
#include "scene.h"

class thread_manager;

// Represents a rectangular tile of the image
typedef struct tile {
    // top left pixel coordinates
    int x, y;
    // width and height of tile
    int w, h;
} tile;

/**
 *  Type signature for a rendering function used by individual threads.
 *  This takes coordinates for a point and scene information and returns the color of the corresponding point.
 *
 *  @param x  X coordinate of the pixel in question
 *  @param y  Y coordinate of the pixel in question
 *  @param w  The width of the image
 *  @param h  The height of the image
 *  @param sc Pointer to scene data
 *
 *  @return Color of the resulting image at the given point
 */
// TODO declare render methods with this signature?
using render_method_func=std::function<color (const double x, const double y, const int w, const int h, scene* sc)>;

/**
 *  Manages a pool of rendering threads. Each thread renders one tile of an image at a time from a queue of tiles. 
 *  If the queue is empty, or is_running is false, the threads sleep.
 */
class thread_manager {
public:
    thread_manager() = delete;

    /**
     *  Creates a rendering thread pool with the given parameters. 
     * 
     *  @param raster_ref            Pointer to pointer to image rasters -- Rendered pixels are rendered to this.
     *  @param sc                    Pointer to scene information
     *  @param antialiasing_per_side Antialiasing samples along each size of a pixel. 
     *                               Thus, the number of total samples is the square of this.
     *                               Resulting pixel color is the average of all samples for that pixel. 
     *  @param nthreads              Number of rendering threads in pool. If 0, uses system information to determine best number of threads. 
     *  @param tilesize              Width and height in pixels of each tile.
     */
    thread_manager(raster **raster_ref, scene *sc, const int antialiasing_per_side = 1, int nthreads = 0, const int tilesize = 32);
    
    /**
     *  Starts the thread pool, and re-sets the tile queue.
     */
    void start();
    
    /**
     *  Tells the thread pool to stop rendering.
     */
    void stop();
    
    /**
     *  Sets the rendering function used by the thread pool.
     *
     *  @param render_method A function that renders a single pixel and returns a color. 
     */
    void set_render_method(render_method_func render_method);
    
private:
    std::vector<std::thread> thread_pool;
    std::queue<tile> tiles;
    std::mutex tile_queue_lock;
    std::condition_variable queue_filled_cv;
    
    bool is_running;
    std::mutex is_running_mutex;
    std::condition_variable is_running_cv;
    
    int tilesize;
    int antialiasing_per_side;
    //a double pointer, in case actual pointer changes when image is resized
    raster **raster_ref;
    scene *sc;
    render_method_func render_method;
    
    void fill_tile_queue();
    tile get_next_tile();
    void render_tiles();
    void render_single_tile(const tile& t);
};


#endif /* defined(__Renderer__renderthreads__) */
