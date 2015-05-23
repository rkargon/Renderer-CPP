//
//  renderthreads.cpp
//  Renderer
//
//  Created by Raphael Kargon on 5/16/15.
//  Copyright (c) 2015 Raphael Kargon. All rights reserved.
//

#include "renderthreads.h"

const int thread_manager::tilesize = 32;

thread_manager::thread_manager(const int nthreads, raster **raster_ref, scene *sc)
:raster_ref(raster_ref), sc(sc){
    //set up thread pool
    thread_pool = std::vector<worker_thread>();
    for(int i=0; i<nthreads; i++){
        thread_pool.push_back(worker_thread(this));
    }
    
    tiles = std::queue<tile>();
    fill_tile_queue();
}

void thread_manager::start(){
    for(worker_thread &w : thread_pool){
        w.start();
    }
}

void thread_manager::stop(){
    for(worker_thread &w : thread_pool){
        w.stop();
    }
}

void thread_manager::set_render_method(color (*render_method)(const int, const int, const int, const int, scene*)){
    this->render_method = render_method;
}

void thread_manager::fill_tile_queue(){
    int w = (*raster_ref)->width();
    int h = (*raster_ref)->height();
    
    //clear tiles
    this->tile_queue_lock.lock();
    std::queue<tile>().swap(tiles);
    
    //fill queue with new tiles
    for(int x = 0; x < w; x += tilesize){
        for(int y = 0; y < h; y += tilesize){
            //std::cout << x << " " << y << " " << std::min(w - x, tilesize) << " " << std::min(h - y, tilesize) << " " << std::endl;
            tiles.push(tile{x, y, std::min(w - x, tilesize), std::min(h - y, tilesize)});
        }
    }
    this->tile_queue_lock.unlock();
}

worker_thread::worker_thread(thread_manager *manager){
    this->manager = manager;
    this->is_running = false;
}

void worker_thread::start() {
    if(!is_running){
        is_running = true;
        render_thread = std::thread(render_tiles, std::ref(manager->tiles), std::ref(manager->tile_queue_lock), std::ref(is_running), manager->raster_ref, manager->sc, manager->render_method);
        render_thread.detach();
    }
}

void worker_thread::stop(){
    if(is_running){
        is_running = false;
    }
}

void render_tiles(std::queue<tile> &tiles, std::mutex &tile_queue_lock, bool &is_running, raster **raster_ref, scene *sc, color (*render_method)(const int, const int, const int, const int, scene*)){
    while(is_running){
        //get next tile (or stop if no tiles left)
        tile current_tile;
        tile_queue_lock.lock();
        if(!tiles.empty()){
            current_tile = tiles.front();
            tiles.pop();
            tile_queue_lock.unlock();
        }
        else{
            tile_queue_lock.unlock();
            return;
        }
        
        int img_width = (*raster_ref)->width();
        int img_height = (*raster_ref)->height();
        //check if tile is valid (ie in case of image resize)
        if(current_tile.x+current_tile.w > img_width || current_tile.y + current_tile.h > img_height){
            continue;
        }
        
        //render current tile
        for(int x = current_tile.x; x < current_tile.x + current_tile.w; x++){
            for(int y = current_tile.y; y < current_tile.y + current_tile.h; y++){
                if(!is_running){
                    return;
                }
                color col = render_method(x, y, img_width, img_height, sc);
                (*raster_ref)->colbuffer[y*img_width + x] = colorToRGB(col);
            }
        }
    }
}