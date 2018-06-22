//
//  renderthreads.cpp
//  Renderer
//
//  Created by Raphael Kargon on 5/16/15.
//  Copyright (c) 2015 Raphael Kargon. All rights reserved.
//

#include "renderthreads.h"

thread_manager::thread_manager(raster *raster_ref, const scene *sc,
                               const render_options &opts)
    : opts(opts), raster_ref(raster_ref), sc(sc) {
  this->is_running = false;
  if (this->opts.threads == 0) {
    // Use one less thread than max to allow other programs to run on system.
    this->opts.threads =
        std::max<unsigned>(1, (std::thread::hardware_concurrency() - 1));
  }
  for (unsigned int i = 0; i < this->opts.threads; i++) {
    this->thread_pool.push_back(
        std::thread(&thread_manager::render_tiles, this));
  }
}

thread_manager::~thread_manager() {
  this->stop();
  for (auto &th : thread_pool) {
    th.join();
  }
}

// TODO worker threads might finish current tile before actually re-starting
void thread_manager::start() {
  if (!is_running) {
    // Signal threads to start running
    {
      std::lock_guard<std::mutex> lk(this->is_running_mutex);
      this->is_running = true;
    }
    this->is_running_cv.notify_all();
    this->fill_tile_queue();
  }
}

void thread_manager::stop() { this->is_running = false; }

void thread_manager::set_render_method(render_method_func render_method) {
  this->render_method = render_method;
}

void thread_manager::fill_tile_queue() {
  unsigned int w = raster_ref->width();
  unsigned int h = raster_ref->height();

  std::unique_lock<std::mutex> lk(this->tile_queue_lock);

  // clear tiles
  std::queue<tile>().swap(tiles);

  // fill queue with new tiles
  for (unsigned int x = 0; x < w; x += opts.tilesize) {
    for (unsigned int y = 0; y < h; y += opts.tilesize) {
      tiles.push(
          tile{x, y, std::min(static_cast<unsigned int>(w - x), opts.tilesize),
               std::min(static_cast<unsigned int>(h - y), opts.tilesize)});
    }
  }
  lk.unlock();
  this->queue_filled_cv.notify_all();
}

tile thread_manager::get_next_tile() {
  // Wait for tile queue to not be empty
  std::unique_lock<std::mutex> lk(this->tile_queue_lock);
  this->queue_filled_cv.wait(lk, [this] { return !this->tiles.empty(); });
  // Get next tile and remove it from queue
  tile current_tile = this->tiles.front();
  this->tiles.pop();
  return current_tile;
  // Mutex lock automatically unlocks when exiting scope
}

void thread_manager::render_single_tile(const tile &t) {
  unsigned int img_width = raster_ref->width();
  unsigned int img_height = raster_ref->height();
  // check if tile is valid (ie in case of image resize)
  if (t.x + t.w > img_width || t.y + t.h > img_height) {
    return;
  }

  // render current tile
  double antialias_delta = 1.0 / opts.antialiasing_per_side;
  double antialias_factor = antialias_delta * antialias_delta;
  for (unsigned int x = t.x; x < t.x + t.w; x++) {
    for (unsigned int y = t.y; y < t.y + t.h; y++) {
      color col = {0, 0, 0};
      for (unsigned int i = 0; i < opts.antialiasing_per_side; i++) {
        for (unsigned int j = 0; j < opts.antialiasing_per_side; j++) {
          if (!is_running) {
            return;
          }
          col += antialias_factor *
                 render_method(x + i * antialias_delta, y + j * antialias_delta,
                               img_width, img_height, *sc, opts);
        }
      }
      raster_ref->colbuffer[y * img_width + x] = color_to_rgb(col);
    }
  }
}

void thread_manager::render_tiles() {
  while (true) {
    // Wait until signalled that thread should be running
    std::unique_lock<std::mutex> lk(this->is_running_mutex);
    this->is_running_cv.wait(lk, [this] { return this->is_running; });
    lk.unlock();
    // Waits until there is a tile available in the queue to render
    tile current_tile = this->get_next_tile();
    this->render_single_tile(current_tile);
  }
}
