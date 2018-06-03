#ifndef __Renderer__sampling__
#define __Renderer__sampling__

#include "geom.h"

#include <random>
#include <utility>

// per-thread RNG
static thread_local std::mt19937 generator;

typedef std::pair<vertex, double> direction_sample;

direction_sample uniform_unit_sphere();

// Returns a vector corresponding to a random point on the unit hemisphere
// in
// the same direction as the given normal vector
direction_sample uniform_unit_hemisphere(const vertex &normal);

#endif // __Renderer__sampling__
