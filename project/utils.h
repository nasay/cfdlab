#pragma once

#include "LBDefinitions.h"
#include <vector>
#include <type_traits>

template <class T>
struct Vector_3D {
  T x;
  T y;
  T z;
  bool operator== (const Vector_3D<T> &other) const {
    return x == other.x && y == other.y && z == other.z;
  }
  Vector_3D<T> operator+ (const Vector_3D<T> &rhs) {
    return { this->x + rhs.x, this->y + rhs.y, this->z + rhs.z };
}
};


typedef Vector_3D<int> Cell;
typedef Vector_3D<float> Velocity;

static const Velocity LATTICEVELOCITIES[Q] = {
                                              {0, -1, -1}, {-1, 0, -1}, {0, 0, -1}, {1, 0, -1}, {0, 1, -1},
                                              {-1, -1, 0}, {0, -1, 0}, {1, -1, 0}, {-1, 0, 0}, {0, 0, 0}, {1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {1, 1, 0},
                                              {0, -1, 1}, {-1, 0, 1}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}
};

static const float LATTICEWEIGHTS[Q] = {
                                        1.0/36, 1.0/36, 2.0/36, 1.0/36, 1.0/36,
                                        1.0/36, 2.0/36, 1.0/36, 2.0/36, 12.0/36, 2.0/36, 1.0/36, 2.0/36, 1.0/36, 
                                        1.0/36, 1.0/36, 2.0/36, 1.0/36, 1.0/36
};

struct Config {
  int length[3];
  int timesteps;
  int timestepsPerPlotting;
  int boundaries[6];
  int r;
  int t;
  int n_threads;
  int i;

  char problem[80];

  float tau;
  float extForces[3];
  float exchange;
  Velocity velocity;
  float ro_in;
  float ro_ref;
};

struct Fields {
  /* Collide field */
  std::vector<float> collide;

  /* Streaming field */
  std::vector<float> stream;

  /* Flag field */
  std::vector<int> flag;

  /* Mass field */
  std::vector<float> mass;

  /* fluid fraction field */
  std::vector<float> fraction;

  Fields (int n[3]) {
    collide.resize (Q * n[0] * n[1] * n[2]);
    stream.resize (Q * n[0] * n[1] * n[2]);
    flag.resize (n[0] * n[1] * n[2]);
    mass.resize (n[0] * n[1] * n[2]);
    fraction.resize (n[0] * n[1] * n[2]);
  }

  Fields (Fields &other) = delete;
  Fields (const Fields &other) = delete;
};
