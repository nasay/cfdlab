#pragma once

#include "LBDefinitions.h"
#include <vector>

struct Cell {
  int x;
  int y;
  int z;
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

  char *problem;

  float tau;
  float extForces[3];
  float exchange;
  float velocity[3];
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
};
