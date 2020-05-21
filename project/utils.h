#pragma once

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


