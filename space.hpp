#ifndef SPACE_H
#define SPACE_H

#include "state.hpp"
#include <vector>
#include <iostream>

class Space {
private:
  int m_;
  int n_;
  int M_;
  std::vector<State> states_;

public:
  double lambdax1;
  double lambdax2;
  double lambday;

  double mu;
  double gamma;

  // Back order costs
  double bx1;
  double bx2;
  double by;

  // Holding costs
  double hx1;
  double hx2;
  double hy;

  // Bulk order
  int Dx1;
  int Dx2;
  int Dy;

  double pn;
  double po;

  double alpha;

  /* Creates a new Markov Space object with the given bounds
   *
   * Args: n = minimum inventory to simulate
   *       m = maximum inventory to simulate
   *       M = maximum server population
   */
  Space(int n, int m, int M);

  // Gets the State at the given indices (by reference)
  State& state(int N, int x1, int x2, int y);

  // Holding cost function for current x, y
  double h(int x1, int x2, int y);

  // Bulk size probability
  double px1(int d);
  double px2(int d);
  double  py(int d);

  // Calculate T??v
  double T1v(int x1, int x2, int y, int N);
  double T2v(int x1, int x2, int y, int N, int dy);
  double T3v(int x1, int x2, int y, int N);
  double T4v(int x1, int x2, int y, int N);

  double T11v(int x1, int x2, int y, int N, int dx1);
  double T12v(int x1, int x2, int y, int N, int dx2);

  // Perform value iteration
  double vi(double thresh, double maxIter);

  // Print out decisions
  void decide(std::string filename);
};

#endif
