#include <iostream>
#include <algorithm>
#include <cmath>
#include <string>
#include <fstream>

#include "space.hpp"

#define plus(x)  ((x > 0) ?  (x) : 0)
#define minus(x) ((x < 0) ? -(x) : 0)
#define indic(x) ((x < 0) ?   1  : 0)

Space::Space(int n, int m, int M) {
  this->m_ = m;
  this->n_ = n;
  this->M_ = M;

  // Initialize all states
  for (int l = 0; l <= M; l++) {
    for (int i = n; i <= m; i++) {
      for (int j = n; j <= m; j++) {
        for (int k = n; k <= m; k++) {
          this->states_.push_back(State(l,i,j,k,0));
        }
      }
    }
  }
}

State& Space::state(int N, int x1, int x2, int y) {
  
  // Check boundary conditions and adjust as needed
  if (x1 < this->n_) {
    x1 = this->n_;
  } else if (x1 > this->m_) {
    x1 = this->m_;
  }

  if (x2 < this->n_) {
    x2 = this->n_;
  } else if (x2 > this->m_) {
    x2 = this->m_;
  }
  
  if (y < this->n_) {
    y = this->n_;
  } else if (y > this->m_) {
    y = this->m_;
  }

  if (N < 0) {
    N = 0;  
  } else if (N > this->M_) {
    N = this->M_;
  }

  // Calculate index
  x1 = x1 - this->n_;
  x2 = x2 - this->n_;
  y  =  y - this->n_;
  int index = N*(m_ + 1)*(m_ + 1)*(m_ + 1) + x1*(m_ - n_ + 1)*(m_ - n_ + 1) + x2*(m_ - n_ + 1) + y;
  return this->states_[index];
}

double Space::h(int x1, int x2, int y) {
  double hVal = 0.0;
  hVal += (x1 < 0) ? -x1*bx1 : x1*hx1;
  hVal += (x2 < 0) ? -x2*bx2 : x2*hx2;
  hVal += (y  < 0) ? -y*by : y*hy;
  return hVal;
}

double Space::px1(int d) {
  return (1.0 / this->Dx1);
}

double Space::px2(int d) {
  return (1.0 / this->Dx2);
}

double Space::py(int d) {
  return (1.0 / this->Dy);
}

double Space::T1v(int x1, int x2, int y, int N) {
  double T1v = 0;
  T1v += (this->lambdax1 / (this->lambdax1 + this->lambdax2))*T11v(x1, x2, y, N, 1);
  T1v += (this->lambdax2 / (this->lambdax1 + this->lambdax2))*T12v(x1, x2, y, N, 1);
  return T1v;
}

double Space::T2v(int x1, int x2, int y, int N, int dy) {
  return this->state(x1, x2, y - dy, N).f;
}

double Space::T3v(int x1, int x2, int y, int N) {
  double T3v = 0;
  T3v += (1 - this->pn - this->po)*this->state(x1, x2, y + 1, N - 1).f;
  T3v += this->pn*T1v(x1, x2, y + 1, N - 1);
  T3v += this->po*T2v(x1, x2, y + 1, N - 1, 1);
  T3v *= (N / this->M_);
  T3v += ((this->M_ - N) / this->M_)*this->state(x1, x2, y, N).f;
  return T3v;
}

double Space::T4v(int x1, int x2, int y, int N) {
  double T4v = 0;
  if (y <= 0) {
    double f1 = this->state(x1, x2, y, N).f;
    double f2 = this->state(x1 + 1, x2, y, N + indic(x1)).f;
    double f3 = this->state(x1, x2 + 1, y, N + indic(x2)).f;

    T4v = fmin(f1, f2);
    T4v = fmin(T4v, f3);
  } else {
    double f1 = this->state(x1, x2, y, N).f;
    double f2 = this->state(x1 + 1, x2, y, N + indic(x1)).f;
    double f3 = this->state(x1, x2 + 1, y, N + indic(x2)).f;
    double f4 = this->state(x1 + 1, x2, y - 1, N + indic(x1)).f;
    double f5 = this->state(x1, x2 + 1, y - 1, N + indic(x2)).f;

    T4v = fmin(f1, f2);
    T4v = fmin(T4v, f3);
    T4v = fmin(T4v, f4);
    T4v = fmin(T4v, f5);
  }
  return T4v;
}

double Space::T11v(int x1, int x2, int y, int N, int dx1) {
  double T11v = 0;
  if (dx1 <= x1) {
    T11v = this->state(x1 - dx1, x2, y, N + dx1).f;

  } else if (x1 < dx1 && dx1 <= (x1 + x2)) {
    double f1 = this->state(x1 - dx1, x2, y, N + plus(x1)).f;
    double f2 = this->state(0, x2 - (dx1 - x1), y, N + dx1).f;
    T11v = fmin(f1, f2);

  } else if (dx1 > (x1 + x2)) {
    double f1 = this->state(x1 - dx1, x2, y, N + plus(x1)).f;
    double f2 = this->state(x1 + x2 - dx1, 0, y, N + plus(dx1 - x1 - x2)).f;
    T11v = fmin(f1, f2);

  }
  return T11v;
}

double Space::T12v(int x1, int x2, int y, int N, int dx2) {
  double T12v = 0;
  if (dx2 <= x2) {
    T12v = this->state(x1, x2 - dx2, y, N + dx2).f;

  } else if (x2 < dx2 && dx2 <= (x1 + x2)) {
    double f1 = this->state(x1, x2 - dx2, y, N + plus(x2)).f;
    double f2 = this->state(x1 - (dx2 - x2), 0, y, N + dx2).f;
    T12v = fmin(f1, f2);

  } else if (dx2 > (x1 + x2)) {
    double f1 = this->state(x1, x2 - dx2, y, N + plus(x2)).f;
    double f2 = this->state(0, x1 + x2 - dx2, y, N + plus(dx2 - x1 - x2)).f;
    T12v = fmin(f1, f2);

  }
  return T12v;
}

double Space::vi(double thresh, double maxIter) {
  double error = thresh + 1;
  double i = 0;
  while (error > thresh && i++ < maxIter) {
    std::cout << i << "\t" << error << std::endl;
    double localError = 0;
    for (State& s: this->states_) {

      int x1 = s.x1;
      int x2 = s.x2;
      int y  = s.y;
      unsigned int N = s.N;
      double old_f = s.f;
      
      double fx1 = 0;
      double fx2 = 0;
      double fy  = 0;

      for (int dx1 = 1; dx1 <= this->Dx1; dx1++) {
        fx1 += px1(dx1);
        fx1 *= T11v(x1, x2, y, N, dx1);
      }

      for (int dx2 = 1; dx2 <= this->Dx2; dx2++) {
        fx2 += px2(i)*T12v(x1, x2, y, N, dx2);
      }

      for (int dy = 1; dy <= this->Dy; dy++) {
        fy += py(dy)*T2v(x1, x2, y, N, dy);
      }

      double f1 = this->M_*gamma*T3v(x1, x2, y, N);
      double f2 = this->mu*T4v(x1, x2, y, N);

      s.f = (h(x1,x2,y) + lambdax1*fx1 + lambdax2*fx2 + lambday*fy + f1 + f2) / 
                                    (alpha + lambdax1 + lambdax2 + lambday + mu + this->M_*gamma);

      localError = std::max(localError, std::fabs(s.f - old_f));
      
    }    
    error = localError;
  }
  return error;
}

/*enum {
  NoProduce = 0, 
  ProduceRaw = 1,
  ProduceUsed = 2
};

void Space::decide(std::string filename) {
  std::ofstream file(filename);

  file << "x\t" << "y\t" << "f-value\t" << "decision" << std::endl;

  for (State& s: this->states_) {
    int x = s.x;
    int y = s.y;

    file << x << "\t" << y << "\t" << s.f << "\t";

    if (x == m_) {
      file << NoProduce << std::endl;
    } else if (y == 0) {
      if (s.f < this->state(x+1,y).f) {
        file << NoProduce << std::endl;
      } else {
        file << ProduceRaw << std::endl;
      }
    } else if (y > 0) {
      if (s.f < this->state(x+1,y).f && s.f < this->state(x+1,y-1).f) {
        file << NoProduce << std::endl;
      } else if (this->state(x+1,y).f < s.f && this->state(x+1,y).f < this->state(x+1,y-1).f) {
        file << ProduceRaw << std::endl;
      } else {
        file << ProduceUsed << std::endl;
      }
    } else {
      std::cerr << "Error: Check x and y!" << std::endl;
    }
  }

  file.close();
}*/
