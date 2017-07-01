#include <iostream>
#include <algorithm>
#include <cmath>
#include <string>
#include <fstream>

#include "space.hpp"

Space::Space(unsigned int m) {
  this->m_ = m;

  // Initialize all states
  for (int i = 0; i <= m; i++) {
    for (int j = 0; j <= m; j++) {
      this->states_.push_back(State(i,j,0));
    }
  }
}

State& Space::state(unsigned int x, unsigned int y) {
  int index = x*(m_ + 1) + y;
  return this->states_[index];
}

double Space::vi(double thresh, double maxIter) {
  double error = thresh + 1;
  double i = 0;
  while (error > thresh && i++ < maxIter) {
    double localError = 0;
    for (State& s: this->states_) {
      unsigned int x = s.x;
      unsigned int y = s.y;
      double old_f = s.f;
      
      double T1v, T2v, T3v, T4v;

      if (x == 0) {
        T1v = this->state(x, y).f + cx;
      } else {
        T1v = this->state(x-1, y).f;
      }
      
      if (y == 0) {
        T2v = this->state(x, y).f + cy;

        if (x == this->m_) {
          T4v = std::min(s.f, this->state(x,y).f);
        } else {
          T4v = std::min(s.f, this->state(x+1,y).f);
        }

      } else {
        T2v = this->state(x, y-1).f;

        if (x == this->m_) {
          T4v = std::min(s.f, this->state(x,y-1).f);
        } else {
          T4v = std::min(s.f, this->state(x+1,y).f);
          T4v = std::min(T4v, this->state(x+1,y-1).f);
        }
      }

      if (y == this->m_) {
        T3v = this->state(x,y).f;
      } else {
        T3v = this->state(x, y+1).f;
      }

      s.f = (hx*x + hy*y + lambdax*T1v + lambday*T2v + delta*T3v + mu*T4v) / (alpha + lambdax + lambday
                                                                              +  mu + delta);

      localError = std::max(localError, std::fabs(s.f - old_f));
    }    
    error = localError;
  }
  return error;
}

enum {
  NoProduce = 0, 
  ProduceRaw = 1,
  ProduceUsed = 2
};

void Space::decide(std::string filename) {
  std::ofstream file(filename);

  file << "x\t" << "y\t" << "f-value\t" << "decision" << std::endl;

  for (State& s: this->states_) {
    unsigned int x = s.x;
    unsigned int y = s.y;

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
}
