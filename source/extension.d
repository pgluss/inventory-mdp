import std.stdio;
import std.math : abs, isNaN;
import std.algorithm.comparison : min,max;
import std.parallelism;
import std.datetime;
import std.range : iota;

enum double lambdax1 = 0.5;
enum double lambdax2 = 0.5;
enum double lambday  = 0.5; 

enum double mu = 1;
enum double gamma = 0.5;
enum double alpha = 0.05;

enum double bx1 = 80;
enum double bx2 = 100;
enum double by  = 50;

enum double hx1 = 2;
enum double hx2 = 2;
enum double hy  = 1;

enum int Dx1 = 3;
enum int Dx2 = 3;
enum int Dy  = 3;

enum double pn = 0.3;
enum double po = 0.3;

enum int n = -20;       // Minimum allowed inventory for stock
enum int m = 20;        // Maximum allowed inventory for stock
enum int Nmin = 1;      // Minumum simulated population
enum int M = 51;        // Maximum simulated population

enum double EN = (lambdax1 + lambdax2 + lambday) / (gamma*(1 - pn)); // Expected value of N

__gshared float[] fvals;

// The main function, uses the below functions to perform value iteration
void main() {
    fvals = new float[(M-Nmin+1)*(m-n+1)^^3];

    for (int i = 0; i < fvals.length; i++) {
        fvals[i] = 0;
    }
    

    auto startTime = Clock.currTime(UTC());
    int numIter = vi(1e-2);
    auto elapsedTime = Clock.currTime(UTC())- startTime;
    
    decide("test-e.csv");

    stdout.writefln("Len: %d, \tNum Iter: %d, \tTotal Time: %s", fvals.length, numIter, elapsedTime);
}

// Takes in a value x and returns x is x is positive and 0 if x is negative
int plus(int x) {
    return ((x > 0) ? (x) : 0);
}

// Takes in a value x and returns -x if x is negative and 0 if x is positive
int minus(int x) {
    return ((x < 0) ? -(x) : 0);
}

// Indicator: returns 1 if x is negative and 0 if x is positive
int indic(int x) {
    return ((x < 0) ? 1 : 0);
}

// Returns by reference the f-value for a given state
// Handles boundary conditions for each input variable
ref float f(int x1, int x2, int y, int N) {
    if (x1 < n) {
        x1 = n;
    } else if (x1 > m) {
        x1 = m;
    }

    if (x2 < n) {
        x2 = n;
    } else if (x2 > m) {
        x2 = m;
    }

    if (y < n) {
        y = n;
    } else if (y > m) {
        y = m;
    }
    
    if (N < Nmin) {
        N = Nmin;
    } else if (N > M) {
        N = M;
    }

    x1 = x1 - n;
    x2 = x2 - n;
    y  = y  - n;
    N  = N  - Nmin;
    int i = x1*(M-Nmin+1)*(m-n+1)^^2 + x2*(M-Nmin+1)*(m-n+1) + y*(M-Nmin+1) + N;
    return fvals[i];
}

// Return the x1 value associated with the given index in the array
// Note that we are representing a 4D matrix as a one dimensional array
int x1(int index) {
    return index / ((M-Nmin+1)*(m-n+1)^^2) + n;
}

// Return the x2 value associated with the given index in the array
// Note that we are representing a 4D matrix as a one dimensional array
int x2(int index) {
    return (index % ((M-Nmin+1)*(m-n+1)^^2)) / ((M-Nmin+1)*(m-n+1)) + n;
}

// Return the y value associated with the given index in the array
// Note that we are representing a 4D matrix as a one dimensional array
int y(int index) {
    return (index % ((M-Nmin+1)*(m-n+1))) / (M-Nmin+1) + n;
}

// Return the N value associated with the given index in the array
// Note that we are representing a 4D matrix as a one dimensional array
int N(int index) {
    return index % (M-Nmin+1) + Nmin;
}

// Calculates the sum of the holding cost and backorder costs for a given state
// It depends on current stock levels
float h(int x1, int x2, int y) {
    float hVal = 0.0;
    hVal += (x1 < 0) ? -x1*bx1 : x1*hx1;
    hVal += (x2 < 0) ? -x2*bx2 : x2*hx2;
    hVal += (y  < 0) ? -y*by : y*hy;
    return hVal;
}

// Calculates the probability of bulk demand of size d
// Currently it's constant depending on maximum bulk demand size simulated
double px1(int d) {
    return 1.0/Dx1;
}

// Calculates the probability of bulk demand of size d
// Currently it's constant depending on maximum bulk demand size simulated
double px2(int d) {
    return 1.0/Dx2;
}

// Calculates the probability of bulk demand of size d
// Currently it's constant depending on maximum bulk demand size simulated
double py(int d) {
    return 1.0/Dy;
}

// Calculates T1v as per the equation found in the paper
// See section "Formulation of the Extended Model"
double T1v(int x1, int x2, int y, int N) {
  double T1v = 0;
  T1v += (lambdax1 / (lambdax1 + lambdax2))*T11v(x1, x2, y, N, 1);
  T1v += (lambdax2 / (lambdax1 + lambdax2))*T12v(x1, x2, y, N, 1);
  return T1v;
}

// Calculates T2v as per the equation found in the paper
// See section "Formulation of the Extended Model"
double T2v(int x1, int x2, int y, int N, int dy) {
  return f(x1, x2, y - dy, N);
}

// Calculates T3v as per the equation found in the paper
// See section "Formulation of the Extended Model"
double T3v(int x1, int x2, int y, int N) {
  double T3v = 0;
  T3v += (1 - pn - po)*f(x1, x2, y + 1, N - 1);
  T3v += pn*T1v(x1, x2, y + 1, N - 1);
  T3v += po*T2v(x1, x2, y + 1, N - 1, 1);
  T3v *= (N / M);
  T3v += ((M - N) / M)*f(x1, x2, y, N);
  return T3v;
}

// Calculates T4v as per the equation found in the paper
// See section "Formulation of the Extended Model"
double T4v(int x1, int x2, int y, int N) {
  double T4v = 0;
  if (y <= 0) {
    double f1 = f(x1, x2, y, N);
    double f2 = f(x1 + 1, x2, y, N + indic(x1));
    double f3 = f(x1, x2 + 1, y, N + indic(x2));

    T4v = min(f1, f2, f3);
  } else {
    double f1 = f(x1, x2, y, N);
    double f2 = f(x1 + 1, x2, y, N + indic(x1));
    double f3 = f(x1, x2 + 1, y, N + indic(x2));
    double f4 = f(x1 + 1, x2, y - 1, N + indic(x1));
    double f5 = f(x1, x2 + 1, y - 1, N + indic(x2));

    T4v = min(f1, f2, f3, f4, f5);
  }
  return T4v;
}

// Calculates T11v as per the equation found in the paper
// See section "Formulation of the Extended Model"
double T11v(int x1, int x2, int y, int N, int dx1) {
  double T11v = 0;
  if (dx1 <= x1) {
    T11v = f(x1 - dx1, x2, y, N + dx1);

  } else if (x1 < dx1 && dx1 <= (x1 + x2)) {
    double f1 = f(x1 - dx1, x2, y, N + plus(x1));
    double f2 = f(0, x2 - (dx1 - x1), y, N + dx1);
    T11v = min(f1, f2);

  } else if (dx1 > (x1 + x2)) {
    double f1 = f(x1 - dx1, x2, y, N + plus(x1));
    double f2 = f(x1 + x2 - dx1, 0, y, N + plus(dx1 - x1 - x2));
    T11v = min(f1, f2);

  }
  return T11v;
}

// Calculates T12v as per the equation found in the paper
// See section "Formulation of the Extended Model"
double T12v(int x1, int x2, int y, int N, int dx2) {
  double T12v = 0;
  if (dx2 <= x2) {
    T12v = f(x1, x2 - dx2, y, N + dx2);

  } else if (x2 < dx2 && dx2 <= (x1 + x2)) {
    double f1 = f(x1, x2 - dx2, y, N + plus(x2));
    double f2 = f(x1 - (dx2 - x2), 0, y, N + dx2);
    T12v = min(f1, f2);

  } else if (dx2 > (x1 + x2)) {
    double f1 = f(x1, x2 - dx2, y, N + plus(x2));
    double f2 = f(0, x1 + x2 - dx2, y, N + plus(dx2 - x1 - x2));
    T12v = min(f1, f2);

  }
  return T12v;
}

// Performs value iteration until the error is below the specified threshold
// or the maximum number of iterations have occurred
int vi(float thresh, uint maxIter = cast(int)1e6) {
    float maxError = thresh + 1;
    float avgError = 0;
    int i = 0;

    while (maxError > thresh && i++ < maxIter) {
      
        auto startTime = Clock.currTime(UTC());

        // Sequential
        //float localError = 0;
        //for (int index = 0; index < fvals.length; index++) {
        //  localError = max(viTask(index), localError);
        //}
        //maxError = localError;


        // Parallel: Uses multiple CPUs
        auto localErrors = taskPool.amap!viTask(iota(cast(int)fvals.length));
        maxError = taskPool.reduce!max(localErrors);
        avgError = taskPool.reduce!"a+b"(localErrors)/localErrors.length;


        auto elapsedTime = Clock.currTime(UTC()) - startTime;

        stdout.writefln("%d \t\t %15f %15f \t\t %s", i, avgError, maxError, elapsedTime);
    }
    return i;
}

// Calculates the new f-value for a given index
float viTask(int index) {
    int x1 = x1(index);
    int x2 = x2(index);
    int y  =  y(index);
    int N  =  N(index);

    float old_f = f(x1,x2, y, N);

    float fx1 = 0;
    float fx2 = 0;
    float fy  = 0;

    for (int dx1 = 1; dx1 <= Dx1; dx1++) {
        fx1 += px1(dx1)*T11v(x1, x2, y, N, dx1);
    }

    for (int dx2 = 1; dx2 <= Dx2; dx2++) {
        fx2 += px2(dx2)*T12v(x1, x2, y, N, dx2);
    }

    for (int dy = 1; dy <= Dy; dy++) {
        fy += py(dy)*T12v(x1, x2, y, N, dy);
    }

    float f1 = M*gamma*T3v(x1, x2, y, N);
    float f2 = mu*T4v(x1, x2, y, N);

    f(x1, x2, y, N) = (h(x1, x2, y) + lambdax1*fx1 + lambdax2*fx2 + lambday*fy + f1 + f2) /
                                      (alpha + lambdax1 + lambdax2 + lambday + mu + M*gamma);

    return abs(f(x1, x2, y, N) - old_f);
}

enum int NoProduce = 0;     // Idle
enum int ProduceRawX1 = 1;  // Produce new product x1 from raw materials
enum int ProduceRawX2 = 2;  // Produce new product x2 from raw materials
enum int ProduceUsedX1 = 3; // Produce new product x1 from used materials
enum int ProduceUsedX2 = 4; // Produce new product x2 from used materials

// Get decisions for each state and print to a file
void decide(string filename) {
    File file = File(filename, "w");

    for (int i = 0; i < fvals.length; i++) {
        int x1 = x1(i);
        int x2 = x2(i);
        int y  =  y(i);
        int N  =  N(i);

        double fval = f(x1, x2, y, N);
        int d;

        if (x1 != m && x2 != m && y > 0) {
            if (f(x1,x2,y,N) < f(x1+1,x2,y,N) && 
                f(x1,x2,y,N) < f(x1,x2+1,y,N) &&
                f(x1,x2,y,N) < f(x1+1,x2,y-1,N) &&
                f(x1,x2,y,N) < f(x1,x2+1,y-1,N)) {
            
                d = NoProduce;
            } else if (f(x1+1,x2,y,N) < f(x1,x2,y,N) &&
                       f(x1+1,x2,y,N) < f(x1,x2+1,y,N) &&
                       f(x1+1,x2,y,N) < f(x1+1,x2,y-1,N) &&
                       f(x1+1,x2,y,N) < f(x1,x2+1,y-1,N)) {
                
                d = ProduceRawX1;
            } else if (f(x1,x2+1,y,N) < f(x1,x2,y,N) &&
                       f(x1,x2+1,y,N) < f(x1+1,x2,y,N) &&
                       f(x1,x2+1,y,N) < f(x1+1,x2,y-1,N) &&
                       f(x1,x2+1,y,N) < f(x1,x2+1,y-1,N)) {
            
                d = ProduceRawX2;
            } else if (f(x1+1,x2,y-1,N) < f(x1,x2,y,N) &&
                       f(x1+1,x2,y-1,N) < f(x1,x2+1,y,N) &&
                       f(x1+1,x2,y-1,N) < f(x1+1,x2,y,N) &&
                       f(x1+1,x2,y-1,N) < f(x1,x2+1,y-1,N)) {

                d = ProduceUsedX1;
            } else if (f(x1,x2+1,y-1,N) < f(x1,x2,y,N) &&
                       f(x1,x2+1,y-1,N) < f(x1+1,x2,y,N) &&
                       f(x1,x2+1,y-1,N) < f(x1+1,x2,y-1,N) &&
                       f(x1,x2+1,y-1,N) < f(x1,x2+1,y,N)) {
                
                d = ProduceUsedX2;
            }
        } else if (x1 != m && x2 != m && y <= 0) {
            if (f(x1,x2,y,N) < f(x1+1,x2,y,N) && 
                f(x1,x2,y,N) < f(x1,x2+1,y,N)) {
            
                d = NoProduce;
            } else if (f(x1+1,x2,y,N) < f(x1,x2,y,N) &&
                       f(x1+1,x2,y,N) < f(x1,x2+1,y,N)) {
                
                d = ProduceRawX1;
            } else if (f(x1,x2+1,y,N) < f(x1,x2,y,N) &&
                       f(x1,x2+1,y,N) < f(x1+1,x2,y,N)) {
            
                d = ProduceRawX2;
            }
        } else if (x1 == m && x2 != m && y > 0){
            if (f(x1,x2,y,N) < f(x1,x2+1,y,N) &&
                f(x1,x2,y,N) < f(x1,x2+1,y-1,N)) {
            
                d = NoProduce;
            } else if (f(x1,x2+1,y,N) < f(x1,x2,y,N) &&
                       f(x1,x2+1,y,N) < f(x1,x2+1,y-1,N)) {
            
                d = ProduceRawX2;
            } else if (f(x1,x2+1,y-1,N) < f(x1,x2,y,N) &&
                       f(x1,x2+1,y-1,N) < f(x1,x2+1,y,N)) {
                
                d = ProduceUsedX2;
            }
        } else if (x1 == m && x2 != m && y <= 0) {
            if (f(x1,x2,y,N) < f(x1,x2+1,y,N)) {
            
                d = NoProduce;
            } else if (f(x1,x2+1,y,N) < f(x1,x2,y,N)) {
            
                d = ProduceRawX2;
            }
        } else if (x2 != m && x2 == m && y > 0) {
            if (f(x1,x2,y,N) < f(x1+1,x2,y,N) && 
                f(x1,x2,y,N) < f(x1+1,x2,y-1,N)) {
            
                d = NoProduce;
            } else if (f(x1+1,x2,y,N) < f(x1,x2,y,N) &&
                       f(x1+1,x2,y,N) < f(x1+1,x2,y-1,N)) {
                
                d = ProduceRawX1;
            } else if (f(x1+1,x2,y-1,N) < f(x1,x2,y,N) &&
                       f(x1+1,x2,y-1,N) < f(x1+1,x2,y,N)) {

                d = ProduceUsedX1;
            }
        } else if (x2 != m && x2 == m && y <= 0) {
            if (f(x1,x2,y,N) < f(x1+1,x2,y,N)) { 
            
                d = NoProduce;
            } else if (f(x1+1,x2,y,N) < f(x1,x2,y,N)) {
                
                d = ProduceRawX1;
            }
        } else {
            d = NoProduce;
        }

        file.writefln("%d, %d, %d, %d, %f, %d", x1, x2, y, N, fval, d);
    }
}
