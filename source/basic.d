import std.stdio : writeln, File;
import std.math : abs;
import std.algorithm.comparison : min,max;
import std.parallelism;
import std.datetime;
import std.range : iota;

enum double lambdax = 0.85;
enum double lambday = 0.85;
enum double mu = 1;
enum double delta = 0.8;
enum double cx = 11.4;
enum double cy = 5.7;
enum double hx = 2;
enum double hy = 1;
enum double alpha = 0.05;
enum int m = 100;

__gshared float[(m+1)^^2] fvals = 0;

// The main function. Runs the value iteration function.
void main() {

    auto startTime = Clock.currTime(UTC());
    float error = vi(1e-10, cast(int)1e10);
    auto elapsedTime = Clock.currTime(UTC())- startTime;

    writeln("v(0,0): ", f(0,0), "\t\t-\tTotal Time: ", elapsedTime);

    decide("test1.csv");
}

// Returns the f-value for a given (x,y) state
ref float f(int x, int y) {
    int i = x*(m+1) + y;
    return fvals[i];
}

// Returns the x value for a given index in the array of f-values
int x(int index) {
    return index / (m+1);
}

// Returns the y value for a givne index in the array of f-values
int y(int index) {
    return index % (m+1);
}

// Performs value iterations until the error is below the specified threshold
// or the max number of iterations have occurred
float vi(float thresh, uint maxIter) {
    float error = thresh + 1;
    int i = 0;
    while (error > thresh && i++ < maxIter) {

        auto startTime = Clock.currTime(UTC());
        error = taskPool.reduce!max(taskPool.amap!viTask(iota(cast(int)fvals.length)));
        auto elapsedTime = Clock.currTime(UTC())- startTime;

        writeln("Error: ", error, "\t\t-\tTime: ", elapsedTime);
    }
    return error;
}

// Does value iteration for a single state
float viTask(int index) {
    int x = x(index);
    int y = y(index);

    float old_f = f(x,y);

    double T1v, T2v, T3v, T4v;

    if (x == 0) {
        T1v = f(x, y) + cx;
    } else {
        T1v = f(x-1, y);
    }

    if (y == 0) {
        T2v = f(x, y) + cy;

        if (x == m) {
            T4v = min(f(x,y), f(x,y));
        } else {
            T4v = min(f(x,y), f(x+1,y));
        }

    } else {
        T2v = f(x, y-1);

        if (x == m) {
            T4v = min(f(x,y), f(x,y-1));
        } else {
            T4v = min(f(x,y), f(x+1,y), f(x+1,y-1));
        }
    }

    if (y == m) {
        T3v = f(x,y);
    } else {
        T3v = f(x, y+1);
    }

    f(x,y) = (hx*x + hy*y + lambdax*T1v + lambday*T2v + delta*T3v + mu*T4v) / (alpha + lambdax + lambday
            +  mu + delta);

    return abs(f(x,y) - old_f);
}

enum int NoProduce = 0;   // Idle
enum int ProduceRaw = 1;  // Produce new product from raw materials
enum int ProduceUsed = 2; // Produce new product by refurbishing used product

// Get decisions for each state and print to a file
void decide(string filename) {
    File file = File(filename, "w");
    file.writefln("x, y, f-value, decision");
  
    for (int i = 0; i < fvals.length; i++) {
        int x = x(i);
        int y = y(i);

        double fval = f(x,y);
        int d;

        if (x == m) {
            d = NoProduce;
        } else if (y == 0) {
            if (f(x,y) < f(x+1,y)) {
                d = NoProduce;
            } else {
                d = ProduceRaw;
            }
        } else if (y > 0) {
            if (f(x,y) < f(x+1,y) && f(x,y) < f(x+1,y-1)) {
                d = NoProduce;
            } else if (f(x+1,y) < f(x,y) && f(x+1,y) < f(x+1,y-1)) {
                d = ProduceRaw;
            } else {
                d = ProduceUsed;
            }
        } else {
            assert(false, "Error: Check x and y!");
        }

        file.writefln("%d, %d, %f, %d", x, y, fval, d);
    }
}
