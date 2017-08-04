import std.stdio;
import std.math : abs;
import std.algorithm.comparison : min,max;
import std.parallelism;
import std.datetime;
import std.range : iota;
import std.string;
import std.conv;
import std.array;

double lambdax, lambday, mu, delta, alpha;
double cx, cy, hx, hy;
int m;

__gshared float[] fvals;

//__gshared float[(m+1)^^2] fvals = 0;

// The main function. Runs the value iteration function.
void main() {
    File params = File("params-basic.csv");
    int batchID = 0;
    foreach (line; params.byLine) {
        if (!line.empty && line[0] != '#') {
            batchID++;
            char[][] tokens = line.split(",");

            lambdax = tokens[0].strip.to!double;
            lambday = tokens[1].strip.to!double;
            mu      = tokens[2].strip.to!double;
            delta   = tokens[3].strip.to!double;
            cx      = tokens[4].strip.to!double;
            cy      = tokens[5].strip.to!double;
            hx      = tokens[6].strip.to!double;
            hy      = tokens[7].strip.to!double;
            alpha   = tokens[8].strip.to!double;
            m       = tokens[9].strip.to!int;

            fvals = new float[(m+1)^^2];

            for (int i = 0; i < fvals.length; i++) {
                fvals[i] = 0;
            }

            auto startTime = Clock.currTime(UTC());
            float error = vi(1e-10, cast(int)1e10);
            auto elapsedTime = Clock.currTime(UTC())- startTime;

            decide("basic-batch-" ~ batchID.to!string ~ ".csv");

            stdout.writefln("Batch %d is complete.", batchID);
            stdout.writefln("v(0,0): %f \t Total Time: %s \n", f(0,0), elapsedTime);
        }
    }
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
