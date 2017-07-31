import std.stdio : writeln;
import std.math : abs;
import std.algorithm.comparison : min,max;
import std.parallelism;
import std.datetime;
import std.range : iota;

enum double lamda1 = 0.5;
enum double lamda2 = 0.5;
enum double mu1 = 1;
enum double mu2 = 1;
enum double b1 = 80;
enum double b2 = 100;
enum double h1 = 3;
enum double h2 = 2;
enum double alpha = 0.05;
enum int n = -1000;
enum int m = 1000;

__gshared float[(m-n+1)^^2] fvals = 0;

void main() {

    auto startTime = Clock.currTime(UTC());
    float error = vi(1e-5, cast(int)1e6);
    auto elapsedTime = Clock.currTime(UTC())- startTime;

    writeln("v(0,0): ", f(0,0), "\t\t-\tTotal Time: ", elapsedTime);
}

ref float f(int x1, int x2) {
    x1 = x1 - n;
    x2 = x2 - n;
    int i = x1*(m-n+1) + x2;
    return fvals[i];
}


int x1(int index) {
    return index / (m-n+1) + n;
}

int x2(int index) {
    return index % (m-n+1) + n;
}

float h(int x1, int x2) {
    float hVal = 0.0;
    hVal += (x1 < 0) ? -x1*b1 : x1*h1;
    hVal += (x2 < 0) ? -x2*b2 : x2*h2;
    return hVal;
}

float vi(float thresh, uint maxIter) {
    float error = thresh + 1;
    int i = 0;
    while (error > thresh && i++ < maxIter) {

        auto startTime = Clock.currTime(UTC());
        error = taskPool.reduce!max(taskPool.amap!viTask(iota(fvals.length)));
        auto elapsedTime = Clock.currTime(UTC())- startTime;

        writeln("Error: ", error, "\t\t-\tTime: ", elapsedTime);
    }
    return error;
}

float viTask(ulong index) {
    int x1 = x1(cast(int)index);
    int x2 = x2(cast(int)index);

    float old_f = f(x1,x2);

    if (x1 == m) {
        f(x1,x2) = f(x1-1, x2) + h1/alpha;
    } else if (x2 == m) {
        f(x1,x2) = f(x1, x2-1) + h2/alpha;
    } else if (x1 == n) {
        f(x1,x2) = f(x1+1, x2) + b1/alpha;
    } else if (x2 == n) {
        f(x1,x2) = f(x1, x2+1) + b2/alpha;
    } else {
        float f1 = f(x1-1, x2);
        float f2 = f(x1, x2-1);
        float f3 = f(x1+1, x2);
        float f4 = f(x1, x2+1);
        f(x1,x2) = (h(x1, x2) + lamda1*f1 + lamda2*f2 + mu1*min(f(x1,x2), f3, f4)) / (mu1+lamda1+lamda2+alpha);
    }

    return abs(f(x1,x2) - old_f);
}

