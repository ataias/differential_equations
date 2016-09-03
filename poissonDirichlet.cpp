#include <iostream>
#include <chrono>
#include <functional>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;

//compile: g++ -std=c++14 poissonDirichlet.cpp -o solver_dirichlet

// Thanks for an answer here, I am using this timer
//http://stackoverflow.com/questions/728068/how-to-calculate-a-time-difference-in-c
class Timer
{
public:
    Timer() : beg_(clock_::now()) {}
    void reset() { beg_ = clock_::now(); }
    double elapsed() const {
        return std::chrono::duration_cast<second_>
            (clock_::now() - beg_).count(); }

private:
    typedef std::chrono::high_resolution_clock clock_;
    typedef std::chrono::duration<double, std::ratio<1> > second_;
    std::chrono::time_point<clock_> beg_;
};


double abs(double x) {
  if(x < 0) return -x;
  return x;
}

// Remember we are considering number of points in x == y
double ** getMesh(int n) {
  // we use calloc to initialize everything with zeros
  double ** mesh = (double **) calloc (n, sizeof(double *));
  for (int i = 0; i < n; i++) {
    mesh[i] = (double *) calloc (n, sizeof(double));
  }

  return mesh;
}

/**
  g is a function of (x,y) that is applied on the boundaries
*/
void computeBoundaries(double ** mesh, int n, function<double(double, double)> g) {
  int i, j;
  double dx = 1 / (n - 1);
  i = 0;
  for (j = 0; j < n; j++) mesh[i][j] = g(i * dx, j * dx);
  i = n - 1;
  for (j = 0; j < n; j++) mesh[i][j] = g(i * dx, j * dx);
  j = 0;
  for (i = 0; i < n; i++) mesh[i][j] = g(i * dx, j * dx);
  j = n - 1;
  for (i = 0; i < n; i++) mesh[i][j] = g(i * dx, j * dx);

  return;
}

/**
 f is a function that expects mesh and i, j and computes the value of
 that point. It also takes the absolute change reference.
 The return of this function is the absolute change
*/
double computeInternalPoints(double ** mesh, int n, function<double(double**, int, int, double*)> f) {
  double absoluteChange = 0.0;
  // only internal points, that's why there is the loop goes from 0 to n - 2
  for (int i = 1; i < n - 1; i++) {
    for (int j = 1; j < n - 1; j++) {
      mesh[i][j] = f(mesh, i, j, &absoluteChange);
    }
  }

  return absoluteChange;
}

void getSolutionDirichlet(double **mesh, int n, function<double(double, double)> g, double epsilon) {

  double absoluteChange = 0.0;
  double dx = 1.0 / (n - 1);

  // nabla f = u
  auto u = [] (double x, double y)  {return x*x + y*y;};

  // f_ij function
  auto f = [u, dx](auto mesh, int i, int j, auto absoluteChange) {
    double old_value = mesh[i][j];
    double c = mesh[i-1][j] + mesh[i+1][j] + mesh[i][j-1] + mesh[i][j+1];
    mesh[i][j] = (c - u(i * dx, j * dx)*dx*dx) / 4.0;
    *absoluteChange += abs(old_value - mesh[i][j]);
    return mesh[i][j];
  };

  // we just need to compute the boundaries once
  computeBoundaries(mesh, n, g);
  int k = 0;
  Timer tmr;
  do {
    absoluteChange = computeInternalPoints(mesh, n, f);
    k++;
  } while (absoluteChange > epsilon);

  cout << "Finished computing! Latest absolute change was: " << absoluteChange << endl;
  cout << "Total number of iterations on mesh: " << k << endl;
  double time_elapsed = tmr.elapsed();
  cout << "Total time spent: " << time_elapsed << " seconds" << endl;
  cout << "Total time per iteration: " << time_elapsed / k << " seconds" << endl;
  return;
}

void writeDataToFile(string filename, double **mesh, int n){
  ofstream file(filename, ios::out|ios::binary);
  file.write((char *) &n, sizeof(int));
  for (int i = 0; i < n; i++) file.write((char *) mesh[i], n*sizeof(double));
  cout << "File size: " << file.tellp() << " bytes" << endl;
  file.close();
  cout << "File successfully saved at: " << filename << endl;
}

void freeArray(double **mesh, int n) {
  for (int i = 0; i < n; i++) free(mesh[i]);
  free(mesh);
}


int main(int argc, char *argv[]) {
  if (argc < 3) {
    cout << "Usage: solver N outputFilename" << endl;
    return 1;
  }

  int n = 200;
  istringstream ss(argv[1]);
  if (!(ss >> n)){
    cerr << "Invalid number: " << argv[1] << '\n';
    return 1;
  }

  string filename(argv[2]);
  // Boundary condition
  auto g = [](double x, double y) {
    if (abs(x) < 0.0001) return (y - 0.5)*(y - 0.5);
    else if (abs(x) > 0.999) return (y - 0.5)*(y - 0.5);
    else if (abs(y) < 0.0001) return (x - 0.5)*(x-0.5);
    else if (abs(y) > 0.999) return (x - 0.5)*(x-0.5);
    return 0.0;
  };


  double epsilon = 0.0001;
  double ** mesh = getMesh(n);
  getSolutionDirichlet(mesh, n, g, epsilon);
  writeDataToFile(filename, mesh, n);
  // freeing memory is important
  freeArray(mesh, n);

  return 0;
}
