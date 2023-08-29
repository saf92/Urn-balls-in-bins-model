#include <iostream>
#include <vector>

#include <ctime>
#include <cmath>
#include <cstdlib>

#include <fstream>
#include <iterator>

#include <random>

#include "urn_binary_search_functions.h"

using namespace std;

/*
Runs code for pure birth process for balls in bins process with feedback/generalised Polya urn model
Saves values from simulations in a .txt file that can be read into python
*/

//  To run:
// g++ -O9 -std=c++11 pure_birth_process.cpp urn_binary_search_functions.cpp -o main.out && ./main.out

// Mersenne twister random numbers
// https://stackoverflow.com/questions/22923551/generating-number-0-1-using-mersenne-twister-c
mt19937 generator (static_cast <unsigned> (time(0)));

//Single run of pure birth process with rates a*w^g corresponding to urn model (balls in bins model).

int pure_birth_process_single(int w, int w_add, double a, double g,
  double jt_max, int w_max){
    double jump_time = 0;
    double holding_time = 0;
    double f = 0;
    f = a*pow(w,g);
    exponential_distribution<double> dis(f);
    holding_time = dis(generator);
    jump_time = jump_time+holding_time;
    while (jump_time<jt_max && w<w_max){
      w=w+w_add;
      f = a*pow(w,g);
      exponential_distribution<double> dis(f);
      holding_time = dis(generator);
      jump_time = jump_time+holding_time;
    }
    return w;
  }
/*
Multiple runs of the pure birth process
with rates a*w^g corresponding to urn model (balls in bins model).
*/
  vector<int>pure_birth_process_agg(int runs, int w0, int w_add, double a, double g,
    double jt_max, int w_max){
      int w;
      vector<int>w_vec(runs);
      for (int i=0; i< runs; i++){
      w=w0;
      w = pure_birth_process_single(w, w_add, a, g, jt_max, w_max);
      w_vec[i]=w;
      }
      return w_vec;
    }

//Single run of pure birth process with rates a*w^g corresponding to urn model.
// Outputs jump time
double pure_birth_process_single_jt(int w, int w_add, double a, double g,
    int w_max){
      double jump_time = 0;
      double holding_time = 0;
      double f = 0;
      f = a*pow(w,g);
      exponential_distribution<double> dis(f);
      holding_time = dis(generator);
      jump_time = jump_time+holding_time;
      while (w<w_max){
        w=w+w_add;
        f = a*pow(w,g);
        exponential_distribution<double> dis(f);
        holding_time = dis(generator);
        jump_time = jump_time+holding_time;
      }
      return jump_time;
    }

// Prints vector of integers
void vector_print_int(vector<int> x){
  for(int i=0; i<x.size(); i++){
    cout << x[i] << endl;
  }
}

// Outputs vector of integers to a .txt file
void output_int_vec_in_txt(string file_path, vector<int> x){
  ofstream output_file(file_path);
  ostream_iterator<int> output_iterator(output_file, "\n");
  copy(x.begin(), x.end(), output_iterator);
}

int main(){

int w0=1;  // initial wealth
int w = 1;
int w_add=1; //wealth packet
double a =1;
double t_bp=0;
double gamma = 2.0;
string gamma_s=to_string(gamma);
cout <<  gamma_s.substr(0, 3) << endl;
double t_g = 1/(gamma-1);
double t = 1;
vector<double> ts ={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};
int w_max = (int) pow(10,4);
int runs = (int) pow(10,4);
vector<int> w_vec(runs,w);
vector<double> exp_times(runs,t);

long int start_s=clock();

// Pure birth process over multiple runs, outputs vector of aggregate of runs
w_vec=pure_birth_process_agg(runs, w0, w_add, a, gamma, t_g, w_max);

// ************************************************************************
// Pure birth process simulated over multiple runs and times, each aggregate saved as .txt file
// for (int i=0; i<ts.size();i++){
// t=ts[i];
// w_vec=pure_birth_process_agg(runs, w0, w_add, a, gamma, t, w_max);
// cout << i << endl;
// string t_s=to_string(t);
// output_int_vec_in_txt("./output_bp/pbp_t"+t_s.substr(0, 3)+
// "g"+gamma_s.substr(0, 3)+".txt", w_vec);
// }
// *************************************************************************


long int stop_s=clock();

// ****************************************************************************
// Pure birth process simulated over multiple runs, output is jump time
// long int start_s=clock();
// for (int i=0; i< runs; i++){
// w=w0;
// t_bp = pure_birth_process_single_jt(w, w_add, a, gamma, w_max);
// exp_times[i]=t_bp;
// }
// long int stop_s=clock();
// *****************************************************************************

cout << "Time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)/60 << endl << endl;

// Save wealth values as .txt file
output_int_vec_in_txt("./output_bp/pbp_g"+gamma_s.substr(0, 3)+".txt", w_vec);

// Save jump times as .txt file
// output_vec_in_txt("./output_bp/pbp_exp_t"+gamma_s.substr(0, 3)+".txt", exp_times);

return 0;

}
