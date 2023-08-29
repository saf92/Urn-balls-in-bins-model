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

// ********************************************************************

// Functions to print out double or long int vector

void vector_print_double(vector<double> x){
  for(int i=0; i<x.size(); i++){
    cout << x[i] << endl;
  }
}

void vector_print_long_int(vector<long int> x){
  for(int i=0; i<x.size(); i++){
    cout << x[i] << endl;
  }
}

// **********************************************************************

//Gets the sum of a double vector

double sum(vector<double>x){
  double t=0;
  for(int i=0; i< x.size();i++){
    t+=x[i];
  }
  return t;
}

//Gets a vector of the cumulative sum of a double vector

vector<double> cumsum(vector<double>x){
  vector<double> Cumsum;
  double c=0;
  for(int i=0; i<x.size();i++){
    c+=x[i];
    Cumsum.push_back(c);
  }

  return Cumsum;
}

// ******************************************************************************

/*

The normalised rates are the categorical probabilities. Here we sample from the
categorical distribution with the normalised rates and return the rate this corresponds to.
*/

double get_rate(vector<double> rates, double u){
  int i=0;
  vector<double> cum_sum = cumsum(rates);
  int l = rates.size()-1;
  if(u >= cum_sum[l]){
    i = l;
  }
  else{
   while(cum_sum[i] <= u){
     i++;
   }
 }
   return rates[i];
}

// *****************************************************************************

// Function to concatenate two vectors

vector<double> concat(vector<double> vector1,vector<double> vector2)
{
  vector1.insert( vector1.end(), vector2.begin(), vector2.end() );
 //https://stackoverflow.com/questions/201718/concatenating-two-stdvectors
  return vector1;
}

// Function to get a vector of zeros of length l

vector<double> zeros_vec(int l){
  vector<double> Zeros_vec(l);
  for(int i=0; i<l; i++){
    Zeros_vec[i]=0;
  }
  return Zeros_vec;
}

/*
Function that adds zeros so it is length 2^M, M natural no.,
 if it is already length 2^M returns the original vector
*/

vector<double> get_vec_with_zeros(vector<double> x){
  vector<double> x_new;
  double k = x.size();
  double m = log2(k);
  double m_floor = floor(m);
  if( m_floor == m){
    x_new = x;
  }
  else{
    double l=pow(2,m_floor+1)-k;
    l = (int) l;
    vector<double> zeros_vec_add = zeros_vec(l);
    x_new = concat(x,zeros_vec_add);
  }
  return x_new;
}

//Function that raises entries of vector to power beta

vector<double> powers (vector<double> w, double beta){
  long int N = w.size();
  vector<double> Powers(N);

  for(long int i=0; i<N; i++){
    Powers[i]=pow(w[i],beta);
  }

return Powers;
}

//***********************************************************************************

//Binary search functions

//Builds binary search tree as a vector, first adds zeros so vector is 2^M

vector<double> r_tree(vector<double> rates){
  rates = get_vec_with_zeros(rates);
  double N=rates.size();
  double M = log2(N);
  double T=pow(2.0,M+1);
  N= (long int) N;
  T= (long int) T;
  vector<double> R_tree(T);

  for (long int i=0; i<N; i++){
    R_tree[N+i]=rates[i];
  }

  for(long int i=N-1; i>0; i--){
    R_tree[i] = R_tree[(i<< 1)] + R_tree[(i<<1) +1];
  }

  return R_tree;
}

/*
Finds the tree position and corresponding rates position where rate is randomly chosen from the categorical
distribution of the rates
*/

vector<long int> find_position(vector<double> R_tree, double u){
  double l = R_tree.size();
  double M = log2(l)-1;
  double N=pow(2.0,M);
  M= (long int) M;
  N= (long int) N;
  long int tree_pos=1;
  double check = R_tree[tree_pos<<1];
  vector<long int> r;
  for(long int i=1; i<=M; i++){
    if(u<check){
      tree_pos=tree_pos<<1;   // go left
      if(tree_pos<N){
      check = R_tree[tree_pos<<1];
    }}
    else{
    u=u-check;
    tree_pos=(tree_pos<<1)+1;      // go right
    if(tree_pos<N){
    check = R_tree[tree_pos<<1];
  }}}
  long int rates_pos=tree_pos-N;
  r.push_back(tree_pos);
  r.push_back(rates_pos);
  return r;
}

/* Updates the binary tree so it is the tree for the new rates (same as previous rates but
with delt added to chosen rate) */

vector<double> r_tree_update(vector<double> R_tree, long int tree_pos, double delta){
  double l = R_tree.size();
  l=(long int) l;

  for(long int i=0; i<l; i++){
    R_tree[tree_pos] = R_tree[tree_pos]+delta;
    tree_pos=tree_pos >>1;
  }
  return R_tree;
}

//*********************************************************************************

// Outputs vector into a text file file_path is path to file appended with name of file
void output_vec_in_txt(string file_path, vector<double> x){
  ofstream output_file(file_path);
  ostream_iterator<double> output_iterator(output_file, "\n");
  copy(x.begin(), x.end(), output_iterator);
}

/*
Function that runs the urn model with binary search and outputs the results into
a text file at the iterations in T
*/
vector<double> run_urn_model(vector<double> individuals_wealth,
  double W, double gamma, vector<long int> T, string file_path){

int l = T.size();
int k = 0;
long int start_s=clock();

// Mersenne twister random numbers
// https://stackoverflow.com/questions/22923551/generating-number-0-1-using-mersenne-twister-c
mt19937 generator (static_cast <unsigned> (time(0)));

string gamma_s=to_string(gamma);

// Initialise wealth to power gamma
vector<double> rates = powers(individuals_wealth, gamma);

// Initialise tree
vector<double> R_tree = r_tree(rates);

for(long int t=1; t<=T[l-1]; t++){

  double S =R_tree[1];
  uniform_real_distribution<double> dis(0.0, S);
  double u = dis(generator);
  vector<long int> x = find_position(R_tree,u);
  long int tree_pos=x[0];
  long int wealth_pos=x[1];
  double before = individuals_wealth[wealth_pos];
  individuals_wealth[wealth_pos]+= W;
  double after =  individuals_wealth[wealth_pos];
  double delta = pow(after,gamma)-pow(before,gamma);
  R_tree = r_tree_update(R_tree,tree_pos,delta);

  if(t==T[k]){
    output_vec_in_txt(file_path+"t"+to_string(k+1)+"g"+gamma_s.substr(0, 3)+".txt",individuals_wealth);
    cout << "Iteration reached: " << t << endl;
    cout << "Gamma: "<< gamma << endl;
    long int stop_s=clock();
    cout << "Time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)/60 << endl << endl;
    k++;
  }

}
return individuals_wealth;
}
