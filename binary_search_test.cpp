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

//  To run:
// g++ -std=c++11 binary_search_test.cpp urn_binary_search_functions.cpp -o main.out && ./main.out


/*
Tests the binary seacrh tree on a chosen rates vector
*/

int main(){

  vector<double> rates = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

  double S=sum(rates);

  cout << "The sum of the rates vector is: " << S << endl << endl;

  cout << "The cumulative sum of the rates vector is: " << endl;

  vector<double> cum_sum = cumsum(rates);

  vector_print_double(cum_sum);

  cout << endl <<  "The output of the rates tree is: " << endl;

  vector<double> R_tree = r_tree(rates);

  vector_print_double(R_tree);

  double S1 = R_tree[1];

  cout << endl << "The first value of the tree is the sum of the rates: " << S1 << endl << endl;

  // Mersenne twister random numbers
  // https://stackoverflow.com/questions/22923551/generating-number-0-1-using-mersenne-twister-c
  mt19937 generator (static_cast <unsigned> (time(0)));

  uniform_real_distribution<double> dis(0.0, S1);
  double u = dis(generator);


  cout << "Choose a random number u uniformly between 0 and the sum of the rates: u= " << u
  << endl << endl;

  vector<long int> positions = find_position(R_tree,u);

  double rate_chosen = rates[positions[1]];

  double rate_chosen1 = get_rate(rates,u);

  cout << "The tree position and rates position for u are " << endl;

  vector_print_long_int(positions);

  cout << endl <<  "Which corresponds to rate " << rate_chosen << endl <<
  "This should be the same output as other simpler function which is " <<  rate_chosen1 << endl;

  long int tree_pos = positions[0];
  double delta = 2;

  rates[positions[1]] +=delta;

  cout << "The rate chosen updated by delta = " << delta << " gives the new rates as "<<endl;

  vector_print_double(rates);

  vector<double> R_tree_update = r_tree_update(R_tree, tree_pos, delta);

  cout << endl <<  "The updated tree is " << endl;

  vector_print_double(R_tree_update);

  return 0;
}
