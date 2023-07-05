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
// g++ -O9 -std=c++11 urn_binary_search.cpp urn_binary_search_functions.cpp -o main.out && ./main.out

int main(){

double N=pow(10.0,3); //number of agents
N = (long int) N;
double W0=1;  // initial wealth
long int t1= (long int) pow(10.0,6); // number of iterations
long int t2= (long int) pow(10.0,7); // number of iterations
long int t3= (long int) pow(10.0,8); // number of iterations
long int t4= (long int) pow(10.0,9); // number of iterations
vector<long int> T = {t1,t2,t3,t4};
vector<double> gammas = {1.0,1.1,1.2,1.3};

double W=1.0; // wealth packet

//Run urn model for gammas in the gamma vector

for(int i=0;i<gammas.size();i++){
// initialise wealth, every agent with W0
vector<double> individuals_wealth(N,W0);
// double s = sum(individuals_wealth);
// cout << s << endl;
individuals_wealth=run_urn_model(individuals_wealth, W, gammas[i], T, "./output/");
}

return 0;

}
