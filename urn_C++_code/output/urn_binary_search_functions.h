
#ifndef URN_BINARY_SEARCH_FUNCTIONS_H
#define URN_BINARY_SEARCH_FUNCTIONS_H


using namespace std;


void vector_print_double(vector<double> x);

void vector_print_long_int(vector<long int> x);

double sum(vector<double>x);

vector<double> cumsum(vector<double>x);

double get_rate(vector<double> rates, double u);

vector<double> concat(vector<double> vector1,vector<double> vector2);

vector<double> zeros_vec(int l);

vector<double> get_vec_with_zeros(vector<double> x);

vector<double> powers (vector<double> w, double beta);

vector<double> r_tree(vector<double> rates);

vector<long int> find_position(vector<double> R_tree, double u);

vector<double> r_tree_update(vector<double> R_tree, long int tree_pos, double delta);

void output_vec_in_txt(string file_path, vector<double> x);

vector<double> run_urn_model(vector<double> individuals_wealth,
  double W, double gamma, vector<long int> T, string file_path);


#endif
