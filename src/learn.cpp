//learn.cpp
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

//' @title Learn
//' @description An internal function that simulates learning using C++.
//' 
//' @param rep_size Learner repertoire size.
//' @param reps Demonstrator repertoires.
//' @param locs Demonstrator locations.
//' @param loc Learner location.
//' @param t_x Demonstrator attractiveness.
//' @param syl_counter Syllable counter.
//' @param num_dems Number of demonstrators.
//' @param a Level of conformity bias.
//' @param innov Innovation rate.
//' @param p_att Proportion of attractive syllables.
//'
//' @return Returns a one matrix: the Simpson's diversity from each burn-in year (column) and iteration (row).
//'
//' @export
// [[Rcpp::export(learn)]]
NumericVector learn(int rep_size, List reps, NumericVector locs, double loc, NumericVector t_x, NumericVector syl_counter, int num_dems, double a, double innov, double p_att){
  //sample demonstrators weighted by location
  NumericVector probs = 1/(Rcpp::abs(loc-locs)+1); //"+1" prevents Inf values and caps probability at 1
  IntegerVector dems = sample(reps.size(), num_dems, false, probs); //using t_x.size() instead of reps.size() because of speed

  //store collective repertoires in ordered multimap along with demonstrator attractiveness
  std::multimap<int, double> syls_4_T; //create empty multimap
  for(int z = 1; z <= num_dems; z++){ //for each demonstrator
    std::vector<int> temp = reps[dems[z-1]-1]; //temporarily store their repertoire as a vector
    for(int y = 1; y <= temp.size(); y++){
      syls_4_T.insert({temp[y-1], t_x[dems[z-1]-1]}); //and then add each syllable in the repertoire as keys, and the demonstrator's attractiveness as repeated values
    }
  }

  //create list of all syl ids
  NumericVector syl_ids; //create empty vector
  for(auto h : syls_4_T){ //for each entry in the multimap
    syl_ids.push_back(h.first); //add the key to the vector
  }
  NumericVector all_syls = syl_ids; //make copy with all before sort unique

  //sort unique syls via https://www.javatpoint.com/cpp-algorithm-unique-function (faster than Rcpp function)
  std::sort(syl_ids.begin(), syl_ids.end());
  auto last = std::unique(syl_ids.begin(), syl_ids.end());
  syl_ids.erase(last, syl_ids.end());

  //calculate average demonstrator attractiveness per syl (T_x)
  NumericVector T_x; //create empty vector
  for(int w = 1; w <= syl_ids.size(); w++){ //for each syllable id
    auto range = syls_4_T.equal_range(syl_ids[w-1]); //get position of key

    //store values of key via https://www.geeksforgeeks.org/multimap-equal_range-in-c-stl/
    std::vector<double> values;
    for(auto v = range.first; v != range.second; v++){
      values.push_back(v -> second);
    }

    //calculate mean via https://stackoverflow.com/questions/28574346/find-average-of-input-to-vector-c (faster than Rcpp function)
    double mean = accumulate(values.begin(), values.end(), 0.0)/values.size();
    T_x.push_back(mean); //add mean attractiveness to vector
  }

  //calculate biased frequencies (F_x)
  NumericVector F_x; //create empty vector
  for(int q = 1; q <= syl_ids.size(); q++){ //for each syllable id
    F_x.push_back(std::pow(syls_4_T.count(syl_ids[q-1]), a)); //count how many syllables there are and raise it to exponent a, and add to vector
  }

  //retrieve syl attractiveness (M_x)
  NumericVector M_x = syl_counter[syl_ids-1];

  //calculate overall probability of adoption (P_x)
  NumericVector P_x; //create empty vector
  for(int u = 1; u <= syl_ids.size(); u++){ //for each syllable id
    P_x.push_back(T_x[u-1]*F_x[u-1]*M_x[u-1]); //multiply T_x, F_x, and M_x and add to vector
  }

  //sample possible syls (given P_x)
  NumericVector new_rep = sample(syl_ids, rep_size, true, P_x); //"replacement = true" allows code to run for individuals who have a repertoire size larger that total reps of dems (only problem when dems = 2)

  //generate innovation
  LogicalVector to_replace = sample(LogicalVector {true, false}, new_rep.size(), true, NumericVector {innov, 1-innov}); //generate vector of T/F with given probability of innovation for each learned syllable
  int sum_replace = sum(to_replace); //sum these
  if(sum_replace > 0){ //if innovation is to occur
    new_rep[to_replace] = sample(syl_counter.size(), sum_replace, false); //replace with random syllables from population repertoire
  }
  
  return new_rep;
}
