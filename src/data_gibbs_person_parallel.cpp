#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

struct Data_Gibbs_Person : public Worker
{
  // source vector
  const RMatrix<double> interactions;
  const RMatrix<double> thresholds;
  const RMatrix<int>    data;
  const RVector<int>    no_categories;

  const int             no_persons;
  const int             iter;
  const int             no_nodes;
  const int             max_no_categories;

  // output
  RMatrix<int> augmented_data;

  // constructors
  Data_Gibbs_Person(
    const NumericMatrix interactions,
    const NumericMatrix thresholds,
    const IntegerMatrix data,
    const IntegerVector no_categories,
    const int           no_persons,
    const int           iter,
    const int           no_nodes,
    const int           max_no_categories,
          IntegerMatrix augmented_data) :
    interactions(interactions),
    thresholds(thresholds),
    data(data),
    no_categories(no_categories),
    no_persons(no_persons),
    iter(iter),
    no_nodes(no_nodes),
    max_no_categories(max_no_categories),
    augmented_data(augmented_data)
  {}

  // sample for the persons I've been asked to
  void operator()(std::size_t begin, std::size_t end)
  {

    // Rcpp::Rcout << "end: " << end << std::endl;
    double rest_score;
    double exponent;
    double cumsum;
    double u;
    int score;
    // RVector<double> probabilities(max_no_categories + 1);
    double probabilities[max_no_categories + 1];
    // NumericVector probabilities(max_no_categories + 1);

    for (std::size_t person = begin; person < end; person++ ) {

      // Fixed starting values -------------------------------------------
      // for(int node = 0; node < no_nodes; node++) {
      //   augmented_data(person, node) = data(person, node);
      // }

      for(int iteration = 0; iteration < iter; iteration++) {
        for(int node = 0; node < no_nodes; node++) {
          rest_score = 0.0;
          for(int vertex = 0; vertex < no_nodes; vertex++) {
            rest_score += augmented_data(person, vertex) * interactions(vertex, node);
          }

          cumsum = 1.0;
          probabilities[0] = 1.0;
          for(int category = 0; category < no_categories[node]; category++) {
            exponent = thresholds(node, category);
            exponent += (category + 1) * rest_score;
            cumsum += std::exp(exponent);
            probabilities[category + 1] = cumsum;
          }

          u = cumsum * R::unif_rand();

          score = 0;
          while (u > probabilities[score]) {
            score++;
          }
          augmented_data(person, node) = score;
        }
      }

    }

  }

};


// [[Rcpp::export]]
void data_gibbs_person_parallel(
    const NumericMatrix interactions,
    const NumericMatrix thresholds,
    const IntegerMatrix data,
    const IntegerVector no_categories,
    const int           no_persons,
    const int           iter,
    const int           no_nodes,
    const int           max_no_categories,
          IntegerMatrix augmented_data
)
{
  Data_Gibbs_Person Data_Gibbs_Person_worker(
      interactions,
      thresholds,
      data,
      no_categories,
      no_persons,
      iter,
      no_nodes,
      max_no_categories,
      augmented_data
    );

  // call parallelFor to do the work
  parallelFor(0, no_persons, Data_Gibbs_Person_worker);

  // output is now modified in place

}

void data_gibbs_person_serial(
    const NumericMatrix interactions,
    const NumericMatrix thresholds,
    const IntegerMatrix data,
    const IntegerVector no_categories,
    const int           no_persons,
    const int           iter,
    const int           no_nodes,
    const int           max_no_categories,
          IntegerMatrix augmented_data
)
{

  double rest_score;
  double exponent;
  double cumsum;
  double u;
  int score;
  NumericVector probabilities(max_no_categories + 1);

  for (std::size_t person = 0; person < no_persons; person++ ) {

    // Fixed starting values -------------------------------------------
    // for(int node = 0; node < no_nodes; node++) {
    //   augmented_data(person, node) = data(person, node);
    // }

    for(int iteration = 0; iteration < iter; iteration++) {
      for(int node = 0; node < no_nodes; node++) {
        rest_score = 0.0;
        for(int vertex = 0; vertex < no_nodes; vertex++) {
          rest_score += augmented_data(person, vertex) * interactions(vertex, node);
        }

        cumsum = 1.0;
        probabilities[0] = 1.0;
        for(int category = 0; category < no_categories[node]; category++) {
          exponent = thresholds(node, category);
          exponent += (category + 1) * rest_score;
          cumsum += std::exp(exponent);
          probabilities[category + 1] = cumsum;
        }

        u = cumsum * R::unif_rand();

        score = 0;
        while (u > probabilities[score]) {
          score++;
        }
        augmented_data(person, node) = score;
      }
    }

  }
}