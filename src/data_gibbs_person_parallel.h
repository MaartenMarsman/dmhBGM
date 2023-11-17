#ifndef DATA_GIBBS_PERSON_PARALLEL
#define DATA_GIBBS_PERSON_PARALLEL

#include <Rcpp.h>
using namespace Rcpp;

void data_gibbs_person_parallel(
    const NumericMatrix interactions,
    const NumericMatrix thresholds,
    const IntegerMatrix data,
    const IntegerVector no_categories,
    const int           no_persons,
    const int           iter,
    const int           no_nodes,
    const int           max_no_categories,
          IntegerMatrix output
);

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
);

#endif // DATA_GIBBS_PERSON_PARALLEL