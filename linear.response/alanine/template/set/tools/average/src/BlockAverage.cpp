#include "BlockAverage.h"
#include <cmath>

void BlockAverage::
calAvg (const std::vector<double > & vec,
	double & my_avg,
	double & my_error_avg)
{
  my_avg = 0;
  for (unsigned i = 0; i < nDataUsed; ++i){
    my_avg += vec[i];
  }
  my_avg /= double(nDataUsed);

  std::vector<double > blockAverage;  
  for (unsigned i = 0; i < nBlock; ++i){
    double tmp = 0;
    for (unsigned j = 0; j < nDataInBlock; ++j){
      tmp += vec[i * nDataInBlock + j];
    }
    tmp /= double (nDataInBlock);
    blockAverage.push_back (tmp);
  }
  
  double sum2 = 0;
  for (unsigned i = 0; i < nBlock; ++i){
    sum2 += (blockAverage[i] - avg) * (blockAverage[i] - avg);
  }

  my_error_avg = sqrt (sum2 / (nBlock * (nBlock-1)));
}

void BlockAverage::
processData (const std::vector<double > & vec,
	     const unsigned & nBlock_)
{
  nBlock = nBlock_;
  unsigned nData = vec.size();
  nDataInBlock = nData / nBlock;
  nDataUsed = nDataInBlock * nBlock;

  calAvg (vec, avg, error_avg);
  
  std::vector<double > var_vec;
  for (unsigned i = 0; i < nDataUsed; ++i){
    double tmp = vec[i] - avg;
    var_vec.push_back (tmp * tmp);
  }
  calAvg (var_vec, var, error_var);

  var *= nDataUsed / (nDataUsed - 1.);
  error_var *= nDataUsed / (nDataUsed - 1.);
}



