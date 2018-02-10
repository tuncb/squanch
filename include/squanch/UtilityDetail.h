#pragma once

namespace squanch {
  namespace detail {
    template <class T> void binomialCoef(Eigen::Matrix<T, -1, -1>& Bin){
      int n, k;
      Bin = Eigen::Matrix<T, -1, -1>::Zero(Bin.rows(), Bin.cols());
      // Setup the first line
      Bin(0, 0) = 1.0;
      for (k = Bin.cols() - 1; k > 0; --k) Bin(0, k) = 0.0;
      // Setup the other lines
      for (n = 0; n < Bin.rows() - 1; n++){
        Bin(n + 1, 0) = 1.0;
        for (k = 1; k < Bin.cols(); k++)
          if (n + 1 < k)
            Bin(n, k) = 0.0;
          else
            Bin(n + 1, k) = Bin(n, k) + Bin(n, k - 1);
      }
    }
  }
}