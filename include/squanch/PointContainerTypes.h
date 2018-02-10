#pragma once
#include <Eigen/StdVector>
#include <vector>
#include <map>
#include <Eigen/core>


namespace squanch {

template <typename T> struct PointContainerTypes
{
  typedef Eigen::Matrix<T, 4, 1> point_type;

  typedef std::map< int, Eigen::Matrix<T, 4, 1>,  std::less<int>, Eigen::aligned_allocator < std::pair<const int, Eigen::Matrix<T, 4, 1>> >>  map_type;

  typedef std::vector<Eigen::Matrix<T, 4, 1>, Eigen::aligned_allocator<Eigen::Matrix<T, 4, 1>>> vector_type;
};

}
