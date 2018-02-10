#pragma once
#include <Eigen\geometry>
#include <squanch/Model.h>
#include <squanch/NurbsCurve.h>
#include <squanch\Nurbs.h>
#include <squanch\PointContainerTypes.h>
#include <squanch\IgaAlgorithms.h>

namespace squanch {

namespace detail {

  template <typename T, typename Fun> void for_each_cp(const squanch::Model<T>& model, const blitz::Array<int, 1>& cpi, const blitz::TinyVector<int, 1>& bounds, Fun fun)
  {
    for (int i = 0; i < bounds[0]; ++i) {
      fun(model.get_point(cpi(i)));
    }
  }

  template <typename T, typename Fun> void for_each_cp(const squanch::Model<T>& model, const blitz::Array<int, 2>& cpi, const blitz::TinyVector<int, 2>& bounds, Fun fun)
  {
    for (int i = 0; i < bounds[0]; ++i) {
      for (int j = 0; j < bounds[1]; ++j) {
        fun(model.get_point(cpi(i, j)));
      }
    }
  }

  template <typename T, typename Fun> void for_each_cp(const squanch::Model<T>& model, const blitz::Array<int, 3>& cpi, const blitz::TinyVector<int, 3>& bounds, Fun fun)
  {
    for (int i = 0; i < bounds[0]; ++i) {
      for (int j = 0; j < bounds[1]; ++j) {
        for (int k = 0; k < bounds[2]; ++k) {
          fun(model.get_point(cpi(i, j, k)));
        }
      }
    }
  }
  
  template <typename T, typename Fun> void for_each_cp(squanch::Model<T>& model, const blitz::Array<int, 1>& cpi, const blitz::TinyVector<int, 1>& bounds, Fun fun)
  {
    for (int i = 0; i < bounds[0]; ++i) {
      fun(model.get_point(cpi(i)));
    }
  }

  template <typename T, typename Fun> void for_each_cp(squanch::Model<T>& model, const blitz::Array<int, 2>& cpi, const blitz::TinyVector<int, 2>& bounds, Fun fun)
  {
    for (int i = 0; i < bounds[0]; ++i) {
      for (int j = 0; j < bounds[1]; ++j) {
        fun(model.get_point(cpi(i, j)));
      }
    }
  }

  template <typename T, typename Fun> void for_each_cp(squanch::Model<T>& model, const blitz::Array<int, 3>& cpi, const blitz::TinyVector<int, 3>& bounds, Fun fun)
  {
    for (int i = 0; i < bounds[0]; ++i) {
      for (int j = 0; j < bounds[1]; ++j) {
        for (int k = 0; k < bounds[2]; ++k) {
          fun(model.get_point(cpi(i, j, k)));
        }
      }
    }
  }
}


template <typename T> struct TransformationTypes
{
  typedef Eigen::Transform<T, 3, Eigen::Affine> transformation;
};

//template <typename T, unsigned int NDIM, typename Fun> void for_each_cp(squanch::Model<T>& model, squanch::Nurbs<T, NDIM>& nurbs, Fun fun)
//{
//  auto&& cpi = nurbs.cpi();
//  auto shape = nurbs.cpi().shape();
//  for (unsigned int i = 0; i < NDIM; ++i) if (nurbs.curves()[i].is_closed()) shape[i] -= 1;
//
//  detail::for_each_cp(model, cpi, shape, fun);
//}

template <typename ModelType, typename T, unsigned int NDIM, typename Fun> void for_each_cp(ModelType& model, squanch::Nurbs<T, NDIM>& nurbs, Fun fun)
{
  auto&& cpi = nurbs.cpi();
  auto shape = nurbs.cpi().shape();
  for (unsigned int i = 0; i < NDIM; ++i) if (nurbs.curves()[i].is_closed()) shape[i] -= 1;

  detail::for_each_cp(model, cpi, shape, fun);
}

template <typename T, unsigned int NDIM> void find_center(const squanch::Model<T>& model, squanch::Nurbs<T, NDIM>& nurbs, typename PointContainerTypes<T>::point_type& center)
{
  center = PointContainerTypes<T>::point_type::Zero();
  for_each_cp(model, nurbs, [&center](const PointContainerTypes<T>::point_type& po) {
    PointContainerTypes<T>::point_type p = po;
    from_homogeneous(p);
    center += p;
  });

  auto&& cpi = nurbs.cpi();
  auto shape = nurbs.cpi().shape();
  for (unsigned int i = 0; i < NDIM; ++i) if (nurbs.curves()[i].is_closed()) shape[i] -= 1;
  size_t n = 1;
  for (unsigned int i = 0; i < NDIM; ++i) n *= shape[i];

  center /= n;
  center.w() = 1;
}

template <typename T, unsigned int NDIM, typename TransType> void transform(squanch::Model<T>& model, squanch::Nurbs<T, NDIM>& nurbs, const TransType& t)
{
  for_each_cp(model, nurbs, [&t](typename PointContainerTypes<T>::point_type& p) { p = t * p; });
}

template <typename T, unsigned int NDIM> void translate(squanch::Model<T>& model, squanch::Nurbs<T, NDIM>& nurbs, T x, T y, T z)
{
  Eigen::Transform<T, 3, Eigen::Affine> t(Eigen::Translation<T, 3>(x,y,z));
  squanch::transform(model, nurbs, t);
}

template <typename T, unsigned int NDIM, typename TransType> void transform_around (squanch::Model<T>& model, squanch::Nurbs<T, NDIM>& nurbs, const TransType& t, T x, T y, T z)
{
  Eigen::Transform<T, 3, Eigen::Affine> to_origin(Eigen::Translation<T, 3>(-x, -y, -z)); 
  Eigen::Transform<T, 3, Eigen::Affine> from_origin(Eigen::Translation<T, 3>(x, y, z)); 
  Eigen::Transform<T, 3, Eigen::Affine> tot;
  tot = from_origin * t * to_origin;

  squanch::transform(model, nurbs, tot);
}

}