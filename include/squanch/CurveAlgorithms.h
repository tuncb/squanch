#pragma once
#include <squanch/UtilityDetail.h>

namespace squanch {

  // Algorithm 3.1 from the NURBS book.
  // bspline control points are nurbs control points with weights 1, at least they are assumed to be.
  template <template<class, unsigned int> class NurbsType, typename T>
  typename squanch::PointContainerTypes<T>::point_type compute_bspline_point(const squanch::Model<T>& model, const NurbsType<T, 1>& bspline, T u, int span, const Eigen::Matrix<T, -1, 1>& n)
  {
    squanch::PointContainerTypes<T>::point_type C;
    C.setZero();

    const auto p = bspline.curves()[0].p();
    for (auto i = 0; i <= p; ++i) {
      C += n[i] * (model.get_point(bspline.cpi()(span - p + i)));
    }

    return C;
  }

  // Algorithm 4.1 from the NURBS book.
  template <template<class, unsigned int> class NurbsType, typename T>
  typename squanch::PointContainerTypes<T>::point_type compute_nurbs_point(const squanch::Model<T>& model, const NurbsType<T, 1>& nurbs, T u, int span, const Eigen::Matrix<T, -1, 1>& n)
  {
    squanch::PointContainerTypes<T>::point_type C = compute_bspline_point(model, nurbs, u, span, n);
    auto w = C[3];
    for (auto i = 0; i < 3; ++i) { C[i] /= w; }
    return C;
  }

  template <template<class, unsigned int> class NurbsType, typename T>
  Eigen::Matrix<T, 4, -1> compute_bspline_derivatives(const squanch::Model<T>& model, const NurbsType<T, 1>& bspline, T u, int span, const Eigen::Matrix<T, -1, 1>& n, const Eigen::Matrix<T, -1, -1>& der)
  {
    const auto num_ders = der.cols();
    const auto p = bspline.curves()[0].p();

    Eigen::Matrix<T, 4, -1> C;
    C = Eigen::MatrixXd::Zero(4, num_ders + 1);
    C.col(0) = compute_bspline_point(model, bspline, u, span, n);

    for (auto k = 1; k <= num_ders; ++k) {
      for (auto i = 0; i <= p; ++i) {
        C.col(k) += der(i, k-1) * (model.get_point(bspline.cpi()(span - p + i)));
      }
    }

    return C;
  }

  // Algorithm 4.2 from the NURBS book.
  template <typename T> void compute_nurbs_derivatives(const Eigen::Matrix<T, 4, -1>& bspline_ders, Eigen::Matrix<T, 4, -1>& CK)
  {
    const auto d = bspline_ders.cols();
    CK = Eigen::Matrix<T, 4, -1>::Zero(4, d);
    CK.setZero();

    Eigen::Matrix<T, -1, -1> Bin(d+1, d+1);
    detail::binomialCoef(Bin);
    Eigen::Matrix<T, 3, 1> v;


    for (auto k = 0; k < d; ++k) {
      v = bspline_ders.block<3, 1>(0, k);
      for (auto i = 1; i <= k; ++i) {
        v -= Bin(k, i) * bspline_ders(3, i) * CK.block<3, 1>(0, k - i);
      }
      CK.block<3, 1>(0, k) = v / bspline_ders(3, 0);
    }
  }

}