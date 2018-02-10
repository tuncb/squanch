#pragma once
#include <squanch/Array2D.h>

namespace squanch {


  template <template<class, unsigned int> class NurbsType, typename T>
	typename squanch::PointContainerTypes<T>::point_type compute_bspline_point(const squanch::Model<T>& model, const NurbsType<T, 2>& bspline, const std::array<T, 2>& u, const std::array<unsigned int, 2>& span, const Eigen::Matrix<T, -1, 1>& n)
  {
    auto p1 = bspline.curves()[0].p();
    auto p2 = bspline.curves()[1].p();
    auto&& cpi = bspline.cpi();

    squanch::PointContainerTypes<T>::point_type point = Eigen::Matrix<T, 4, 1>::Zero();

    unsigned int c = 0;
    for (unsigned int i = 0; i <= p2; ++i) {
      for (unsigned int j = 0; j <= p1; ++j) {
        point += model.get_point(cpi((int)(span[0] - p1 + j), (int)(span[1] - p2 + i))) * n[c++];
      }
    }

    return point;
  }


  template <template<class, unsigned int> class NurbsType, typename T>
  typename squanch::PointContainerTypes<T>::point_type compute_nurbs_point(const squanch::Model<T>& model, const NurbsType<T, 2>& nurbs, const std::array<T, 2>& u, const std::array<unsigned int, 2>& span, const Eigen::Matrix<T, -1, 1>& n)
  {
    squanch::PointContainerTypes<T>::point_type C = compute_bspline_point(model, nurbs, u, span, n);
    auto w = C[3];
    for (auto i = 0; i < 3; ++i) { C[i] /= w; }
    return C;
  }

  template <template<class, unsigned int> class NurbsType, typename T>
  Eigen::Matrix<T, 4, -1> compute_bspline_derivatives_on_curve
		(squanch::Model<T>& model, const NurbsType<T, 2>& bspline, const std::array<T, 2>& u, unsigned int constant_dim, const std::array<unsigned int, 2>& span, const Eigen::Matrix<T, -1, 1>& n_2d, const Eigen::Matrix<T, -1, -1>& der)
  {
    unsigned int curve_dim;
    if (constant_dim == 0) curve_dim = 1;
    else curve_dim = 0;

    const auto num_ders = der.cols();
    const auto p = bspline.curves()[curve_dim].p();

    Eigen::Matrix<T, 4, -1> C;
    C = Eigen::MatrixXd::Zero(4, num_ders + 1);
    C.col(0) = compute_bspline_point(model, bspline, u, span, n_2d);

    blitz::TinyVector<int, 2> index;
    index[constant_dim] = span[constant_dim];

    for (auto k = 1; k <= num_ders; ++k) {
      for (auto i = 0; i <= p; ++i) {
        index[curve_dim] = span[curve_dim] - p + i;
        C.col(k) += der(i, k - 1) * model.get_point(bspline.cpi()(index));
      }
    }

    return C;
  }

  template <template<class, unsigned int> class NurbsType, typename T, unsigned int NDIM>   Eigen::Matrix<T, 4, -1> compute_bspline_derivative_points
    (squanch::Model<T>& model, const NurbsType<T, NDIM>& spline, typename const NurbsType<T, NDIM>::ParametricPointType& u, typename const NurbsType<T, NDIM>::SpanType& span, const Eigen::Matrix<T, -1, -1>& der)
  {
    const auto num_ders = der.cols();
    auto nodeids = gather_node_ids(spline, span);

    Eigen::Matrix<T, 4, -1> C;
    C = Eigen::MatrixXd::Zero(4, num_ders);
    for (auto k = 1; k <= num_ders; ++k) {
      for (auto&& id : nodeids) {
        C.col(k) += der(i, k) * model.get_point(id);
      }
    }
    return C;
  }

  // Algorithm A3.6 from the NURBS Book.
  template <template<class, unsigned int> class NurbsType, typename T> void 
    compute_bspline_surface_derivatives(const squanch::Model<T>& model, const NurbsType<T, 2>& spline, unsigned int nder, const std::array<T, 2>& pos, PointArray2D<T>& SKL)
  {
    const auto d = nder;
    const auto u = pos[0];
    const auto v = pos[1];

    const auto n = spline.curves()[0].n();
    const auto p = spline.curves()[0].p();
    
    const auto m = spline.curves()[1].n();
    const auto q = spline.curves()[1].p();

    // this part is different from the book
    // the book zeros only the derivatives that will be zero.
    // I think it misses some partial derivative combinations
    // here I zero everything.
    const auto du = std::min(d, p);
    const auto dv = std::min(d, q);

    for (auto k = 0u; k <= d; ++k) {
      for (auto l = 0u; l <= d; ++l) {
        SKL(k, l).setZero();
      }
    }

    // notice that Nu and Nv matrices are the transpose of their book definition
    Eigen::Matrix<T, -1, 1> Nu_basis(p + 1);
    Eigen::Matrix<T, -1, -1> Nu_der(p + 1, du);
    Eigen::Matrix<T, -1, -1> Nu(p + 1, du + 1);

    auto uspan = find_span(spline.curves()[0].knots(), spline.curves()[0].p(), u);
    compute_b_ders_basis(spline.curves()[0], uspan, u, du, Nu_basis, Nu_der);
    Nu.col(0) = Nu_basis;
    Nu.block(0, 1, p + 1, du) = Nu_der;

    Eigen::Matrix<T, -1, 1> Nv_basis(q + 1);
    Eigen::Matrix<T, -1, -1> Nv_der (q + 1, dv);
    Eigen::Matrix<T, -1, -1> Nv(q + 1, dv + 1);

    auto vspan = find_span(spline.curves()[1].knots(), spline.curves()[1].p(), v);
    compute_b_ders_basis(spline.curves()[1], vspan, v, dv, Nv_basis, Nv_der);
    Nv.col(0) = Nv_basis;
    Nv.block(0, 1, q + 1, dv) = Nv_der;
    
    typename PointContainerTypes<T>::vector_type temp(q + 1);

    for (auto k = 0; k <= du; ++k) {
      for (auto s = 0; s <= q; ++s) {
        temp[s].setZero();
        for (auto r = 0; r <= p; ++r) {
          auto&& pnt = model.get_point(spline.cpi()((int)(uspan - p + r), (int)(vspan - q + s)));
          temp[s] += Nu(r, k) * pnt;
        }
      }
      auto dd = std::min(d - k, dv);
      for (auto l = 0; l <= dd; ++l) {
        SKL(k, l).setZero();
        for (auto s = 0; s <= q; ++s){
          SKL(k, l) += Nv(s, l) * temp[s];
        }
      }
    }
  }

  // Algorithm A4.4 from the NURBS Book.
  template <typename T> void
    compute_nurbs_surface_derivatives(const squanch::Model<T>& model, PointArray2D<T>& bspline_SKL, PointArray2D<T>& SKL)
  {
    auto shape = bspline_SKL.shape();
    SKL.resize(shape);
    const auto d = shape[0] - 1; 
    for (auto i = 0u; i < shape[0]; ++i) {
      for (auto j = 0u; j < shape[1]; ++j) {
        SKL(i, j).setZero();
      }
    }

    Eigen::Matrix<T, -1, -1> Bin(d + 1, d + 1);
    detail::binomialCoef(Bin);
    Eigen::Matrix<T, 3, 1> v;
    Eigen::Matrix<T, 3, 1> v2;

    for (auto k = 0; k <= d; ++k) {
      for (auto l = 0; l <= d - k; ++l) {
        v = bspline_SKL(k, l).block<3, 1>(0, 0);
        for (auto j = 0; j <= l; ++j) {
          v -= Bin(l, j) * bspline_SKL(0, j)(3, 0) * SKL(k, l - j).block<3,1>(0, 0);
        }
        
        for (auto i = 0; i <= k; ++i) {
          v -= Bin(k, i) * bspline_SKL(i, 0)(3, 0) * SKL(k - i, l).block<3,1>(0,0);
          v2.setZero();
          for (auto j = 0; j <= l; ++j) {
            v2 += Bin(l, j) * bspline_SKL(i, j)(3, 0) * SKL(k - i, l - j).block<3,1>(0,0);
          }
          v -= Bin(k, i) * v2;
        }

        SKL(k, l).block<3,1>(0, 0) = v / bspline_SKL(0, 0)(3, 0);
      }
    }
  }

}