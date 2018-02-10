#pragma once
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <Eigen/core>
#include <squanch\Nurbs.h>
#include <squanch\PointContainerTypes.h>

namespace squanch {

/// <summary>
/// Finds unique knots from specified knots.
/// </summary>
/// <param name="knots">The knots.</param>
/// <returns>unique knot vector</returns>
template <typename T> std::vector<T> unique_knots(const std::vector<T>& knots)
{
  using namespace std;
  auto iter = unique(begin(knots), end(knots()), [](T a, T b) {return fabs(a - b) < (T)0.00000001;});
  return vector<T>(begin(knots), iter);
}

 /// <summary>
/// converts a non-homogeneous point to a homogeneous one.
/// </summary>
/// <param name="vec">The non-homogeneous point.</param>
template <typename T> void to_homogeneous(Eigen::Matrix<T, 4, 1>& vec)
{
  auto w = vec.w();
  vec *= w;
  vec[3] = w;
}

/// <summary>
/// converts all points in a map to homogeneous ones.
/// </summary>
/// <param name="map">The map.</param>
template <typename MapType> void map_to_homogeneous(MapType& map) 
{
  for (auto&& pair : map) {
    to_homogeneous(pair.second);
  }
}

template <typename MapType> void map_from_homogeneous(MapType& map) 
{
  for (auto&& pair : map) {
    from_homogeneous(pair.second);
  }
}


/// <summary>
/// converts a homogeneous point to a non-homogeneous one.
/// </summary>
/// <param name="vec">The homogeneous point.</param>
template <typename T> void from_homogeneous(Eigen::Matrix<T, 4, 1>& vec)
{
  auto w = vec.w();
  vec *= (T)1/w;
  vec[3] = w;
}


template <template<class, unsigned int> class NurbsType, typename T> std::vector<int> gather_node_ids(const NurbsType<T, 1>& spline, const unsigned int span)
{
  auto p1 = nurbs.curves()[0].p();
  auto&& cpi = nurbs.cpi();
  std::vector<int> nodeids(p + 1);

  for (unsigned int i = 0; i <= p; ++i) {
    nodeids[i] = cp.find(cpi(span - p + i));
  }
  return nodeids;
}

template <template<class, unsigned int> class NurbsType, typename T> std::vector<int> gather_node_ids(const NurbsType<T, 2>& spline, const std::array<unsigned int, 2>& span)
{
  auto p1 = nurbs.curves()[0].p();
  auto p2 = nurbs.curves()[1].p();
  auto numnodes = (p1 + 1) * (p2 + 1);
  auto&& cpi = nurbs.cpi();
  std::vector<int> nodeids(numnodes);

  unsigned int c = 0;
  for (unsigned int i = 0; i <= p2; ++i) {
    for (unsigned int j = 0; j <= p1; ++j) {
      nodeids[c++] = cpi((int)(span[0] - p1 + j), (int)(span[1] - p2 + i));
    }
  }

  return nodeids;
}

template <typename T, template<class, unsigned int> class NurbsType> std::vector<int> gather_node_ids(const NurbsType<T, 3>& spline, const std::array<unsigned int, 3>&  span)
{
  auto p1 = nurbs.curves()[0].p();
  auto p2 = nurbs.curves()[1].p();
  auto p3 = nurbs.curves()[2].p();
  auto numnodes = (p1 + 1) * (p2 + 1) * (p3 + 1);
  auto&& cpi = nurbs.cpi();
  std::vector<int> nodeids(numnodes);

  unsigned int c = 0;
  for (unsigned int k = 0; k <= p3; ++k) {
    for (unsigned int i = 0; i <= p2; ++i) {
      for (unsigned int j = 0; j <= p1; ++j) {
        nodeids[c++] = cpi((int)(span[0] - p1 + j), (int)(span[1] - p2 + i), (int)(span[2] - p3 + k));
      }
    }
  }
  return nodeids;
}


/// <summary>
/// Gathers weights for the specified span for a 1D NURBS.
/// </summary>
/// <param name="nurbs">The nurbs.</param>
/// <param name="cp">Control Point map.</param>
/// <param name="span">The span.</param>
/// <param name="weights">The weights.</param>
template <typename T, template<class, unsigned int> class NurbsType> void gather_weights(const squanch::Model<T>& model, const NurbsType<T, 1>& nurbs, const unsigned int span, Eigen::Matrix<T, -1, 1>& weights) {
  auto p = nurbs.curves()[0].p();
  auto&& cpi = nurbs.cpi();
  
  weights.resize(p+1);
  for (unsigned int i = 0; i <= p; ++i) {
    weights[i] = model.get_point(cpi(span - p + i)).w();
  }
}

/// <summary>
/// Gathers weights for the specified span for a 2D NURBS.
/// </summary>
/// <param name="nurbs">The nurbs.</param>
/// <param name="cp">Control Point map.</param>
/// <param name="span">The span.</param>
/// <param name="weights">The weights.</param>>
template <typename T, template<class, unsigned int> class NurbsType> void gather_weights(const squanch::Model<T>& model, const NurbsType<T, 2>& nurbs,
                                         const std::array<unsigned int, 2> span, Eigen::Matrix<T, -1, 1>& weights)
{
  auto p1 = nurbs.curves()[0].p();
  auto p2 = nurbs.curves()[1].p();
  auto&& cpi = nurbs.cpi();
  
  weights.resize((p1+1)*(p2+1));

  unsigned int c = 0;
  for (unsigned int i = 0; i <= p2; ++i) {
    for (unsigned int j = 0; j <= p1; ++j) {
      weights[c++] = model.get_point(cpi((int)(span[0] - p1 + j),(int)(span[1] - p2 + i))).w();
    }
  }
}

/// <summary>
/// Gathers weights for the specified span for a 3D NURBS.
/// </summary>
/// <param name="nurbs">The nurbs.</param>
/// <param name="cp">Control Point map.</param>
/// <param name="span">The span.</param>
/// <param name="weights">The weights.</param>
template <typename T, template<class, unsigned int> class NurbsType> void gather_weights(const squanch::Model<T>& model, const NurbsType<T, 3>& nurbs,
                                         const std::array<unsigned int, 3> span, Eigen::Matrix<T, -1, 1>& weights)
{
  auto p1 = nurbs.curves()[0].p();
  auto p2 = nurbs.curves()[1].p();
  auto p3 = nurbs.curves()[2].p();
  auto&& cpi = nurbs.cpi();
  
  weights.resize((p1+1)*(p2+1)*(p3+1));

  unsigned int c = 0;
  for (unsigned int k = 0; k <= p3; ++k) {
    for (unsigned int i = 0; i <= p2; ++i) {
      for (unsigned int j = 0; j <= p1; ++j) {
        weights[c++] = model.get_point(cpi((int)(span[0] - p1 + j),(int)(span[1] - p2 + i),(int)(span[2] - p3 + k))).w();
      }
    }
  }
}

/// <summary>
/// Computes model coordinates from basis functions of a parametric coordinate for the given 1D NURBS.
/// </summary>
/// <param name="nurbs">The nurbs.</param>
/// <param name="cp">Control Point map.</param>
/// <param name="span">The span of the parametric coordinate.</param>
/// <param name="r">Basis functions. For homogeoneous control points b-splibe basis functions should be given. 
/// Otherwise NURBS basis functions should be given.</param>
/// <param name="point">The calculated point.</param>
template <template<class, unsigned int> class NurbsType, typename T>  void compute_model_coordinate(const squanch::Model<T>& model, const NurbsType<T, 1>& nurbs,
                                         const unsigned int span, const Eigen::Matrix<T, -1, 1>& r, Eigen::Matrix<T, 4, 1>& point)
{  
  auto p = nurbs.curves()[0].p();
  auto&& cpi = nurbs.cpi();

  point = Eigen::Matrix<T, 4, 1>::Zero();

  for (unsigned int i = 0; i <= p; ++i) {
    point += (model.get_point(cpi(span - p + i))) * r[i];
  }
}

/// <summary>
/// Computes model coordinates from basis functions of a parametric coordinate for the given 2D NURBS.
/// </summary>
/// <param name="nurbs">The nurbs.</param>
/// <param name="cp">Control Point map.</param>
/// <param name="span">The span of the parametric coordinate.</param>
/// <param name="r">Basis functions. For homogeoneous control points b-splibe basis functions should be given. 
/// Otherwise NURBS basis functions should be given.</param>
/// <param name="point">The calculated point.</param>>
template <template<class, unsigned int> class NurbsType, typename T> void compute_model_coordinate(const squanch::Model<T>& model, const NurbsType<T, 2>& nurbs,
                                         const std::array<unsigned int, 2>& span, const Eigen::Matrix<T, -1, 1>& r, Eigen::Matrix<T, 4, 1>& point)
{  
  auto p1 = nurbs.curves()[0].p();
  auto p2 = nurbs.curves()[1].p();
  auto&& cpi = nurbs.cpi();

  point = Eigen::Matrix<T, 4, 1>::Zero();

  unsigned int c = 0;
  for (unsigned int i = 0; i <= p2; ++i) {
    for (unsigned int j = 0; j <= p1; ++j) {
      point += (model.get_point(cpi((int)(span[0] - p1 + j),(int)(span[1] - p2 + i)))) * r[c++];
    }
  } 
}

/// <summary>
/// Computes model coordinates from basis functions of a parametric coordinate for the given 3D NURBS.
/// </summary>
/// <param name="nurbs">The nurbs.</param>
/// <param name="cp">Control Point map.</param>
/// <param name="span">The span of the parametric coordinate.</param>
/// <param name="r">Basis functions. For homogeoneous control points b-splibe basis functions should be given. 
/// Otherwise NURBS basis functions should be given.</param>
/// <param name="point">The calculated point.</param>
template <template<class, unsigned int> class NurbsType, typename T> void compute_model_coordinate(const squanch::Model<T>& model, const NurbsType<T, 3>& nurbs,
                                         const std::array<unsigned int, 3>& span, const Eigen::Matrix<T, -1, 1>& r, Eigen::Matrix<T, 4, 1>& point)
{  
  auto p1 = nurbs.curves()[0].p();
  auto p2 = nurbs.curves()[1].p();
  auto p3 = nurbs.curves()[2].p();
  auto&& cpi = nurbs.cpi();

  point = Eigen::Matrix<T, 4, 1>::Zero();

  unsigned int c = 0;
  for (unsigned int k = 0; k <= p3; ++k) {
    for (unsigned int i = 0; i <= p2; ++i) {
      for (unsigned int j = 0; j <= p1; ++j) {
        point += (model.get_point(cpi((int)(span[0] - p1 + j),(int)(span[1] - p2 + i),(int)(span[2] - p3 + k)))) * r[c++];
      }
    }
  }

}

/// <summary>
/// Find the span of a parametric point u.
/// </summary>
/// <param name="knots">The knot vector.</param>
/// <param name="p">Degree of the curve.</param>
/// <param name="u">The parametric point u.</param>
/// <returns></returns>
template <typename T> unsigned int find_span(const std::vector<T>& knots, unsigned int p, T u)
{
  if (u == knots.back()) return knots.size() - p - 2;

  auto iter = std::upper_bound(knots.cbegin(), knots.cend(), u);
  return std::distance(knots.cbegin(), iter) - 1;
}

template <template<class, unsigned int> class NurbsType, typename T> unsigned int find_span(const NurbsType<T, 1>& spline, T u)
{
  auto&& knots = spline.curves()[0].knots();
  auto p = spline.curves()[0].p();

  if (u == knots.back()) return knots.size() - p - 2;

  auto iter = std::upper_bound(knots.cbegin(), knots.cend(), u);
  return std::distance(knots.cbegin(), iter) - 1;
}

template <template<class, unsigned int> class NurbsType, typename T> std::array<unsigned int, 2> find_span(const NurbsType<T, 2>& spline, const std::array<T, 2>& u)
{
  const auto ndim = 2u;
  std::array<unsigned int, ndim> span;
  for (auto i = 0u; i < ndim; ++i) {
    auto&& curve = spline.curves()[i];
    span[i] = find_span(curve.knots(), curve.p(), u[i]);
  }
  return span;
}

template <template<class, unsigned int> class NurbsType, typename T> std::array<unsigned int, 3> find_span(const NurbsType<T, 3>& spline, const std::array<T, 3>& u)
{
  const auto ndim = 3u;
  std::array<unsigned int, ndim> span;
  for (auto i = 0u; i < ndim; ++i) {
    auto&& curve = spline.curves()[i];
    span[i] = find_span(curve.knots(), curve.p(), u[i]);
  }
  return span;
}


template <typename T> struct ComputeBDers1BasisTemporary
{
  ComputeBDers1BasisTemporary() {}
  ComputeBDers1BasisTemporary(const squanch::NurbsCurve<T>& curve)
  {
    this->reset(curve);
  }

  void reset(const squanch::NurbsCurve<T>& curve)
  {
    auto p = curve.p();

    left.resize(p + 1);
    right.resize(p + 1);
    ndu.resize(p + 1, p + 1);
    a.resize(p + 1, p + 1);
  }

  Eigen::Matrix<T, -1, 1> left;
  Eigen::Matrix<T, -1, 1> right;

  Eigen::Matrix<T, -1, -1> ndu;
  Eigen::Matrix<T, -1, -1> a;
};

/// <summary>
/// Computes 1st derivatives and B-Spline basis functions at a parametric point for the given curve.
/// </summary>
/// <param name="curve">The curve.</param>
/// <param name="span">The span of the parametric point u.</param>
/// <param name="u">The parametric point.</param>
/// <param name="r">OUT Vector(p+1) holding B-Spline basis functions.</param>
/// <param name="der">OUT Matrix(p+1, 1)  1st derivatives of the B-Spline basis functions @u.</param>
template <typename T> void compute_b_ders1_basis(const squanch::NurbsCurve<T>& curve, const unsigned int span, T u, ComputeBDers1BasisTemporary<T>& temporary, Eigen::Matrix<T, -1, 1>& r, Eigen::Matrix<T, -1, 1>& der)
{
  /* Modified from Algorithm A2.3 from the NURBS book */

  auto p = curve.p();
  auto&& knots = curve.knots();

  auto&& left = temporary.left;
  auto&& right = temporary.right;
  auto&& ndu = temporary.ndu;
  auto&& a = temporary.a;

  ndu(0,0) = (T)1;
  for (unsigned int i = 1; i <= p; ++i) {
    left[i]  = u - knots[span + 1 - i];
    right[i] = knots[span + i] - u;    
    T saved = (T)0;
    for (unsigned int j = 0; j < i; ++j) {
      ndu(i,j)  = right[j+1] + left[i-j];
      auto temp = ndu(j,i-1) / ndu(i,j);
      ndu(j,i) = saved + right[j+1]*temp;
      saved=left[i-j]*temp;      
    }
    ndu(i,i) = saved;
  }

  r = ndu.col(p);

  const unsigned int k = 1;
  
  for( int i=0; i<=p; i++ ) {
    unsigned int s1=0, s2=1;   
    a(0,0) = (T)1;

    T d = (T)0;
    int rk = i-k;
    int pk = p-k;

    if( i >= k ) {
      a(s2,0)= a(s1,0) / ndu(pk+1,rk);
      d = a(s2,0)*ndu(rk,pk);
    }

    auto j1 = rk >= -1 ? 1 : -rk;
    auto j2 = (i-1<=pk) ? k-1 : p-i;

    for( unsigned int j=j1; j<=j2; j++) {
      a(s2,j) = (a(s1,j)- a(s1,j-1)) / ndu(pk+1,rk+j);
      d += a(s2,j)*ndu(rk+j,pk);
    }
    if(i <=pk) {
      a(s2,k) = -a(s1,k-1)/ndu(pk+1,i);
      d += a(s2,k)*ndu(i,pk);
    }
    der[i] = d;
    s1=s2; s2=j2;  
  }

  der *= p;
}

/// <summary>
/// Computes 1st derivatives and B-Spline basis functions at a parametric point for the given curve.
/// </summary>
/// <param name="curve">The curve.</param>
/// <param name="span">The span of the parametric point u.</param>
/// <param name="u">The parametric point.</param>
/// <param name="r">OUT Vector(p+1) holding B-Spline basis functions.</param>
/// <param name="der">OUT Matrix(p+1, num_ders)  ist derivatives of the B-Spline basis functions @u.</param>
template <typename T> void compute_b_ders_basis(const squanch::NurbsCurve<T>& curve, const unsigned int span, T u, ComputeBDers1BasisTemporary<T>& temporary, unsigned int num_ders,  Eigen::Matrix<T, -1, 1>& r, Eigen::Matrix<T, -1, -1>& der)
{
  /* Modified from Algorithm A2.3 from the NURBS book */


  auto p = curve.p();
  auto&& knots = curve.knots();
  
  auto&& left = temporary.left;
  auto&& right = temporary.right;
  auto&& ndu = temporary.ndu;
  auto&& a = temporary.a;

  ndu(0, 0) = (T)1;
  for (unsigned int i = 1; i <= p; ++i) {
    left[i] = u - knots[span + 1 - i];
    right[i] = knots[span + i] - u;
    T saved = (T)0;
    for (unsigned int j = 0; j < i; ++j) {
      ndu(i, j) = right[j + 1] + left[i - j];
      auto temp = ndu(j, i - 1) / ndu(i, j);
      ndu(j, i) = saved + right[j + 1] * temp;
      saved = left[i - j] * temp;
    }
    ndu(i, i) = saved;
  }

  r = ndu.col(p);

  for (auto ri = 0; ri <= p; ++ri) {
    unsigned int s1 = 0, s2 = 1;
    a(0, 0) = (T)1;

    for (auto k = 1; k <= num_ders; ++k) {
      T d = (T)0;
      int rk = ri - k;
      int pk = p - k;

      if (ri >= k) {
        a(s2, 0) = a(s1, 0) / ndu(pk + 1, rk);
        d = a(s2, 0)*ndu(rk, pk);
      }

      int j1, j2;
      if (rk >= -1) {
        j1 = 1;
      }
      else {
        j1 = -rk;
      }
      if (ri - 1 <= pk) {
        j2 = k - 1;
      }
      else {
        j2 = p - ri;
      }

      for (unsigned int j = j1; j <= j2; j++) {
        a(s2, j) = (a(s1, j) - a(s1, j - 1)) / ndu(pk + 1, rk + j);
        d += a(s2, j)*ndu(rk + j, pk);
      }
      if (ri <= pk) {
        a(s2, k) = -a(s1, k - 1) / ndu(pk + 1, ri);
        d += a(s2, k)*ndu(ri, pk);
      }
      der(ri, k-1) = d;
      /*j = s1;*/ s1 = s2; s2 = j2;
    }
  }


  auto coeff = p;
  for (auto k = 1; k <= num_ders; ++k) {
    der.col(k-1) *= coeff;
    coeff *= p - k;
  }
}

template <typename T> void compute_b_ders_basis(const squanch::NurbsCurve<T>& curve, const unsigned int span, T u, unsigned int num_ders, Eigen::Matrix<T, -1, 1>& r, Eigen::Matrix<T, -1, -1>& der)
{
  ComputeBDers1BasisTemporary<T> temporary(curve);
  compute_b_ders_basis(curve, span, u, temporary, num_ders, r, der);
}

template <typename T> void compute_b_ders1_basis(const squanch::NurbsCurve<T>& curve, const unsigned int span, T u, Eigen::Matrix<T, -1, 1>& r, Eigen::Matrix<T, -1, 1>& der)
{
  ComputeBDers1BasisTemporary<T> temporary(curve);
  compute_b_ders1_basis(curve, span, u, temporary, r, der);
}



template <typename T> struct compute_b_basis_temp {
	compute_b_basis_temp() {}
	compute_b_basis_temp(const squanch::NurbsCurve<T>& curve)
	{
		auto n = curve.p() + 1;
		left.resize(n);
		right.resize(n);
	}
	Eigen::Matrix< T, -1, 1> left;
	Eigen::Matrix< T, -1, 1> right;
};


template <typename T, unsigned int NDIM> struct compute_b_basis_multi_temp {
	std::array<compute_b_basis_temp<T>, NDIM> curve_temps;
	std::array<Eigen::Matrix< T, -1, 1>, NDIM> r_temps;

	template<template<typename, unsigned int> class NurbsType> compute_b_basis_multi_temp(const NurbsType<T, NDIM>& nurbs)
	{
		for (auto i = 0u; i < NDIM; ++i) {
			auto n = nurbs.curves()[i].p() + 1;
			curve_temps[i].left.resize(n);
			curve_temps[i].right.resize(n);
			r_temps[i].resize(n);
		}
	}
};

template <typename T> void compute_b_basis(const squanch::NurbsCurve<T>& curve, const size_t span, T u, Eigen::Matrix<T, -1, 1>& r, compute_b_basis_temp<T>& temp)
{
	auto p = curve.p();
	auto&& knots = curve.knots();

	auto&& left = temp.left;
	auto&& right = temp.right;

	r[0] = (T)1.0;
	for (unsigned int i = 1; i <= p; ++i) {
		T saved = (T)0;
		left[i - 1] = u - knots[span + 1 - i];
		right[i - 1] = knots[span + i] - u;
		for (unsigned int j = 0; j < i; ++j) {
			auto tmp = r[j] / (right[j] + left[i - j - 1]);
			r[j] = saved + right[j] * tmp;
			saved = left[i - j - 1] * tmp;
		}
		r[i] = saved;
	}
}

template <typename T> void compute_b_basis(const squanch::NurbsCurve<T>& curve, const size_t span, T u, Eigen::Matrix<T, -1, 1>& r)
{
	compute_b_basis_temp<T> temp(curve);
	compute_b_basis(curve, span, u, r, temp);
}

template <template<class, unsigned int> class NurbsType, typename T> void compute_b_basis(const NurbsType<T, 1>& nurbs, const unsigned int span, const T u, Eigen::Matrix<T, -1, 1>& r, compute_b_basis_multi_temp<T, 1>& temp)
{
	compute_b_basis(nurbs.curves()[0], span, u, r, temp.curve_temps[0]);
}

template <template<class, unsigned int> class NurbsType, typename T> void compute_b_basis(const NurbsType<T, 1>& nurbs, const unsigned int span, const T u, Eigen::Matrix<T, -1, 1>& r)
{
  compute_b_basis(nurbs.curves()[0], span, u, r);
}

/// <summary>
/// Computes B-Spline Basis function @u for the specified 2D NURBS.
/// </summary>
/// <param name="nurbs">The nurbs.</param>
/// <param name="span">The span of the parametric coordinate u.</param>
/// <param name="u">The parametric coordinate.</param>
/// <param name="r">OUT Vector([p1+1][p2+1]) B-Spline basis functions.</param>
template <template<class, unsigned int> class NurbsType, typename T> void compute_b_basis(const NurbsType<T, 2>& nurbs, const std::array<unsigned int, 2>& span, const std::array<T, 2>& u, 
	Eigen::Matrix<T, -1, 1>& r, compute_b_basis_multi_temp<T, 2>& mtemp)
{
  auto p1 = nurbs.curves()[0].p();
  auto p2 = nurbs.curves()[1].p();

	auto&& rn = mtemp.r_temps[0];
	auto&& rm = mtemp.r_temps[1];

  compute_b_basis(nurbs.curves()[0], span[0], u[0], rn, mtemp.curve_temps[0]);
  compute_b_basis(nurbs.curves()[1], span[1], u[1], rm, mtemp.curve_temps[1]);

  auto num_n = (p1+1)*(p2+1);

  unsigned int c = 0;
  for (unsigned int i = 0; i <= p2; ++i) {
    for (unsigned int j = 0; j <= p1; ++j) {
      r[c]      = rn[j]   * rm[i];
      ++c;
    }
  } 
}

template <template<class, unsigned int> class NurbsType, typename T> void compute_b_basis(const NurbsType<T, 2>& nurbs, const std::array<unsigned int, 2>& span, const std::array<T, 2>& u,
	Eigen::Matrix<T, -1, 1>& r)
{
	compute_b_basis_multi_temp<T, 2> mtemp(nurbs);
	compute_b_basis(nurbs, span, u, r, mtemp);
}

/// <summary>
/// Computes B-Spline Basis function @u for the specified 3D NURBS.
/// </summary>
/// <param name="nurbs">The nurbs.</param>
/// <param name="span">The span of the parametric coordinate u.</param>
/// <param name="u">The parametric coordinate.</param>
/// <param name="r">OUT Vector([p1+1][p2+1][p3+1]) B-Spline basis functions.</param>>
template <template<class, unsigned int> class NurbsType, typename T> void compute_b_basis(const NurbsType<T, 3>& nurbs, const std::array<unsigned int, 3>& span, const std::array<T, 3>& u, 
	Eigen::Matrix<T, -1, 1>& r, compute_b_basis_multi_temp<T, 3>& mtemp)
{
  auto p1 = nurbs.curves()[0].p();
  auto p2 = nurbs.curves()[1].p();
  auto p3 = nurbs.curves()[2].p();

	auto&& n = mtemp.r_temps[0];
	auto&& m = mtemp.r_temps[1];
	auto&& l = mtemp.r_temps[2];

  compute_b_basis(nurbs.curves()[0], span[0], u[0], n, mtemp.curve_temps[0]);
  compute_b_basis(nurbs.curves()[1], span[1], u[1], m, mtemp.curve_temps[1]);
  compute_b_basis(nurbs.curves()[2], span[2], u[2], l, mtemp.curve_temps[2]);

  unsigned int c = 0;
  for (unsigned int k = 0; k <= p3; ++k) {
    for (unsigned int i = 0; i <= p2; ++i) {
      for (unsigned int j = 0; j <= p1; ++j) {
        r[c] = n[j] * m[i] * l[k];
        ++c;
      }
    }
  }
}

template <template<class, unsigned int> class NurbsType, typename T> void compute_b_basis(const NurbsType<T, 3>& nurbs, const std::array<unsigned int, 3>& span, const std::array<T, 3>& u,
	Eigen::Matrix<T, -1, 1>& r)
{
	compute_b_basis_multi_temp<T, 3> mtemp(nurbs);
	compute_b_basis(nurbs, span, u, r, mtemp);
}

/// <summary>
/// Computes NURBS basis functions @u for the specified 1D NURBS.
/// </summary>
/// <param name="nurbs">The nurbs.</param>
/// <param name="span">The span of the parametric coordinate u.</param>
/// <param name="u">The parametric coordinate</param>
/// <param name="cp">The control point map.</param>
/// <param name="r">OUT Vector(p+1) NURBS basis functions</param>
template <template<class, unsigned int> class NurbsType, typename T> void compute_nurbs_basis(const NurbsType<T, 1>& nurbs, const unsigned int span, const T u, const Eigen::Matrix<T, -1, 1>& weights, Eigen::Matrix<T, -1, 1>& r)
{
  auto p1 = nurbs.curves()[0].p();

  compute_b_basis(nurbs.curves()[0], span, u, r);
  r = r.cwiseProduct(weights);
  T w = r.sum();
  T winv = (T)1 / w;
  r *= winv;
}
/// <summary>
/// Computes NURBS basis functions @u for the specified 2D NURBS.
/// </summary>
/// <param name="nurbs">The nurbs.</param>
/// <param name="span">The span of the parametric coordinate u.</param>
/// <param name="u">The parametric coordinate</param>
/// <param name="cp">The control point map.</param>
/// <param name="r">OUT Vector([p1+1][p2+1]) NURBS basis functions</param>
template <template<class, unsigned int> class NurbsType, typename T> void compute_nurbs_basis(const NurbsType<T, 2>& nurbs, const std::array<unsigned int, 2>& span, const std::array<T, 2>& u, const Eigen::Matrix<T, -1, 1>& weights, Eigen::Matrix<T, -1, 1>& r)
{
  auto p1 = nurbs.curves()[0].p();
  auto p2 = nurbs.curves()[1].p();

  Eigen::Matrix<T, -1,  1> rn(p1+1);
  Eigen::Matrix<T, -1,  1> rm(p2+1);
  compute_b_basis(nurbs.curves()[0], span[0], u[0], rn);
  compute_b_basis(nurbs.curves()[1], span[1], u[1], rm);

  auto num_n = (p1+1)*(p2+1);

  unsigned int c = 0;
  for (unsigned int i = 0; i <= p2; ++i) {
    for (unsigned int j = 0; j <= p1; ++j) {
      r[c]      = rn[j]   * rm[i]   * weights[c];
      ++c;
    }
  } 

  T w = r.sum();
  T winv = (T)1 / w;
  r *= winv;
}
/// <summary>
/// Computes NURBS basis functions @u for the specified 3D NURBS.
/// </summary>
/// <param name="nurbs">The nurbs.</param>
/// <param name="span">The span of the parametric coordinate u.</param>
/// <param name="u">The parametric coordinate</param>
/// <param name="cp">The control point map.</param>
/// <param name="r">OUT Vector([p1+1][p2+1][p3+1]) NURBS basis functions</param>
template <template<class, unsigned int> class NurbsType, typename T> void compute_nurbs_basis(const NurbsType<T, 3>& nurbs, const std::array<unsigned int, 3>& span, const std::array<T, 3>& u, const Eigen::Matrix<T, -1, 1>& weights, Eigen::Matrix<T, -1, 1>& r)
{
  auto p1 = nurbs.curves()[0].p();
  auto p2 = nurbs.curves()[1].p();
  auto p3 = nurbs.curves()[2].p();

  Eigen::Matrix<T, -1, 1> n(p1+1);
  Eigen::Matrix<T, -1, 1> m(p2+1);
  Eigen::Matrix<T, -1, 1> l(p3+1);
  compute_b_basis(nurbs.curves()[0], span[0], u[0], n);
  compute_b_basis(nurbs.curves()[1], span[1], u[1], m);
  compute_b_basis(nurbs.curves()[2], span[2], u[2], l);

  unsigned int c = 0;
  for (unsigned int k = 0; k <= p3; ++k) {
    for (unsigned int i = 0; i <= p2; ++i) {
      for (unsigned int j = 0; j <= p1; ++j) {
        r[c] = n[j] * m[i] * l[k] * weights[c];
        ++c;
      }
    }
  }

  T w = r.sum();
  T winv = (T)1 / w;
  r   *= winv;
}

template <typename T> class ComputeNurbsDers1BasisTemporary
{
private:
  template <typename T> struct ComputeNurbsDers1BasisTemporaryPerCurve
  {
    Eigen::Matrix<T, -1, 1> basis;
    Eigen::Matrix<T, -1, 1> derivative;
    ComputeBDers1BasisTemporary<T> b_ders1basis_temporary;
  };
public:
  template <unsigned int NDIM> ComputeNurbsDers1BasisTemporary(const squanch::Nurbs<T, NDIM>& nurbs)
  {
    curve_temporary.resize(NDIM);
    for (auto i = 0u; i < NDIM; ++i) {
      auto p = nurbs.curves()[i].p();
      curve_temporary[i].basis.resize(p + 1);
      curve_temporary[i].derivative.resize(p + 1);
      curve_temporary[i].b_ders1basis_temporary.reset(nurbs.curves()[i]);
    }
  }

  std::vector<ComputeNurbsDers1BasisTemporaryPerCurve<T>> curve_temporary;
};


/// <summary>
/// Computes NURBS basis functions and their 1st derivatives @u for the specified 1D NURBS.
/// </summary>
/// <param name="nurbs">The nurbs.</param>
/// <param name="span">The span of the parametric coordinate u.</param>
/// <param name="u">The parametric coordinate</param>
/// <param name="cp">The control point map.</param>
/// <param name="r">OUT Vector(p+1) NURBS basis functions</param>>
/// <param name="der">OUT Matrix(p+1,1)The 1st derivatives of the NURBS basis functions.</param>
template <typename T> void compute_nurbs_ders1_basis(const squanch::Nurbs<T, 1>& nurbs, const unsigned int span, const T u, const Eigen::Matrix<T, -1, 1>& weights,
                                                     Eigen::Matrix<T, -1, 1>& r, Eigen::Matrix<T, -1, 1>& der)
{
  /* uses Equation 2.31 from IGA book to compute rational basis values and 1st derivatives. */

  compute_b_ders1_basis(nurbs.curves()[0], span, u, r, der);

  r = r.cwiseProduct(weights);
  der.col(0) = der.col(0).cwiseProduct(weights);
  
  T w  = r.sum();
  T winv = (T)1 / w;
  T wd = der.col(0).sum();

  der.col(0) -= r*wd*winv;
  r *= winv;
  der *= winv;
}

/// <summary>
/// Computes NURBS basis functions and their 1st derivatives @u for the specified 2D NURBS.
/// </summary>
/// <param name="nurbs">The nurbs.</param>
/// <param name="span">The span of the parametric coordinate u.</param>
/// <param name="u">The parametric coordinate</param>
/// <param name="cp">The control point map.</param>
/// <param name="r">OUT Vector([p1+1][p2+1]) NURBS basis functions</param>>
/// <param name="der">OUT Matrix([p1+1][p2+1],2)The 1st derivatives of the NURBS basis functions.</param>>
template <typename T> void compute_nurbs_ders1_basis(const squanch::Nurbs<T, 2>& nurbs, const std::array<unsigned int, 2>& span, const std::array<T, 2>& u, const Eigen::Matrix<T, -1, 1>& weights,
                                                     ComputeNurbsDers1BasisTemporary<T>& temporary,
                                                     Eigen::Matrix<T, -1, 1>& r, Eigen::Matrix<T, -1, 2>& der)
{
  auto p1 = nurbs.curves()[0].p();
  auto p2 = nurbs.curves()[1].p();

  auto&& rn = temporary.curve_temporary[0].basis;
  auto&& dn = temporary.curve_temporary[0].derivative;
  auto&& rm = temporary.curve_temporary[1].basis;
  auto&& dm = temporary.curve_temporary[1].derivative;

  compute_b_ders1_basis(nurbs.curves()[0], span[0], u[0], temporary.curve_temporary[0].b_ders1basis_temporary, rn, dn);
  compute_b_ders1_basis(nurbs.curves()[1], span[1], u[1], temporary.curve_temporary[1].b_ders1basis_temporary, rm, dm);

  auto num_n = (p1+1)*(p2+1);

  unsigned int c = 0;
  for (unsigned int i = 0; i <= p2; ++i) {
    for (unsigned int j = 0; j <= p1; ++j) {
      r[c]      = rn[j]   * rm[i]   * weights[c];
      der.col(0)[c] = dn(j,0) * rm[i]   * weights[c];
      der.col(1)[c] = rn[j]   * dm(i,0) * weights[c];
      ++c;
    }
  } 


  T w = r.sum();
  T winv = (T)1 / w;
  T wn = der.col(0).sum();
  T wm = der.col(1).sum();

  der.col(0) -= r.transpose() * wn * winv;
  der.col(1) -= r.transpose() * wm * winv;
  r *= winv;
  der *= winv;
}

/// <summary>
/// Computes NURBS basis functions and their 1st derivatives @u for the specified 2D NURBS.
/// </summary>
/// <param name="nurbs">The nurbs.</param>
/// <param name="span">The span of the parametric coordinate u.</param>
/// <param name="u">The parametric coordinate</param>
/// <param name="cp">The control point map.</param>
/// <param name="r">OUT Vector([p1+1][p2+1]) NURBS basis functions</param>>
/// <param name="der">OUT Matrix([p1+1][p2+1],2)The 1st derivatives of the NURBS basis functions.</param>>
template <typename T> void compute_nurbs_ders1_basis(const squanch::Nurbs<T, 2>& nurbs, const std::array<unsigned int, 2>& span, const std::array<T, 2>& u, const Eigen::Matrix<T, -1, 1>& weights,
                                                     Eigen::Matrix<T, -1, 1>& r, Eigen::Matrix<T, -1, 2>& der)
{
  ComputeNurbsDers1BasisTemporary<T> temporary(nurbs);
  compute_nurbs_ders1_basis(nurbs, span, u, weights, temporary, r, der);
}

/// <summary>
/// Computes NURBS basis functions and their 1st derivatives @u for the specified 3D NURBS.
/// </summary>
/// <param name="nurbs">The nurbs.</param>
/// <param name="span">The span of the parametric coordinate u.</param>
/// <param name="u">The parametric coordinate</param>
/// <param name="cp">The control point map.</param>
/// <param name="r">OUT Vector([p1+1][p2+1][p3+1]) NURBS basis functions</param>>
/// <param name="der">OUT Matrix([p1+1][p2+1][p3+1],3)The 1st derivatives of the NURBS basis functions.</param>>
template <typename T> void compute_nurbs_ders1_basis(const squanch::Nurbs<T, 3>& nurbs, const std::array<unsigned int, 3>& span, const std::array<T, 3>& u, const Eigen::Matrix<T, -1, 1>& weights,
                                                     ComputeNurbsDers1BasisTemporary<T>& temporary, Eigen::Matrix<T, -1, 1>& r, Eigen::Matrix<T, -1, 3>& der)
{
  auto p1 = nurbs.curves()[0].p();
  auto p2 = nurbs.curves()[1].p();
  auto p3 = nurbs.curves()[2].p();

  auto&& n = temporary.curve_temporary[0].basis;
  auto&& dn = temporary.curve_temporary[0].derivative;
  auto&& m = temporary.curve_temporary[1].basis;
  auto&& dm = temporary.curve_temporary[1].derivative;
  auto&& l = temporary.curve_temporary[2].basis;
  auto&& dl = temporary.curve_temporary[2].derivative;

  compute_b_ders1_basis(nurbs.curves()[0], span[0], u[0], temporary.curve_temporary[0].b_ders1basis_temporary, n, dn);
  compute_b_ders1_basis(nurbs.curves()[1], span[1], u[1], temporary.curve_temporary[1].b_ders1basis_temporary, m, dm);
  compute_b_ders1_basis(nurbs.curves()[2], span[2], u[2], temporary.curve_temporary[2].b_ders1basis_temporary, l, dl);

  unsigned int c = 0;
  for (unsigned int k = 0; k <= p3; ++k) {
    for (unsigned int i = 0; i <= p2; ++i) {
      for (unsigned int j = 0; j <= p1; ++j) {
        r[c] = n[j] * m[i] * l[k] * weights[c];
        der.col(0)[c] = dn(j,0) * m[i]    * l[k]    * weights[c];
        der.col(1)[c] = n[j]    * dm(i,0) * l[k]    * weights[c];
        der.col(2)[c] = n[j]    * m[i]    * dl(k,0) * weights[c];
        ++c;
      }
    }
  }

  T w = r.sum();
  T winv = (T)1 / w;
  T wn = der.col(0).sum();
  T wm = der.col(1).sum();
  T wl = der.col(2).sum();

  der.col(0) -= r.transpose() * wn * winv;
  der.col(1) -= r.transpose() * wm * winv;
  der.col(2) -= r.transpose() * wl * winv;
  r   *= winv;
  der *= winv;

}

/// <summary>
/// Computes NURBS basis functions and their 1st derivatives @u for the specified 3D NURBS.
/// </summary>
/// <param name="nurbs">The nurbs.</param>
/// <param name="span">The span of the parametric coordinate u.</param>
/// <param name="u">The parametric coordinate</param>
/// <param name="cp">The control point map.</param>
/// <param name="r">OUT Vector([p1+1][p2+1][p3+1]) NURBS basis functions</param>>
/// <param name="der">OUT Matrix([p1+1][p2+1][p3+1],3)The 1st derivatives of the NURBS basis functions.</param>>
template <typename T> void compute_nurbs_ders1_basis(const squanch::Nurbs<T, 3>& nurbs, const std::array<unsigned int, 3>& span, const std::array<T, 3>& u, const Eigen::Matrix<T, -1, 1>& weights,
                                                     Eigen::Matrix<T, -1, 1>& r, Eigen::Matrix<T, -1, 3>& der)
{
  ComputeNurbsDers1BasisTemporary<T> temporary(nurbs);
  compute_nurbs_ders1_basis(nurbs, span, u, weights, temporary, r, der);
}

}