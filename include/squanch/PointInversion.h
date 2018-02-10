#pragma once
#include <iostream>
#include <squanch/CurveAlgorithms.h>
#include <squanch/SurfaceAlgorithms.h>
#include <squanch/IgaAlgorithms.h>
#include <squanch/Array2D.h>

namespace squanch
{

  namespace detail_pi {

    template <typename T> void enforce_point_condition(T& val, bool is_closed)
    {
      if (is_closed) {
        if (val < 0) val += 1;
        if (val > 1) val -= 1;
      }
      else {
        if (val < 0) val = 0;
        if (val > 1) val = 1;
      }
    };

    template <typename T, typename PointType> T compute_squared_norm(const PointType& p){
      return p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
    };

    template <typename T, typename PointType> T compute_norm(const PointType& p){
      return sqrt(compute_squared_norm<T>(p));
    };

    template <typename T> T find_smallest_span(const NurbsCurve<T>& curve)
    {
      auto&& knots = curve.knots();
      const auto e = (T)0.0000001;
      T min_span_length = std::numeric_limits<T>::max();

      for (auto i = 1; i < knots.size(); ++i) {
        auto span_length = knots[i] - knots[i - 1];
        if (span_length > e) {
          if (span_length < min_span_length) min_span_length = span_length;
        }
      }

      return min_span_length;
    }

    template <template<class, unsigned int> class NurbsType, typename T>  std::vector<T> find_available_points(const NurbsType<T, 1>& spline, T guess)
    {
      const auto e = 0.0000001;
			const auto num_parts_per_span = 20;
      std::vector<T> points;

      auto&& curve = spline.curves()[0];
      auto&& knots = curve.knots();

      T min_span_length = find_smallest_span(curve);
			auto du = min_span_length / num_parts_per_span;
      unsigned int n = (T)1.0 / du + 1;
      for (auto i = 0u; i < n; ++i) {
        points.push_back(du * i);
      }

			auto guess_delta = min_span_length / (num_parts_per_span * 2);
			for (auto i = 0; i < num_parts_per_span; ++i) {
				auto p = guess - (num_parts_per_span / 2) * guess_delta + i * guess_delta;
				if (p < 0) p = 0;
				if (p > 1) p = 1;

				points.push_back(p);
			}

      return points;
    }

    template <template<class, unsigned int> class NurbsType, typename T>  std::vector<std::array<T, 2>> find_available_points(const NurbsType<T, 2>& spline, const std::array<T, 2>& guess)
    {
      const auto e = 0.0000001;
      const auto ndim = 2u;
			const auto num_parts_per_span = 20;

      std::vector<std::array<T, ndim>> points;
      std::array<T, ndim> du;
			for (auto i = 0u; i < ndim; ++i) du[i] = find_smallest_span(spline.curves()[i]) / num_parts_per_span;

      std::array<unsigned int, ndim> n;
      for (auto i = 0u; i < ndim; ++i) n[i] = (T)1 / du[i] + 1;

      for (auto i = 0u; i < n[0]; ++i) {
        auto dx = du[0] * i;
        for (auto j = 0u; j < n[1]; ++j) {
          points.push_back({ { dx, du[1] * j } });
        }
      }

			auto guess_deltax = du[0] / (num_parts_per_span * 2);
			auto guess_deltay = du[1] / (num_parts_per_span * 2);

			for (auto i = 0; i < num_parts_per_span; ++i) {
				auto px = guess[0] - (num_parts_per_span / 2) * guess_deltax + i * guess_deltax;
				if (px < 0) px = 0;
				if (px > 1) px = 1;
				for (auto j = 0; j < num_parts_per_span; ++j) {
					auto py = guess[1] - (num_parts_per_span / 2) * guess_deltay + i * guess_deltay;
					if (py < 0) py = 0;
					if (py > 1) py = 1;
					points.push_back({ {px, py} });
				}
			}

      return points;
    }

    template <template<class, unsigned int> class NurbsType, typename T, unsigned int NDIM> typename NurbsType<T, NDIM>::ParametricPointType 
      find_starting_point(const squanch::Model<T>& model, const NurbsType<T, NDIM>& nurbs, typename const squanch::PointContainerTypes<T>::point_type& point, const std::vector<typename NurbsType<T, NDIM>::ParametricPointType>& available_points)
    {
      NurbsType<T, NDIM>::ParametricPointType closest;
      if (available_points.empty()) return NurbsType<T, NDIM>::ParametricPointType();

      closest = available_points[0];
      T closest_norm = std::numeric_limits<T>::max();

      squanch::PointContainerTypes<T>::point_type point_nh = point;
      auto numnodes = num_nodes_per_span(nurbs.curves());

      Eigen::Matrix<T, -1, 1> n(numnodes);
      PointContainerTypes<T>::point_type p_new;
			compute_b_basis_multi_temp<T, NDIM> mtemp(nurbs);

      for (auto u : available_points) {
        auto span = find_span(nurbs, u);
        compute_b_basis(nurbs, span, u, n, mtemp);
        p_new = compute_nurbs_point(model, nurbs, u, span, n);
        auto norm = compute_squared_norm<T>(p_new - point_nh);
        if (norm < closest_norm) {
          closest = u;
          closest_norm = norm;
        }
      }

      return closest;
    }

  }




  template <template<class, unsigned int> class NurbsType, typename T> T inverse_point_nr(const squanch::Model<T>& model, const NurbsType<T, 1>& nurbs, T guess, typename const squanch::PointContainerTypes<T>::point_type& point)
  {
    const T small_number = (T) 0.000001;
    const auto num_max_iter = 100;

		T u0 = guess;

    auto&& curve = nurbs.curves()[0];
    const auto p = curve.p();
    auto&& knots = curve.knots();

    bool is_converged = false;
    auto iter = 0;

    Eigen::Matrix<T, -1, 1> n(p + 1);
    Eigen::Matrix<T, -1, -1> der(p + 1, 2);
    Eigen::Matrix<T, 4, -1> Bder(4, 3);
    Eigen::Matrix<T, 4, -1> Cder(4, 3);
    squanch::PointContainerTypes<T>::point_type C;

    auto u = u0;
    while (!is_converged && iter < num_max_iter) {
      auto span = find_span(knots, p, u);

      compute_b_ders_basis(curve, span, u, 2, n, der);
      Bder = compute_bspline_derivatives(model, nurbs, u, span, n, der);
      compute_nurbs_derivatives(Bder, Cder);
			C = Cder.col(0);

			auto err = detail_pi::compute_norm<T>(C - point);
      if (err <= small_number) return u;
			err = abs(Cder.col(1).dot(C - point)) / (detail_pi::compute_norm<T>(Cder.col(1)) * detail_pi::compute_norm<T>(Cder.col(1) - point));
			if (err <= small_number) {
				std::cout << "Warning, Point Inversion only converged with " << detail_pi::compute_norm<T>(C - point) << "\n";
				return u;
			}

      auto der1_norm_square = detail_pi::compute_squared_norm<T>(Cder.col(1));
      u = u - Cder.col(1).dot(C - point) / (Cder.col(2).dot(C - point) + der1_norm_square);

      detail_pi::enforce_point_condition(u, curve.is_closed());
    }
		std::cerr << "Warning Point reversal did not converge \n";
    return u;
  }

	template <template<class, unsigned int> class NurbsType, typename T> T inverse_point_nr(const squanch::Model<T>& model, const NurbsType<T, 1>& nurbs, typename const squanch::PointContainerTypes<T>::point_type& point)
	{
		auto&& available_points = detail_pi::find_available_points(nurbs, 0.0);
		T u0 = detail_pi::find_starting_point(model, nurbs, point, available_points);

		return inverse_point_nr(model, nurbs, u0, point);
	}


  template <template<class, unsigned int> class NurbsType, typename T> std::array<T, 2> inverse_point_nr(const squanch::Model<T>& model, const NurbsType<T, 2>& nurbs, const std::array<T, 2>& guess,
		typename const squanch::PointContainerTypes<T>::point_type& point)
  {
    auto invert_matrix = [](Eigen::Matrix<T, 2, 2>& mat){
      auto a = mat(0, 0);
      auto b = mat(0, 1);
      auto c = mat(1, 0);
      auto d = mat(1, 1);
      auto coeff = 1 / (a*d - b*c);
      
      mat(0, 0) = coeff *  d;
      mat(0, 1) = coeff * -b;
      mat(1, 0) = coeff * -c;
      mat(1, 1) = coeff *  a;
    };


    const T small_number = (T) 0.000001;
    const auto num_max_iter = 1000;

		auto u0 = guess;
    auto iter = 0;
    
    PointArray2D<T> bspline_SKL(2, 2);
    PointArray2D<T> SKL(2, 2);

    squanch::PointContainerTypes<T>::point_type r;
    
    Eigen::Matrix<T, 2, 1> k;
    Eigen::Matrix<T, 2, 1> d;
    Eigen::Matrix<T, 2, 2> j;
    Eigen::Matrix<T, 2, 2> ji;

    // guarantee no convergence at the beginning
    d[0] = std::numeric_limits<double>::max();
    d[1] = std::numeric_limits<double>::max();

    auto u = u0;
    while (iter < num_max_iter) {
      auto span = find_span(nurbs, u);

      compute_bspline_surface_derivatives(model, nurbs, 1, u, bspline_SKL);
      compute_nurbs_surface_derivatives(model, bspline_SKL, SKL);

      auto&& S   = SKL(0, 0);
      auto&& Su  = SKL(1, 0);
      auto&& Sv  = SKL(0, 1);
      auto&& Suv = SKL(1, 1);
      
      r = S - point;

			//std::cout << r.head<3>().norm() << " and " << (d[0] * Su + d[1] * Sv).norm() <<  "\n";

      // check convergence
      if (detail_pi::compute_norm<double>(r) <= small_number) return u;
			auto conv_norm = detail_pi::compute_norm<double>(d[0] * Su + d[1] * Sv);
			if (conv_norm <= small_number) {
				std::cout << "Warning, Point Inversion only converged with " << detail_pi::compute_norm<double>(r) << "\n";
				return u;
			}
      //if (std::abs(Su.dot(r)) / (Su.norm() * r.norm()) <= small_number) return u;
      //if (std::abs(Sv.dot(r)) / (Sv.norm() * r.norm()) <= small_number) return u;

      // compute new point
      k[0] = -r.dot(Su);
      k[1] = -r.dot(Sv);

      j(0, 0) = Su.squaredNorm() + r.dot(Su);
      j(0, 1) = Su.dot(Sv) + r.dot(Suv);
      j(1, 0) = Su.dot(Sv) + r.dot(Suv);
      j(1, 1) = Sv.squaredNorm() + r.dot(Sv);

      invert_matrix(j);
      d.noalias() = j * k;

      u[0] += d[0];
      u[1] += d[1];

      // check new point criteria
      detail_pi::enforce_point_condition(u[0], nurbs.curves()[0].is_closed());
      detail_pi::enforce_point_condition(u[1], nurbs.curves()[1].is_closed());
      ++iter;
    }

		std::cerr << "Warning Point reversal did not converge \n";
    return u;
  }

	template <template<class, unsigned int> class NurbsType, typename T> std::array<T, 2> inverse_point_nr(const squanch::Model<T>& model, const NurbsType<T, 2>& nurbs, typename const squanch::PointContainerTypes<T>::point_type& point)
	{
		std::array<double, 2> guess = { { 0.0, 0.0 } };
		auto&& available_points = detail_pi::find_available_points(nurbs, guess);
		auto u0 = detail_pi::find_starting_point(model, nurbs, point, available_points);

		return inverse_point_nr(model, nurbs, u0, point);
	}



}