#include <boost/test/unit_test.hpp>
#include <squanch/Model.h>
#include <squanch/Line.h>
#include <squanch/Refinement.h>
#include <squanch/Extrude.h>
#include "igafem.h"

template <typename T, unsigned int NDIM> struct ParametricPointComputation {};

template <typename T> struct ParametricPointComputation<T, 1> {
	typedef T ValueType;

	template <typename FunType> void compute_at_points(std::array<T, 1> dt, FunType fun) {
		unsigned int n = (1 / dt[0]) + 1;
		for (auto i = 0u; i < n; ++i) {
			auto u = i * dt[0];
			fun(u);
		}
	}
};

template <typename T> struct ParametricPointComputation<T, 2> {
	typedef std::array<T, 2> ValueType;

	template <typename FunType> void compute_at_points(std::array<T, 2> dt, FunType fun) {
		ValueType u;
		unsigned int n = (1 / dt[0]) + 1;
		unsigned int n2 = (1 / dt[1]) + 1;
		for (auto i = 0u; i < n; ++i) {
			u[0] = i * dt[0];
			for (auto j = 0u; j < n2; ++j) {
				u[1] = j * dt[1];
				fun(u);
			}
		}
	}
};

template <typename T> struct ParametricPointComputation<T, 3> {
	typedef std::array<T, 3> ValueType;

	template <typename FunType> void compute_at_points(std::array<T, 3> dt, FunType fun) {
		ValueType u;
		unsigned int n = (1 / dt[0]) + 1;
		unsigned int n2 = (1 / dt[1]) + 1;
		unsigned int n3 = (1 / dt[2]) + 1;
		for (auto i = 0u; i < n; ++i) {
			u[0] = i * dt[0];
			for (auto j = 0u; j < n2; ++j) {
				u[1] = j * dt[1];
				for (auto k = 0u; k < n3; ++k) {
					u[2] = j * dt[2];
					fun(u);
				}
			}
		}
	}
};

auto make_vector = [](const Eigen::VectorXd& eig) {
	std::vector<double> vec(eig.rows());
	for (auto i = 0u; i < vec.size(); ++i) vec[i] = eig[i];
	return vec;
};

void igafem_compute(double u, const squanch::Nurbs<double, 1>& nurbs, const Eigen::VectorXd& weights, Eigen::VectorXd& r, Eigen::VectorXd& der)
{
	igafem::nurbs1DBasisDers(u, nurbs.curves()[0].p(), nurbs.curves()[0].knots(), make_vector(weights), r, der);
}

void igafem_compute(std::array<double, 2> u, const squanch::Nurbs<double, 2>& nurbs, const Eigen::VectorXd& weights, Eigen::VectorXd& r, Eigen::Matrix<double, -1, 2>& der)
{
	auto ncp = nurbs.num_cp();
	Eigen::VectorXd rx(ncp);
	Eigen::VectorXd ry(ncp);

	igafem::nurbs2DBasisDers(u[0], u[1], nurbs.curves()[0].p(), nurbs.curves()[1].p(), nurbs.curves()[0].knots(), nurbs.curves()[1].knots(), make_vector(weights), r, rx, ry);
	der.col(0) = rx;
	der.col(1) = ry;
}

void igafem_compute(std::array<double, 3> u, const squanch::Nurbs<double, 3>& nurbs, const Eigen::VectorXd& weights, Eigen::VectorXd& r, Eigen::Matrix<double, -1, 3>& der)
{
	auto ncp = nurbs.num_cp();
	Eigen::VectorXd rx(ncp);
	Eigen::VectorXd ry(ncp);
	Eigen::VectorXd rz(ncp);

	igafem::nurbs3DBasisDers(u[0], u[1], u[2], nurbs.curves()[0].p(), nurbs.curves()[1].p(), nurbs.curves()[2].p(),
		nurbs.curves()[0].knots(), nurbs.curves()[1].knots(), nurbs.curves()[2].knots(), make_vector(weights), r, rx, ry, rz);

	der.col(0) = rx;
	der.col(1) = ry;
	der.col(2) = rz;
}

template <unsigned int NDIM> void check_matrix(const Eigen::Matrix<double, -1, NDIM>& m1, const Eigen::Matrix<double, -1, NDIM>& m2)
{
	for (auto i = 0u; i < m1.rows(); ++i) {
		for (auto j = 0u; j < m1.cols(); ++j) {
			auto val = std::abs(m1(i, j) - m2(i, j));
			bool is_close = val < 1e-8;
			BOOST_TEST(is_close);
		}
	}
}

template <unsigned int NDIM> void check_iga_computations(const squanch::Model<double>& model, const squanch::Nurbs<double, NDIM>& nurbs, std::array<double, NDIM> dt) 
{
	auto ncp = nurbs.num_cp();
	Eigen::VectorXd weights(ncp);
	Eigen::VectorXd r(ncp);
	Eigen::VectorXd r2(ncp);
	Eigen::Matrix<double, -1, NDIM> ders(ncp, NDIM);
	Eigen::Matrix<double, -1, NDIM> ders2(ncp, NDIM);

	using PpcType = ParametricPointComputation<double, NDIM>;
	PpcType ppc;

	ppc.compute_at_points(dt, [&](PpcType::ValueType u){
		auto span = squanch::find_span(nurbs, u);
		squanch::gather_weights(model, nurbs, span, weights);

		igafem_compute(u, nurbs, weights, r2, ders2);
		squanch::compute_nurbs_ders1_basis(nurbs, span, u, weights, r, ders);

		//std::cout << r.transpose() << "\n\n";
		//std::cout << r2.transpose() << "\n\n";
		//std::cout << ders.transpose() << "\n\n";
		//std::cout << ders2.transpose() << "\n\n";

		check_matrix<1>(r, r2);
		check_matrix<NDIM>(ders, ders2);
	});
}

squanch::Nurbs<double, 1>& create_curve(squanch::Model<double>& model)
{
	squanch::Line<double, Eigen::Vector3d> line;
	line.p1 = Eigen::Vector3d(0, 0, 0);
	line.p2 = Eigen::Vector3d(10, 0, 0);
	auto&& curve = line.to_nurbs(model);
	squanch::p_refine(model, curve, 2);

	auto set_node = [&model, &curve](int id, double x, double y, double z, double w){
		auto&& p = model.get_point(id);
		p[0] = x;
		p[1] = y;
		p[2] = z;
		p[3] = w;
	};

	set_node(curve.cpi()(0), 5, 0, 0, 1);
	set_node(curve.cpi()(1), 5.5, 0.3, 0, 1.2);
	set_node(curve.cpi()(2), 5.7, 0.6, 0, 1.5);
	set_node(curve.cpi()(3), 5, 1, 1, 1);

	return curve;
}


BOOST_AUTO_TEST_SUITE(dornisch)

BOOST_AUTO_TEST_CASE(curve)
{
  squanch::Model<double> model;

	auto&& curve = create_curve(model);
	squanch::p_refine(model, curve,  3 );
	squanch::h_refine_elements(model, curve, 10);

	check_iga_computations<1>(model, curve, { { 0.05 } });

}

BOOST_AUTO_TEST_CASE(surface)
{
	squanch::Model<double> model;

	auto&& curve = create_curve(model);
	Eigen::Vector4d ex(0, 1, 0, 0);
	auto&& surface = squanch::extrude(model, curve, ex);

	squanch::p_refine(model, surface, { { 3, 0 } });
	squanch::h_refine_elements(model, surface, 10,  1);

	check_iga_computations<2>(model, surface, { { 0.05, 0.05 } });


}

BOOST_AUTO_TEST_CASE(volume)
{
	squanch::Model<double> model;

	auto&& curve = create_curve(model);
	Eigen::Vector4d ex(0, 1, 0, 0);
	auto&& surface = squanch::extrude(model, curve, ex);
	Eigen::Vector4d ex2(0, 0, 1, 0);
	auto&& volume = squanch::extrude(model, surface, ex);

  squanch::p_refine(model, volume, { {3, 0, 0} });
	squanch::h_refine_elements(model, volume, 10, 1, 1);

	const auto dt = .05;
	check_iga_computations<3>(model, volume, { { dt, dt, dt } });
}

BOOST_AUTO_TEST_SUITE_END()