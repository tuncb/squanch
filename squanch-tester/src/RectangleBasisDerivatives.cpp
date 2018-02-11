#include "catch.hpp"
#include <squanch/Model.h>
#include <squanch/Nurbs.h>
#include <squanch/Line.h>
#include <squanch/Extrude.h>
#include <squanch/Refinement.h>
#include <squanch/IgaAlgorithms.h>
#include "DebugTools.h"
#include "igafem.h"


inline squanch::Nurbs<double, 2>& create_rectangle
(squanch::Model<double>& model, unsigned int lx, unsigned int ly, unsigned int hx, unsigned int hy, unsigned int px, unsigned int py)
{
	using namespace squanch;
	using namespace Eigen;

	Line<double, Vector3d> line;
	line.p1 = Vector3d(0, 0, 0);
	line.p2 = Vector3d(hx, 0, 0);

	auto&& nline = line.to_nurbs(model);
	Eigen::Vector4d vec(0, hy, 0, 0);
	auto&& rect = squanch::extrude(model, nline, vec);

	std::array<unsigned int, 2> p_arr = { { px, py } };
	squanch::p_refine(model, rect, p_arr);

	squanch::h_refine_elements(model, rect, hx, hy);

	return rect;
}

inline void check_nurbs(const squanch::Model<double>& model, const squanch::Nurbs<double, 2>& nurbs, unsigned int nx, unsigned int ny)
{
	using namespace squanch;

	double dx = 1.0 / nx;
	double dy = 1.0 / ny;

	double x = 0;
	double y = 0;

	auto num_cp = nurbs.num_cp();

	std::vector<double> weight_vec(num_cp);
	Eigen::VectorXd weights(num_cp);
	Eigen::VectorXd R(num_cp);
	Eigen::VectorXd du(num_cp);
	Eigen::VectorXd dv(num_cp);

	Eigen::VectorXd Ru_ozp(num_cp);
	Eigen::VectorXd Rv_ozp(num_cp);
	Eigen::VectorXd Ru(num_cp);
	Eigen::VectorXd Rv(num_cp);

	Eigen::VectorXd r_ozp(num_cp);
	Eigen::Matrix<double, -1, 2> der_ozp(num_cp, 2);

	std::array<unsigned int, 2> span;
	std::array<double, 2> u;

	auto equate_weights = [](const Eigen::VectorXd& w1, std::vector<double>& w2) {
		for (auto i = 0; i < w1.rows(); ++i) w2[i] = w1[i];
	};

	for (auto i = 0u; i <= ny; ++i) {
		y = i * dy;
		if (y > 1.0) y = 1.0;
		span[1] = find_span(nurbs.curves()[1].knots(), nurbs.curves()[1].p(), y);
		u[1] = y;
		for (auto j = 0u; j <= nx; ++j) {
			x = dx * j;
			if (x > 1.0) x = 1.0;
			u[0] = x;
			span[0] = find_span(nurbs.curves()[0].knots(), nurbs.curves()[0].p(), x);
			gather_weights(model, nurbs, span, weights);

			equate_weights(weights, weight_vec);

			/*    igafem::BasisFuns(span[0], u[0], nurbs.curves()[0].p(), nurbs.curves()[0].knots(), Ru);
					igafem::BasisFuns(span[1], u[1], nurbs.curves()[1].p(), nurbs.curves()[1].knots(), Rv);

					compute_b_ders1_basis(nurbs.curves()[0], span[0], u[0], Ru_ozp, du);
					compute_b_ders1_basis(nurbs.curves()[1], span[1], u[1], Rv_ozp, dv);*/

			/*    std::cout << "----------------------------------" << std::endl;
					std::cout << Ru_ozp.transpose() << std::endl;
					std::cout << Ru.transpose() << std::endl;
					std::cout << "----------------------------------" << std::endl;
					std::cout << Rv_ozp.transpose() << std::endl;
					std::cout << Rv.transpose() << std::endl;*/

			igafem::nurbs2DBasisDers(x, y, nurbs.curves()[0].p(), nurbs.curves()[1].p(),
				nurbs.curves()[0].knots(), nurbs.curves()[1].knots(), weight_vec, R, du, dv);

			compute_nurbs_ders1_basis(nurbs, span, u, weights, r_ozp, der_ozp);

			/*   std::cout << "----------------------------------" << std::endl;
				 std::cout << r_ozp.transpose() << std::endl;
				 std::cout << R.transpose() << std::endl;

				 std::cout << "----------------------------------" << std::endl;
				 std::cout << der_ozp.transpose() << std::endl;
				 std::cout << du.transpose() << std::endl;
				 std::cout << dv.transpose() << std::endl;*/

			for (auto k = 0; k < num_cp; ++k) {
				REQUIRE(r_ozp[k] == Approx(R[k]).margin(1e-6));
				REQUIRE(der_ozp.col(0)[k] == Approx(du[k]).margin(1e-6));
				REQUIRE(der_ozp.col(1)[k] == Approx(dv[k]).margin(1e-6));
			}
		}
	}
}

//BOOST_AUTO_TEST_SUITE(RectangleBasisDerivatives)
//
//BOOST_AUTO_TEST_CASE(RectangleBasisDerivatives_special) {
//  const auto lx = 10.0;
//  const auto ly = 1.0;
//  
//  const auto px = 1u;
//  const auto py = 0u;
//  
//  const auto hx = 4u;
//  const auto hy = 0u;
//
//  squanch::PointContainerTypes<double>::map_type cp;
//  auto rect = create_rectangle(cp, lx, ly, hx, hy, px, py);
//
//  auto num_cp = rect.num_cp();
//
//  std::vector<double> weight_vec(num_cp);
//  Eigen::VectorXd weights(num_cp);
//  Eigen::VectorXd R(num_cp);
//  Eigen::VectorXd du(num_cp);
//  Eigen::VectorXd dv(num_cp);
//
//  Eigen::VectorXd r_ozp(num_cp);
//  Eigen::MatrixXd der_ozp(num_cp, 2);
//
//  std::array<size_t, 2> span;
//  std::array<double, 2> u = { {0.447169, 0} };
//  
//  using namespace squanch;
//
//  span[0] = find_span(rect.curves()[0].knots(), rect.curves()[0].p(), u[0]);
//  span[1] = find_span(rect.curves()[1].knots(), rect.curves()[1].p(), u[1]);
//  gather_weights(rect, cp, span, weights);
//  gather_weights(rect, cp, span, weight_vec);
//
//
//  igafem::nurbs2DBasisDers(u[0], u[1], rect.curves()[0].p(), rect.curves()[1].p(),
//    rect.curves()[0].knots(), rect.curves()[1].knots(), weight_vec, R, du, dv);
//
//  compute_nurbs_ders1_basis(rect, span, u, weights, cp, r_ozp, der_ozp);
//
//     std::cout << "----------------------------------" << std::endl;
//  std::cout << r_ozp.transpose() << std::endl;
//  std::cout << R.transpose() << std::endl;
//
//  std::cout << "----------------------------------" << std::endl;
//  std::cout << der_ozp.transpose() << std::endl;
//  std::cout << du.transpose() << std::endl;
//  std::cout << dv.transpose() << std::endl;
//
//  
//}
//

TEST_CASE("Rectangle derivatives", "[derivatives]")
{
  SECTION("RectangleBasisDerivatives_main") {

    squanch::Model<double> model;

    const auto lx = 25.0;
    const auto ly = 5.0;

    const auto nx = 20u;
    const auto ny = 20u;

    std::vector<unsigned int> p_vec = { 0, 1, 2, 3 };
    std::vector<unsigned int> h_vec = { 15, 25 };

    for (auto ip = 0u; ip < p_vec.size(); ++ip) {
      for (auto ih = 0u; ih < h_vec.size(); ++ih) {
        for (auto jp = 0u; jp < p_vec.size(); ++jp) {
          for (auto jh = 0u; jh < h_vec.size(); ++jh) {
            auto hx = h_vec[ih];
            auto hy = h_vec[jh];
            auto px = p_vec[ip];
            auto py = p_vec[jp];
            squanch::PointContainerTypes<double>::map_type cp;
            auto&& rect = create_rectangle(model, lx, ly, hx, hy, px, py);
            check_nurbs(model, rect, nx, ny);
          }
        }
      }
    }
  }


  SECTION("RectangleBasisDerivatives_prefine_basic") {

    squanch::Model<double> model;

    const auto lx = 25.0;
    const auto ly = 5.0;

    const auto nx = 20u;
    const auto ny = 20u;

    squanch::PointContainerTypes<double>::map_type cp;
    auto&& rect = create_rectangle(model, lx, ly, 2, 2, 1, 1);
    check_nurbs(model, rect, nx, ny);
  }

}