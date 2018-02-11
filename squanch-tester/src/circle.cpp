#include "catch.hpp"
#include <squanch/Model.h>
#include <squanch/IgaAlgorithms.h>
#include <squanch/Circle.h>
#include "DebugTools.h"

#include "igafem.h"

void create_circle(squanch::Model<double>& model, squanch::Nurbs<double, 1>& nurbs)
{
  using namespace squanch;
  
  const double pi = 3.14159265359;
  
  nurbs.curves()[0].knots().push_back(0.0);
  nurbs.curves()[0].knots().push_back(0.0);
  nurbs.curves()[0].knots().push_back(0.0);
  nurbs.curves()[0].knots().push_back(pi/2.0);
  nurbs.curves()[0].knots().push_back(pi/2.0);
  nurbs.curves()[0].knots().push_back(pi);
  nurbs.curves()[0].knots().push_back(pi);
  nurbs.curves()[0].knots().push_back(3.0*pi/2.0);
  nurbs.curves()[0].knots().push_back(3.0*pi/2.0);
  nurbs.curves()[0].knots().push_back(2.0*pi);
  nurbs.curves()[0].knots().push_back(2.0*pi);
  nurbs.curves()[0].knots().push_back(2.0*pi);

  
  int id0, id1, id2, id3, id4, id5, id6, id7, id8;
  
  std::tie(id0, std::ignore) = model.new_point(1.0, 0.0, 0.0, 1.0);
  std::tie(id1, std::ignore) = model.new_point(1.0, 1.0, 0.0, 0.70710678118);
  std::tie(id2, std::ignore) = model.new_point(0.0, 1.0, 0.0, 1.0);
  std::tie(id3, std::ignore) = model.new_point(-1.0, 1.0, 0.0, 0.70710678118);
  std::tie(id4, std::ignore) = model.new_point(-1.0, 0.0, 0.0, 1.0);
  std::tie(id5, std::ignore) = model.new_point(-1.0, -1.0, 0.0, 0.70710678118);
  std::tie(id6, std::ignore) = model.new_point(0.0, -1.0, 0.0, 1.0);
  std::tie(id7, std::ignore) = model.new_point(1.0, -1.0, 0.0, 0.70710678118);
  std::tie(id8, std::ignore) = model.new_point(1.0, 0.0, 0.0, 1.0);

  nurbs.cpi().resize(9);

  nurbs.cpi()(0) = id0;
  nurbs.cpi()(1) = id1;
  nurbs.cpi()(2) = id2;
  nurbs.cpi()(3) = id3;
  nurbs.cpi()(4) = id4;
  nurbs.cpi()(5) = id5;
  nurbs.cpi()(6) = id6;
  nurbs.cpi()(7) = id7;
  nurbs.cpi()(8) = id8;

  nurbs.curves()[0].set_p(2);
}

TEST_CASE("Circle", "[circle]") 
{
  SECTION("Manual Cricle Creation")
  {
    using namespace squanch;

    squanch::Model<double> model;
    auto&& nurbs = model.new_nurbs<1>();

    const double pi = 3.14159265359;
    create_circle(model, nurbs);

    Eigen::VectorXd r(3);
    Eigen::VectorXd r2(3);
    Eigen::VectorXd ders(3);
    Eigen::VectorXd ders2(3);

    Eigen::Vector4d point;
    std::vector<double> weights_vec;
    Eigen::VectorXd weights, weights2(3);

    const double du = 0.1;
    double u = 0.0;
    while (u < 2 * pi - 0.75*du) {
      auto span = find_span(nurbs.curves()[0].knots(), 2, u);
      gather_weights(model, nurbs, span, weights);
      gather_weights(model, nurbs, span, weights2);

      for (auto i = 0; i < weights.rows(); ++i) weights_vec.push_back(weights[i]);

      igafem::nurbs1DBasisDers(u, 2, nurbs.curves()[0].knots(), weights_vec, r2, ders2);
      compute_nurbs_ders1_basis(nurbs, span, u, weights2, r, ders);

      for (size_t i = 0; i <= 2; ++i) {
        REQUIRE(r[i] == Approx(r2[i]));
        REQUIRE(ders(i, 0) == Approx(ders2[i]));
      }

      REQUIRE(r.sum() == Approx(1.0));
      compute_model_coordinate(model, nurbs, span, r, point);
      double distance2 = point.x() * point.x() + point.y() * point.y();
      REQUIRE(distance2 == Approx(1.0));

      u += pi / 12;
    }
  }

  SECTION("Create with the Circle class")
  {
    using namespace squanch;

    squanch::Model<double> model;

    const double pi = 3.14159265359;

    Circle<double, Eigen::Vector3d> circle;
    circle.start_angle = 0.0;
    circle.end_angle = 2 * pi;
    circle.N1 = Eigen::Vector3d(1, 0, 0);
    circle.N2 = Eigen::Vector3d(0, 1, 0);
    circle.origin = Eigen::Vector3d(0, 0, 0);
    circle.r = 1.0;
    auto&& nurbs = circle.to_nurbs(model);

    //print_cps(nurbs, model.control_points());

    Eigen::VectorXd weights(nurbs.curves()[0].p() + 1);
    Eigen::VectorXd r(nurbs.curves()[0].p() + 1);
    Eigen::Vector4d point;
    const double du = 0.1;
    double u = 0.0;
    while (u <= 1.0) {
      auto span = find_span(nurbs.curves()[0].knots(), nurbs.curves()[0].p(), u);
      gather_weights(model, nurbs, span, weights);

      compute_b_basis(nurbs, span, u, r);

      REQUIRE(r.sum() == Approx(1.0));
      compute_model_coordinate(model, nurbs, span, r, point);
      from_homogeneous(point);
      double distance2 = point.x() * point.x() + point.y() * point.y();
      REQUIRE(distance2 == Approx(1.0));

      u += du;
    }

    REQUIRE(nurbs.curves()[0].is_closed() == true);
    REQUIRE(nurbs.cpi()(0) == nurbs.cpi()(8));
  }
}
