#include "catch.hpp"
#include <squanch/Model.h>
#include <squanch/Line.h>
#include <squanch/Circle.h>
#include <squanch/IgaAlgorithms.h>
#include <squanch/CurveAlgorithms.h>

TEST_CASE("Curve Creation", "[curve]")
{
  SECTION("Line Creation")
  {
    using namespace squanch;
    squanch::Model<double> model;

    Eigen::Vector3d p1(0, 0, 0);
    Eigen::Vector3d p2(10, 0, 0);

    Line<double, Eigen::Vector3d> line(p1, p2);
    auto&& nurbs = line.to_nurbs(model);
    auto&& knots = nurbs.curves()[0].knots();
    auto p = nurbs.curves()[0].p();

    const auto num_spaces = 20;
    const auto num_points = num_spaces + 1;
    const auto dx = 10.0 / num_spaces;
    const auto du = 1.0 / num_spaces;
    Eigen::Vector4d point(0, 0, 0, 1);
    Eigen::Vector4d calculated_point(0, 0, 0, 1);
    double u = 0.0;
    Eigen::VectorXd n(p + 1);
    for (auto i = 0; i < num_points; ++i) {
      point[0] = i * dx;
      u = i * du;
      auto span = find_span(knots, p, u);
      compute_b_basis(nurbs.curves()[0], span, u, n);
      calculated_point = compute_bspline_point(model, nurbs, u, span, n);

      for (auto j = 0; j < 3; ++j) {
        REQUIRE(point[j] == Approx(calculated_point[j]));
      }
    }
  }

  SECTION("Circle Creation")
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

    auto&& knots = nurbs.curves()[0].knots();
    auto p = nurbs.curves()[0].p();

    const auto num_spaces = 20;
    const auto num_points = num_spaces + 1;
    const auto du = 1.0 / num_spaces;
    const auto dangle = 2 * pi / num_spaces;

    Eigen::Vector4d point(0, 0, 0, 1);
    Eigen::Vector4d calculated_point(0, 0, 0, 1);
    double u = 0.0;
    Eigen::VectorXd n(p + 1);
    for (auto i = 0; i < num_points; ++i) {
      u = i * du;

      auto span = find_span(knots, p, u);
      compute_b_basis(nurbs.curves()[0], span, u, n);

      calculated_point = compute_nurbs_point(model, nurbs, u, span, n);

      double distance2 = calculated_point.x() * calculated_point.x() + calculated_point.y() * calculated_point.y();
      REQUIRE(distance2 == Approx(1.0));
    }
  }

}