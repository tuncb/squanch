#include "catch.hpp"
#include <Eigen/core>
#include <squanch/IgaAlgorithms.h>
#include <squanch/Line.h>
#include <squanch/Extrude.h>
#include "DebugTools.h"

TEST_CASE("exture", "[extrude]")
{
  SECTION("rectangle")
  {
    using namespace squanch;
    using namespace Eigen;
    // create a line
    Line<double, Vector3d> line;
    line.p1 = Vector3d(5, 0, 0);
    line.p2 = Vector3d(10, 0, 0);

    squanch::Model<double> model;
    Eigen::VectorXd r = Eigen::VectorXd::Zero(4);

    auto&& nline = line.to_nurbs(model);
    Eigen::Vector4d vec(0, 10, 0, 0);
    auto&& rect = squanch::extrude(model, nline, vec);

    std::array<double, 2> u;
    std::array<unsigned int, 2> span;
    Eigen::VectorXd weights(4);
    Eigen::Vector4d p;
    for (size_t i = 0; i < 11; ++i) {
      u[0] = i * 0.1;
      span[0] = find_span(rect.curves()[0].knots(), rect.curves()[0].p(), u[0]);
      for (size_t j = 0; j < 11; ++j) {
        u[1] = j * 0.1;
        span[1] = find_span(rect.curves()[1].knots(), rect.curves()[1].p(), u[1]);
        gather_weights(model, rect, span, weights);
        compute_nurbs_basis(rect, span, u, weights, r);
        compute_model_coordinate(model, rect, span, r, p);

        REQUIRE(p.x() == Approx(5.0 + u[0] * 5));
        REQUIRE(p.y() == Approx(u[1] * 10));
        REQUIRE(p.z() == Approx(0));
        REQUIRE(p.w() == Approx(1));
      }
    }
  }

  SECTION("cube")
  {
    using namespace squanch;
    using namespace Eigen;
    // create a line
    Line<double, Vector3d> line;
    line.p1 = Vector3d(5, 0, 0);
    line.p2 = Vector3d(10, 0, 0);

    int id_s = 0;
    squanch::Model<double> model;
    Eigen::VectorXd r = Eigen::VectorXd::Zero(8);
    Eigen::VectorXd weights(8);

    auto&& nline = line.to_nurbs(model);
    Eigen::Vector4d vec(0, 10, 0, 0);
    auto&& rect = squanch::extrude(model, nline, vec);
    Eigen::Vector4d vec2(0, 0, 7.5, 0);
    auto&& cube = squanch::extrude(model, rect, vec2);

    std::array<double, 3> u;
    std::array<unsigned int, 3> span;
    Eigen::Vector4d p;
    for (size_t i = 0; i < 11; ++i) {
      u[0] = i * 0.1;
      span[0] = find_span(cube.curves()[0].knots(), cube.curves()[0].p(), u[0]);
      for (size_t j = 0; j < 11; ++j) {
        u[1] = j * 0.1;
        span[1] = find_span(cube.curves()[1].knots(), cube.curves()[1].p(), u[1]);
        for (size_t k = 0; k < 11; ++k) {
          u[2] = k * 0.1;
          span[2] = find_span(cube.curves()[2].knots(), cube.curves()[2].p(), u[2]);

          gather_weights(model, cube, span, weights);
          compute_nurbs_basis(cube, span, u, weights, r);
          compute_model_coordinate(model, cube, span, r, p);

          REQUIRE(p.x() == Approx(5.0 + u[0] * 5));
          REQUIRE(p.y() == Approx(u[1] * 10));
          REQUIRE(p.z() == Approx(u[2] * 7.5));
          REQUIRE(p.w() == Approx(1));
        }
      }
    }
  }
}