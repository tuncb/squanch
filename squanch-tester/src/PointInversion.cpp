#include "catch.hpp"
#include <squanch/Model.h>
#include <squanch/Line.h>
#include <squanch/PointInversion.h>
#include <squanch/Circle.h>
#include <squanch/Extrude.h>

TEST_CASE("Point inversion", "[pointInversion]")
{
  SECTION("line") {

    //boost::debug::break_memory_alloc(1440);
    using namespace squanch;
    squanch::Model<double> model;

    Eigen::Vector3d p1(0, 0, 0);
    Eigen::Vector3d p2(10, 0, 0);

    Line<double, Eigen::Vector3d> line(p1, p2);
    auto&& nurbs = line.to_nurbs(model);

    const auto num_spaces = 20;
    const auto num_points = num_spaces + 1;
    const auto dx = 10.0 / num_spaces;
    const auto du = 1.0 / num_spaces;
    Eigen::Vector4d p(0, 0, 0, 1);
    double u = 0.0;
    for (auto i = 0; i < num_points; ++i) {
      p[0] = i * dx;
      u = i * du;
      auto u_calculated = inverse_point_nr(model, nurbs, 0.0, p);
      REQUIRE(u == Approx(u_calculated));
    }
  }

  SECTION("circle") {
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

    std::vector<double> u_real = { 0, 0.25, 0.5, 0.75 };
    PointContainerTypes<double>::vector_type p_real;
    p_real.push_back({ 1, 0, 0, 1 });
    p_real.push_back({ 0, 1, 0, 1 });
    p_real.push_back({ -1, 0, 0, 1 });
    p_real.push_back({ 0, -1, 0, 1 });

    auto num_points = u_real.size();

    for (auto i = 0; i < num_points; ++i) {
      auto u_calculated = inverse_point_nr(model, nurbs, p_real[i]);
      REQUIRE(u_real[i] == Approx(u_calculated));
    }
  }

  SECTION("rectangular_surface_1d")
  {
    using namespace squanch;
    squanch::Model<double> model;

    Eigen::Vector3d p1(0, 0, 0);
    Eigen::Vector3d p2(10, 0, 0);

    Line<double, Eigen::Vector3d> line(p1, p2);
    auto&& nurbs_line = line.to_nurbs(model);

    PointContainerTypes<double>::point_type dir(0, 20, 0, 1);

    auto&& nurbs = extrude(model, nurbs_line, dir);
    NurbsSkeleton<double, 1> skeleton(nurbs, 1, LineEnds::End);

    const auto num_spaces = 20;
    const auto num_points = num_spaces + 1;
    const auto dx = 20.0 / num_spaces;
    const auto du = 1.0 / num_spaces;
    Eigen::Vector4d p(10, 0, 0, 1);
    double u = 0.0;
    for (auto i = 0; i < num_points; ++i) {
      p[1] = i * dx;
      u = i * du;
      auto u_calculated = inverse_point_nr(model, skeleton, 0.0, p);
      REQUIRE(u == Approx(u_calculated));
    }
  }

  SECTION("rectangular_surface_2d")
  {
    using namespace squanch;
    squanch::Model<double> model;

    Eigen::Vector3d p1(0, 0, 0);
    Eigen::Vector3d p2(10, 0, 0);

    Line<double, Eigen::Vector3d> line(p1, p2);
    auto&& nurbs_line = line.to_nurbs(model);
    PointContainerTypes<double>::point_type dir(0, 20, 0, 1);

    auto&& nurbs = extrude(model, nurbs_line, dir);

    const auto num_spaces = 20;
    const auto num_points = num_spaces + 1;

    const auto dx = 10.0 / num_spaces;
    const auto dy = 20.0 / num_spaces;

    const auto du = 1.0 / num_spaces;
    Eigen::Vector4d p(0, 0, 0, 1);
    std::array<double, 2> u = { { 0.0, 0.0 } };
    for (auto i = 0; i < num_points; ++i) {
      p[0] = i * dx;
      u[0] = i * du;
      for (auto j = 0u; j < num_points; ++j) {
        p[1] = j * dy;
        u[1] = j * du;
        auto u_calculated = inverse_point_nr(model, nurbs, p);
        REQUIRE(u[0] == Approx(u_calculated[0]));
        REQUIRE(u[1] == Approx(u_calculated[1]));
      }
      u[1] = 0.0;
      p[1] = 0.0;
    }
  }

  SECTION("cylinder") {
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
    auto&& nurbs_c = circle.to_nurbs(model);

    PointContainerTypes<double>::point_type dir(0, 0, 20, 1);
    auto&& nurbs = extrude(model, nurbs_c, dir);

    const auto num_spaces = 20;
    const auto num_points = num_spaces + 1;
    const auto dz = 20.0 / num_spaces;
    const auto du = 1.0 / num_spaces;

    std::vector<double> u_real = { 0, 0.25, 0.5, 0.75 };
    PointContainerTypes<double>::vector_type p_real;
    p_real.push_back({ 1, 0, 0, 1 });
    p_real.push_back({ 0, 1, 0, 1 });
    p_real.push_back({ -1, 0, 0, 1 });
    p_real.push_back({ 0, -1, 0, 1 });

    Eigen::Vector4d p(0, 0, 0, 1);
    std::array<double, 2> u = { { 0.0, 0.0 } };
    for (auto i = 0; i < num_points; ++i) {
      u[1] = du * i;
      for (auto j = 0u; j < u_real.size(); ++j) {
        u[0] = u_real[j];
        auto p = p_real[j];
        p[2] = dz * i;
        auto u_calculated = inverse_point_nr(model, nurbs, p);
        // as the circle is closed we can go to zero by getting close to 1
        if (u[0] == 0) {
          if (u_calculated[0] > 0.5) u_calculated[0] -= 1;
          // avoid division by zero
          REQUIRE(u[0] + 0.1 == Approx(u_calculated[0] + 0.1));
        }
        else {
          REQUIRE(u[0] == Approx(u_calculated[0]));
        }
        REQUIRE(u[1] == Approx(u_calculated[1]));
      }
    }
  }
}