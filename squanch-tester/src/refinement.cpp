#include "catch.hpp"
#include <squanch/Nurbs.h>
#include <iostream>
#include <squanch/Model.h>
#include <squanch/Refinement.h>
#include <squanch/Circle.h>
#include <squanch/Line.h>
#include <squanch/Extrude.h>
#include <squanch\RuledSurface.h>
#include <algorithm>
#include "DebugTools.h"


TEST_CASE("line_refinement", "[refinement]")
{
  SECTION("lr_insert_knot")
  {
    using namespace squanch;

    Eigen::Vector3f p1(0, 0, 0);
    Eigen::Vector3f p2(10, 0, 0);

    squanch::Model<float> model;
    Line<float, Eigen::Vector3f> line(p1, p2);

    auto&& nurbs = line.to_nurbs(model);

    std::vector<float> u;
    float du = 0.1f;
    for (size_t i = 1; i < 10; ++i) { u.push_back(i*du); }

    squanch::h_refine(model, nurbs, u);

    for (size_t i = 0; i < nurbs.cpi().size(); ++i) {
      auto&& p = model.get_point(nurbs.cpi()(i));
      REQUIRE(p.x() == Approx(i * 1.0f));
      REQUIRE(p.y() == Approx(0.0f));
      REQUIRE(p.z() == Approx(0.0f));
      REQUIRE(p.w() == Approx(1.0f));
    }

  }

  SECTION("lr_insert_knot_multiple_spans")
  {
    using namespace squanch;
    squanch::Model<double> model;

    Eigen::Vector3d p1(0, 0, 0);
    Eigen::Vector3d p2(10, 0, 0);

    Line<double, Eigen::Vector3d> line(p1, p2);

    auto&& nurbs = line.to_nurbs(model);

    std::vector<double> u;
    u.push_back(0.25); u.push_back(0.50); u.push_back(0.75);
    squanch::h_refine(model, nurbs, u);

    u.clear();
    u.push_back(0.125); u.push_back(0.375); u.push_back(0.625); u.push_back(0.875);
    squanch::h_refine(model, nurbs, u);

    for (size_t i = 0; i < nurbs.cpi().size(); ++i) {
      auto&& p = model.get_point(nurbs.cpi()(i));
      REQUIRE(p.x() == Approx(i * 1.25));
      REQUIRE(p.y() == Approx(0.0));
      REQUIRE(p.z() == Approx(0.0));
      REQUIRE(p.w() == Approx(1.0));
    }
  }

  SECTION("line_h_refine_elements")
  {
    using namespace squanch;
    squanch::Model<double> model;

    Eigen::Vector3d p1(0, 0, 0);
    Eigen::Vector3d p2(10, 0, 0);

    Line<double, Eigen::Vector3d> line(p1, p2);

    auto&& nurbs = line.to_nurbs(model);

    squanch::h_refine_elements(model, nurbs, 3);

    REQUIRE(nurbs.curves()[0].knots().size() == 6);
    REQUIRE(nurbs.curves()[0].n() == 4);

  }

  SECTION("p_refine_line")
  {
    using namespace squanch;
    squanch::Model<double> model;

    Eigen::Vector3d p1(0, 0, 0);
    Eigen::Vector3d p2(10, 0, 0);

    Line<double, Eigen::Vector3d> line(p1, p2);

    auto&& nurbs = line.to_nurbs(model);

    p_refine(model, nurbs, 3);

    auto check_line = [](const squanch::Model<double>& model, Nurbs<double, 1>& nrb)
    {
      const size_t n = 20;
      Eigen::Vector4d p;
      auto vec_size = (nrb.curves()[0].p() + 1);
      Eigen::VectorXd r = Eigen::VectorXd::Zero(vec_size);
      std::array<size_t, 1> span;
      double du = 1.0 / n;
      for (size_t i = 0; i < n; ++i) {
        double u = i * du;
        span[0] = find_span(nrb.curves()[0].knots(), nrb.curves()[0].p(), u);
        compute_b_basis(nrb, span[0], u, r);
        compute_model_coordinate(model, nrb, span[0], r, p);
        from_homogeneous(p);

        REQUIRE(p.x() == Approx(10.0 * i * du));
        REQUIRE(p.y() == Approx(0.0));
        REQUIRE(p.z() == Approx(0.0));
        REQUIRE(p.w() == Approx(1.0));
      }
    };

    check_line(model, nurbs);
  }

  SECTION("lr_h_refine_circle")
  {
    using namespace squanch;
    squanch::Model<double> model;

    auto&& nurbs = model.new_nurbs<1>();
    const double pi = 3.14159265359;

    Circle<double, Eigen::Vector3d> circle;
    circle.start_angle = 0.0;
    circle.end_angle = 2 * pi;
    circle.N1 = Eigen::Vector3d(1, 0, 0);
    circle.N2 = Eigen::Vector3d(0, 1, 0);
    circle.origin = Eigen::Vector3d(0, 0, 0);
    circle.r = 1.0;
    nurbs = circle.to_nurbs(model);

    const double piq = 0.125;

    std::vector<double> k;
    k.push_back(piq); k.push_back(3.0*piq); k.push_back(5.0*piq); k.push_back(7.0*piq);
    squanch::h_refine(model, nurbs, k);

    Eigen::VectorXd weights;
    Eigen::VectorXd r(3);
    Eigen::Vector4d point;
    const double du = 0.05;
    double u = 0.0;
    while (u < 1.0 - du) {
      auto span = find_span(nurbs.curves()[0].knots(), 2, u);
      gather_weights(model, nurbs, span, weights);

      compute_b_basis(nurbs, span, u, r);

      REQUIRE(r.sum() == Approx(1.0));
      compute_model_coordinate(model, nurbs, span, r, point);
      from_homogeneous(point);
      double distance2 = point.x() * point.x() + point.y() * point.y();
      REQUIRE(distance2 == Approx(1.0));

      u += pi / 12;
    }

    //for (size_t i = 0; i < nurbs.cpi().num_elements(); ++i) {
    //  auto&& p = cp[nurbs.cpi()[i]];
    //  std::cout << p.transpose() << std::endl;
    //}
  }
  SECTION("p_refine_circle")
  {
    using namespace squanch;

    squanch::Model<double> model;
    auto&& nurbs = model.new_nurbs<1>();

    const double pi = 3.14159265359;

    Circle<double, Eigen::Vector3d> circle;
    circle.start_angle = 0.0;
    circle.end_angle = 2 * pi;
    circle.N1 = Eigen::Vector3d(1, 0, 0);
    circle.N2 = Eigen::Vector3d(0, 1, 0);
    circle.origin = Eigen::Vector3d(0, 0, 0);
    circle.r = 1.0;
    nurbs = circle.to_nurbs(model);

    //for (auto&& pair : cp) {
    //  std::cout << pair.second << std::endl;
    //}

    squanch::p_refine(model, nurbs, 1);
    squanch::p_refine(model, nurbs, 1);

    Eigen::VectorXd r(nurbs.curves()[0].p() + 1);
    Eigen::Vector4d point;
    const double du = 0.1;
    double u = 0.0;
    while (u <= 1.0) {
      auto span = find_span(nurbs.curves()[0].knots(), nurbs.curves()[0].p(), u);

      compute_b_basis(nurbs, span, u, r);

      REQUIRE(r.sum() == Approx(1.0));
      compute_model_coordinate(model, nurbs, span, r, point);
      from_homogeneous(point);
      double distance2 = point.x() * point.x() + point.y() * point.y();
      REQUIRE(distance2 == Approx(1.0));

      u += du;
    }

    //std::cout << "------------------------------------------" << std::endl;

    //for (size_t i = 0; i < nurbs.curves()[0].knots().size(); ++i) {
    //  std::cout << nurbs.curves()[0].knots()[i] << std::endl;
    //}
    //
    //std::cout << "------------------------------------------" << std::endl;

    //for (size_t i = 0; i < nurbs.cpi().num_elements(); ++i) {
    //  auto&& p = cp[nurbs.cpi()[i]];
    //  std::cout << p.transpose() << std::endl;
    //}
  }

  SECTION("h_p_refine_rectangle")
  {
    using namespace squanch;
    using namespace Eigen;
    // create rectangle
    Line<double, Vector3d> line;
    line.p1 = Vector3d(5, 0, 0);
    line.p2 = Vector3d(10, 0, 0);

    squanch::Model<double> model;

    auto&& nline = line.to_nurbs(model);
    Eigen::Vector4d vec(0, 10, 0, 0);
    auto&& rect = squanch::extrude(model, nline, vec);

    auto check_rect = [](const squanch::Model<double>& model, Nurbs<double, 2>& nrb)
    {
      std::array<double, 2> u;
      std::array<unsigned int, 2> span;
      Eigen::Vector4d p;

      auto vec_size = (nrb.curves()[0].p() + 1) * (nrb.curves()[1].p() + 1);
      Eigen::VectorXd r = Eigen::VectorXd::Zero(vec_size);
      for (size_t i = 0; i < 11; ++i) {
        u[0] = i * 0.1;
        span[0] = find_span(nrb.curves()[0].knots(), nrb.curves()[0].p(), u[0]);
        for (size_t j = 0; j < 11; ++j) {
          u[1] = j * 0.1;
          span[1] = find_span(nrb.curves()[1].knots(), nrb.curves()[1].p(), u[1]);
          compute_b_basis(nrb, span, u, r);
          compute_model_coordinate(model, nrb, span, r, p);
          from_homogeneous(p);

          REQUIRE(p.x() == Approx(5.0 + u[0] * 5));
          REQUIRE(p.y() == Approx(u[1] * 10));
          REQUIRE(p.z() == Approx(0));
          REQUIRE(p.w() == Approx(1));
        }
      }
    };


    auto print_cps = [](const squanch::Model<double>& model, Nurbs<double, 2>& nrb)
    {
      blitz::TinyVector < int, 2 > index;
      auto arr = nrb.cpi().shape();
      for (index[0] = 0; index[0] < arr[0]; ++index[0]) {
        for (index[1] = 0; index[1] < arr[1]; ++index[1]) {
          auto id = nrb.cpi()(index);
          auto& p = model.get_point(id);
          std::cout << p.transpose() << std::endl;
        }
      }
    };

    check_rect(model, rect);

    std::vector<double> x;
    float du = 0.1f;
    for (size_t i = 1; i < 10; ++i) { x.push_back(i*du); }

    std::array<unsigned int, 2> pr;
    pr[0] = 1; pr[1] = 2;

    squanch::p_refine(model, rect, pr);
    //print_cps(rect, cp);
    check_rect(model, rect);
    squanch::h_refine(model, rect, x, x);
    check_rect(model, rect);
  }

  SECTION("h_p_refine_cube")
  {
    using namespace squanch;
    using namespace Eigen;
    // create a line
    Line<double, Vector3d> line;
    line.p1 = Vector3d(5, 0, 0);
    line.p2 = Vector3d(10, 0, 0);

    squanch::Model<double> model;

    auto&& nline = line.to_nurbs(model);
    Eigen::Vector4d vec(0, 10, 0, 0);
    auto&& rect = squanch::extrude(model, nline, vec);
    Eigen::Vector4d vec2(0, 0, 7.5, 0);
    auto&& cube = squanch::extrude(model, rect, vec2);


    auto check_cube = [](const squanch::Model<double>& model, Nurbs<double, 3>& cube)
    {
      std::array<double, 3> u;
      std::array<unsigned int, 3> span;
      Eigen::Vector4d p;
      auto vec_size = (cube.curves()[0].p() + 1) * (cube.curves()[1].p() + 1) * (cube.curves()[2].p() + 1);
      Eigen::VectorXd r = Eigen::VectorXd::Zero(vec_size);
      for (size_t i = 0; i < 11; ++i) {
        u[0] = i * 0.1;
        span[0] = find_span(cube.curves()[0].knots(), cube.curves()[0].p(), u[0]);
        for (size_t j = 0; j < 11; ++j) {
          u[1] = j * 0.1;
          span[1] = find_span(cube.curves()[1].knots(), cube.curves()[1].p(), u[1]);
          for (size_t k = 0; k < 11; ++k) {
            u[2] = k * 0.1;
            span[2] = find_span(cube.curves()[2].knots(), cube.curves()[2].p(), u[2]);

            compute_b_basis(cube, span, u, r);
            compute_model_coordinate(model, cube, span, r, p);
            from_homogeneous(p);

            REQUIRE(p.x() == Approx(5.0 + u[0] * 5));
            REQUIRE(p.y() == Approx(u[1] * 10));
            REQUIRE(p.z() == Approx(u[2] * 7.5));
            REQUIRE(p.w() == Approx(1));
          }
        }
      }
    };

    check_cube(model, cube);

    std::vector<double> x;
    float du = 0.1f;
    for (size_t i = 1; i < 10; ++i) { x.push_back(i*du); }

    std::array<unsigned int, 3> pr;
    pr[0] = 1; pr[1] = 2; pr[2] = 3;

    squanch::p_refine(model, cube, 1, 3);
    check_cube(model, cube);
    squanch::h_refine(model, cube, x, x, x);
    check_cube(model, cube);
  }

  SECTION("ThreeDCurvedBeamRefinement")
  {
    using namespace squanch;
    typedef float ValType;
    squanch::Model<ValType> model;

    auto create_qarc = [&model](ValType r) -> squanch::Nurbs<ValType, 1>& {
      auto pi = (ValType)(3.14159265);
      Circle<ValType, Eigen::Matrix<ValType, 3, 1>> c;
      c.origin << 0, 0, 0;
      c.start_angle = 0;
      c.end_angle = pi / 2;
      c.N1 << 1, 0, 0;
      c.N2 << 0, 1, 0;
      c.r = r;

      return c.to_nurbs(model);
    };

    auto&& arci = create_qarc((ValType)4.12);
    auto&& arco = create_qarc((ValType)4.32);
    auto&& surface = create_ruled_surface(model, arci, arco);
    Eigen::Matrix<ValType, 4, 1> vec(0, 0, 0.1, 0);
    auto&& prism = extrude(model, surface, vec);
    h_refine_elements(model, prism, 5, 2, 2);

    //print_cps(model, prism);
  }
}