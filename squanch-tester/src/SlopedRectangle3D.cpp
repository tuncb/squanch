#include "catch.hpp"
#include <squanch/Model.h>
#include <squanch/IgaAlgorithms.h>
#include "igafem.h"


template <typename T> void gather_weights(const squanch::Nurbs<T, 3>& nurbs, const typename squanch::PointContainerTypes<T>::map_type& cp, 
                                         const std::array<unsigned int, 3> span, std::vector<double>& weights)
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
        weights[c++] = cp.find(cpi((int)(span[0] - p1 + j),(int)(span[1] - p2 + i),(int)(span[2] - p3 + k)))->second.w();
      }
    }
  }
}

void create_sloped_rect(squanch::Model<double>& model, squanch::Nurbs<double, 3>& nurbs)
{
  using namespace squanch;
  
  nurbs.curves()[0].set_p(1);
  nurbs.curves()[0].knots().push_back(0.0);
  nurbs.curves()[0].knots().push_back(0.0);
  nurbs.curves()[0].knots().push_back(1.0);
  nurbs.curves()[0].knots().push_back(1.0);
  
  nurbs.curves()[1].set_p(2);
  nurbs.curves()[1].knots().push_back(0.0);
  nurbs.curves()[1].knots().push_back(0.0);
  nurbs.curves()[1].knots().push_back(0.0);
  nurbs.curves()[1].knots().push_back(0.5);
  nurbs.curves()[1].knots().push_back(1.0);
  nurbs.curves()[1].knots().push_back(1.0);
  nurbs.curves()[1].knots().push_back(1.0);

  nurbs.curves()[2].set_p(1);
  nurbs.curves()[2].knots().push_back(0.0);
  nurbs.curves()[2].knots().push_back(0.0);
  nurbs.curves()[2].knots().push_back(1.0);
  nurbs.curves()[2].knots().push_back(1.0);

  int id0, id1, id2, id3, id4, id5, id6, id7, id8, id9, id10, id11, id12, id13, id14, id15;
  
  std::tie(id0, std::ignore) = model.new_point(-2.0 , 0.0, 0.0, 1.0);
  std::tie(id1, std::ignore) = model.new_point(-1.0, 0.0, 0.0, 1.0);
  std::tie(id2, std::ignore) = model.new_point(-2.0, 1.0, 0.0, 1.0);
  std::tie(id3, std::ignore) = model.new_point(-1.0, 0.5, 0.0, 1.0);
  std::tie(id4, std::ignore) = model.new_point(-0.75, 2.0, 0.0, 1.0);
  std::tie(id5, std::ignore) = model.new_point(-0.5, 1.0, 0.0, 1.0);
  std::tie(id6, std::ignore) = model.new_point(0.0, 2.0, 0.0, 1.0);
  std::tie(id7, std::ignore) = model.new_point(0.0, 1.0, 0.0, 1.0);

  std::tie(id8, std::ignore) = model.new_point(-2.0, 0.0, 1.0, 1.0);
  std::tie(id9, std::ignore) = model.new_point(-1.0, 0.0, 1.0, 1.0);
  std::tie(id10, std::ignore) = model.new_point(-2.0, 1.0, 1.0, 1.0);
  std::tie(id11, std::ignore) = model.new_point(-1.0, 0.5, 1.0, 1.0);
  std::tie(id12, std::ignore) = model.new_point(-0.75, 2.0, 1.0, 1.0);
  std::tie(id13, std::ignore) = model.new_point(-0.5, 1.0, 1.0, 1.0);
  std::tie(id14, std::ignore) = model.new_point(0.0, 2.0, 1.0, 1.0);
  std::tie(id15, std::ignore) = model.new_point(0.0, 1.0, 1.0, 1.0);
  
  nurbs.cpi().resize(2,4,2);

  nurbs.cpi()(1,0,0) = id0;
  nurbs.cpi()(0,0,0) = id1;
  nurbs.cpi()(1,1,0) = id2;
  nurbs.cpi()(0,1,0) = id3;
  nurbs.cpi()(1,2,0) = id4;
  nurbs.cpi()(0,2,0) = id5;
  nurbs.cpi()(1,3,0) = id6;
  nurbs.cpi()(0,3,0) = id7;

  nurbs.cpi()(1,0,1) = id8;
  nurbs.cpi()(0,0,1) = id9;
  nurbs.cpi()(1,1,1) = id10;
  nurbs.cpi()(0,1,1) = id11;
  nurbs.cpi()(1,2,1) = id12;
  nurbs.cpi()(0,2,1) = id13;
  nurbs.cpi()(1,3,1) = id14;
  nurbs.cpi()(0,3,1) = id15;
}

TEST_CASE("sloped_rect3D", "[sloped_rect]")
{
  SECTION("sloped_rect3D")
  {
    using namespace squanch;
    squanch::Model<double> model;
    auto&& nurbs = model.new_nurbs<3>();
    create_sloped_rect(model, nurbs);

    const double du = 0.2;
    const double dv = 0.2;
    const double dk = 0.2;

    std::array<double, 3> u;
    u.fill(0.0);
    std::array<unsigned int, 3> span;

    Eigen::VectorXd r(12);
    Eigen::Matrix<double, -1, 3> der(12, 3);
    Eigen::Vector4d point;

    Eigen::VectorXd r2(12);
    Eigen::VectorXd ders2_0(12);
    Eigen::VectorXd ders2_1(12);
    Eigen::VectorXd ders2_2(12);

    std::vector<double> weights;
    Eigen::VectorXd weights2(12);

    //std::cout << "--------------------------------------------------------------" << std::endl;

    while (u[0] < 1.0) {
      u[1] = 0.0;
      while (u[1] < 1.0) {
        u[2] = 0.0;
        while (u[2] < 1.0) {
          span[0] = find_span(nurbs.curves()[0].knots(), nurbs.curves()[0].p(), u[0]);
          span[1] = find_span(nurbs.curves()[1].knots(), nurbs.curves()[1].p(), u[1]);
          span[2] = find_span(nurbs.curves()[2].knots(), nurbs.curves()[2].p(), u[2]);

          gather_weights(model, nurbs, span, weights2);
          for (auto i = 0; i < weights2.rows(); ++i) weights.push_back(weights2[i]);

          compute_nurbs_ders1_basis(nurbs, span, u, weights2, r, der);

          igafem::nurbs3DBasisDers(u[0], u[1], u[2], nurbs.curves()[0].p(), nurbs.curves()[1].p(), nurbs.curves()[2].p(),
            nurbs.curves()[0].knots(), nurbs.curves()[1].knots(), nurbs.curves()[2].knots(),
            weights, r2, ders2_0, ders2_1, ders2_2);

          REQUIRE(r.sum() == Approx(1.0));
          for (size_t i = 0; i < 6; ++i) {
            if (abs(r[i]) > 0.00000001) REQUIRE(r[i] == Approx(r2[i]));
            if (abs(der(i, 0)) > 0.00000001) REQUIRE(der.col(0)[i] == Approx(ders2_0[i]));
            if (abs(der(i, 1)) > 0.00000001) REQUIRE(der.col(1)[i] == Approx(ders2_1[i]));
            if (abs(der(i, 2)) > 0.00000001) REQUIRE(der.col(2)[i] == Approx(ders2_2[i]));
          }

          //compute_model_coordinate(nurbs, cp, span, r, point);
          //std::cout << "u = " << u[0] << " : " << u[1] << " : " << u[2] << " point = " <<  point.transpose() << std::endl;
          u[2] += dk;
        }
        u[1] += dv;
      }
      u[0] += du;
    }
  }
}
