#include <boost/test/unit_test.hpp>
#include <squanch/Model.h>
#include <squanch/IgaAlgorithms.h>

#include "igafem.h"

template <typename T> void gather_weights(const squanch::Nurbs<T, 2>& nurbs, const typename squanch::PointContainerTypes<T>::map_type& cp, 
                                         const std::array<unsigned int, 2> span, std::vector<double>& weights)
{
  auto p1 = nurbs.curves()[0].p();
  auto p2 = nurbs.curves()[1].p();
  auto&& cpi = nurbs.cpi();
  
  weights.resize((p1+1)*(p2+1));

  unsigned int c = 0;
  for (unsigned int i = 0; i <= p2; ++i) {
    for (unsigned int j = 0; j <= p1; ++j) {
      weights[c++] = cp.find(cpi((int)(span[0] - p1 + j),(int)(span[1] - p2 + i)))->second.w();
    }
  }
}

void create_sloped_rect(squanch::Model<double>& model, squanch::Nurbs<double, 2>& nurbs)
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

  int id0, id1, id2, id3, id4, id5, id6, id7;
  
  std::tie(id0, std::ignore) = model.new_point(-2.0 , 0.0, 0.0, 1.0);
  std::tie(id1, std::ignore) = model.new_point(-1.0, 0.0, 0.0, 1.0);
  std::tie(id2, std::ignore) = model.new_point(-2.0, 1.0, 0.0, 1.0);
  std::tie(id3, std::ignore) = model.new_point(-1.0, 0.5, 0.0, 1.0);
  std::tie(id4, std::ignore) = model.new_point(-0.75, 2.0, 0.0, 1.0);
  std::tie(id5, std::ignore) = model.new_point(-0.5, 1.0, 0.0, 1.0);
  std::tie(id6, std::ignore) = model.new_point(0.0, 2.0, 0.0, 1.0);
  std::tie(id7, std::ignore) = model.new_point(0.0, 1.0, 0.0, 1.0);

  nurbs.cpi().resize(2,4);

  nurbs.cpi()(1,0) = id0;
  nurbs.cpi()(0,0) = id1;
  nurbs.cpi()(1,1) = id2;
  nurbs.cpi()(0,1) = id3;
  nurbs.cpi()(1,2) = id4;
  nurbs.cpi()(0,2) = id5;
  nurbs.cpi()(1,3) = id6;
  nurbs.cpi()(0,3) = id7;
}

BOOST_AUTO_TEST_SUITE(sloped_rect)
BOOST_AUTO_TEST_CASE(sloped_rect)
{
  using namespace squanch;

  squanch::Model<double> model;

  auto&& nurbs = model.new_nurbs<2>();
  create_sloped_rect(model, nurbs);

  const double du = 0.1; 
  const double dv = 0.1; 

  std::array<double, 2> u;
  u.fill(0.0);
  std::array<unsigned int, 2> span;

  Eigen::VectorXd r(6);
  Eigen::VectorXd r2(6);
  Eigen::VectorXd ders2_0(6);
  Eigen::VectorXd ders2_1(6);
  Eigen::Matrix<double, -1, 2> der(6, 2);
  Eigen::Vector4d point;

  std::vector<double> weights_vec;
  Eigen::VectorXd weights, weights2(6);

  while (u[0] < 1.0) {
    u[1] = 0.0;
    while (u[1] < 1.0) {
      span[0] = find_span(nurbs.curves()[0].knots(), nurbs.curves()[0].p(), u[0]);
      span[1] = find_span(nurbs.curves()[1].knots(), nurbs.curves()[1].p(), u[1]);

      gather_weights(model, nurbs, span, weights);
      gather_weights(model, nurbs, span, weights2);
			
			for (auto i = 0; i < weights2.rows(); ++i) weights_vec.push_back(weights[i]);

      compute_nurbs_ders1_basis(nurbs, span, u, weights2, r, der);
      igafem::nurbs2DBasisDers(u[0], u[1], nurbs.curves()[0].p(), nurbs.curves()[1].p(), nurbs.curves()[0].knots(), nurbs.curves()[1].knots(),
        weights_vec, r2, ders2_0, ders2_1);

      //std::cout << der.col(0).transpose() << std::endl;
      //std::cout << ders2_0.transpose() << std::endl;

      BOOST_CHECK_CLOSE( r.sum(), 1.0, 0.00000001);
      for (size_t i = 0; i < 6; ++i) {
        if (abs(r[i]) > 0.00000001) BOOST_CHECK_CLOSE(r[i], r2[i], 0.00000001);
        if (abs(der.col(0)[i]) > 0.00000001) BOOST_CHECK_CLOSE(der.col(0)[i], ders2_0[i], 0.00001);
        if (abs(der.col(1)[i]) > 0.00000001) BOOST_CHECK_CLOSE(der.col(1)[i], ders2_1[i], 0.00001);
      }


      //std::cout << r.transpose() << std::endl;
      //compute_model_coordinate(nurbs, cp, span, r, point);
      //std::cout << "u = " << u[0] << " : " << u[1] << " point = " <<  point.transpose() << std::endl;

      u[1] += dv;
    }
    u[0] += du;
  }


}

BOOST_AUTO_TEST_SUITE_END()