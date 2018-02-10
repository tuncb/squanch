#include <boost/test/unit_test.hpp>
#include <squanch/Model.h>
#include <squanch/Line.h>

BOOST_AUTO_TEST_SUITE(model)

BOOST_AUTO_TEST_CASE(model_clear) {
  using namespace squanch;
  using namespace Eigen;
  
  squanch::Model<double> model;

  Line<double, Vector3d> line;
  line.p1 = Vector3d(0, 0, 0);
  line.p2 = Vector3d(10, 0, 0);

  auto&& nline1 = line.to_nurbs(model);
  auto&& nline2 = line.to_nurbs(model);
  auto&& nline3 = line.to_nurbs(model);

  BOOST_CHECK_EQUAL(model.patches().size(), 3);
  BOOST_CHECK_EQUAL(model.num_nodes(), 6);

  model.clear();

  BOOST_CHECK_EQUAL(model.patches().size(), 0);
  BOOST_CHECK_EQUAL(model.num_nodes(), 0);
}

BOOST_AUTO_TEST_SUITE_END()
