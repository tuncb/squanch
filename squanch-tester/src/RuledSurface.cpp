#include <boost/test/unit_test.hpp>
#include <Eigen/core>
#include <squanch/Model.h>
#include <squanch/IgaAlgorithms.h>
#include <squanch/Circle.h>
#include <squanch/RectangleXy.h>
#include <squanch/RuledSurface.h>

BOOST_AUTO_TEST_SUITE(ruled_surface)

BOOST_AUTO_TEST_CASE(circle_inside_out)
{
  using namespace squanch;
  using namespace Eigen;

  squanch::Model<double> model;

  const double pi = 3.14159265359;

  Circle<double, Eigen::Vector3d> circle1;
  circle1.start_angle = 0.0;
  circle1.end_angle   = 2*pi;
  circle1.N1 = Eigen::Vector3d(1,0,0);
  circle1.N2 = Eigen::Vector3d(0,1,0);
  circle1.origin = Eigen::Vector3d(0,0,0);
  circle1.r = 5.0;
  auto&& ci = circle1.to_nurbs(model);

  Circle<double, Eigen::Vector3d> circle2;
  circle2.start_angle = 0.0;
  circle2.end_angle   = 2*pi;
  circle2.N1 = Eigen::Vector3d(1,0,0);
  circle2.N2 = Eigen::Vector3d(0,1,0);
  circle2.origin = Eigen::Vector3d(0,0,0);
  circle2.r = 10.0;
  auto&& co = circle2.to_nurbs(model);

  auto&& area = create_ruled_surface(model, ci, co);

}

BOOST_AUTO_TEST_CASE(circle_inside_rectangle)
{
  using namespace squanch;
  using namespace Eigen;

  squanch::Model<double> model;

  const double pi = 3.14159265359;

  Circle<double, Eigen::Vector3d> circle1;
  circle1.start_angle = 0.0;
  circle1.end_angle   = 2*pi;
  circle1.N1 = Eigen::Vector3d(1,0,0);
  circle1.N2 = Eigen::Vector3d(0,1,0);
  circle1.origin = Eigen::Vector3d(10,10,0);
  circle1.r = 5.0;
  auto&& ci = circle1.to_nurbs(model);

  RectangleXy<double> rect(20,20);

  auto&& co = rect.to_nurbs(model);

  auto&& area = create_ruled_surface(model,ci, co);

}

BOOST_AUTO_TEST_SUITE_END()