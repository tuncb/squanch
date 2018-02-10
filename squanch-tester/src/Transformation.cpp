#include <iostream>
#include <boost/test/unit_test.hpp>
#include <squanch/Model.h>
#include <squanch/Circle.h>
#include <squanch\Transformation.h>


BOOST_AUTO_TEST_SUITE(transformation)

BOOST_AUTO_TEST_CASE(circle_center)
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
  circle1.origin = Eigen::Vector3d(1,2,3);
  circle1.r = 5.0;
  auto&& ci = circle1.to_nurbs(model);

  Eigen::Vector4d center;
  find_center(model, ci, center);

  BOOST_CHECK_CLOSE( circle1.origin.x(), center.x(), 0.000001);
  BOOST_CHECK_CLOSE( circle1.origin.y(), center.y(), 0.000001);
  BOOST_CHECK_CLOSE( circle1.origin.z(), center.z(), 0.000001);
}

BOOST_AUTO_TEST_CASE(circle_rotate)
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
  auto&& c1 = circle1.to_nurbs(model);

  TransformationTypes<double>::transformation q;
  q = Eigen::AngleAxis<double>(pi/2, Eigen::Vector3d::UnitX());

  transform(model, c1, q);
  
  Circle<double, Eigen::Vector3d> circle2;
  circle2.start_angle = 0.0;
  circle2.end_angle   = 2*pi;
  circle2.N1 = Eigen::Vector3d(1,0,0);
  circle2.N2 = Eigen::Vector3d(0,0,1);
  circle2.origin = Eigen::Vector3d(0,0,0);
  circle2.r = 5.0;
  auto&& c2 = circle2.to_nurbs(model);

  Eigen::Vector4d center1;
  find_center(model, c1, center1);

  Eigen::Vector4d center2;
  find_center(model, c2, center2);

  // +1s are to get rid of division by zeros of BOOST_CHECK_CLOSE

  BOOST_CHECK_CLOSE( center1.x()+1, center2.x()+1, 0.000001);
  BOOST_CHECK_CLOSE( center1.y()+1, center2.y()+1, 0.000001);
  BOOST_CHECK_CLOSE( center1.z()+1, center2.z()+1, 0.000001);

  for (size_t i = 0; i < c1.curves()[0].n(); ++i) {
    Eigen::Vector4d p1 = model.get_point(c1.cpi()(i));
    Eigen::Vector4d p2 = model.get_point(c2.cpi()(i));

    BOOST_CHECK_CLOSE( p1.x()+1, p2.x()+1, 0.000001);
    BOOST_CHECK_CLOSE( p1.y()+1, p2.y()+1, 0.000001);
    BOOST_CHECK_CLOSE( p1.z()+1, p2.z()+1, 0.000001);
    BOOST_CHECK_CLOSE( p1.w()+1, p2.w()+1, 0.000001);
  }
}


BOOST_AUTO_TEST_CASE(circle_translate)
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
  circle1.origin = Eigen::Vector3d(1,2,3);
  circle1.r = 5.0;
  auto&& ci = circle1.to_nurbs(model);

  translate(model, ci, 3.0,4.0,5.0);

  Eigen::Vector4d center;
  find_center(model, ci, center);

  BOOST_CHECK_CLOSE( circle1.origin.x() + 3, center.x(), 0.000001);
  BOOST_CHECK_CLOSE( circle1.origin.y() + 4, center.y(), 0.000001);
  BOOST_CHECK_CLOSE( circle1.origin.z() + 5, center.z(), 0.000001);
}

BOOST_AUTO_TEST_CASE(circle_rotate_origin)
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
  circle1.origin = Eigen::Vector3d(1,1,1);
  circle1.r = 5.0;
  auto&& c1 = circle1.to_nurbs(model);

  TransformationTypes<double>::transformation q;
  q = Eigen::AngleAxis<double>(pi/2, Eigen::Vector3d::UnitX());

  transform_around(model, c1, q, 1.0, 1.0, 1.0);
  
  Circle<double, Eigen::Vector3d> circle2;
  circle2.start_angle = 0.0;
  circle2.end_angle   = 2*pi;
  circle2.N1 = Eigen::Vector3d(1,0,0);
  circle2.N2 = Eigen::Vector3d(0,0,1);
  circle2.origin = Eigen::Vector3d(1,1,1);
  circle2.r = 5.0;
  auto&& c2 = circle2.to_nurbs(model);

  Eigen::Vector4d center1;
  find_center(model, c1, center1);

  Eigen::Vector4d center2;
  find_center(model, c2, center2);

  // +1s are to get rid of division by zeros of BOOST_CHECK_CLOSE

  BOOST_CHECK_CLOSE( center1.x()+1, center2.x()+1, 0.000001);
  BOOST_CHECK_CLOSE( center1.y()+1, center2.y()+1, 0.000001);
  BOOST_CHECK_CLOSE( center1.z()+1, center2.z()+1, 0.000001);

  for (size_t i = 0; i < c1.curves()[0].n(); ++i) {
    Eigen::Vector4d p1 = model.get_point(c1.cpi()(i));
    Eigen::Vector4d p2 = model.get_point(c2.cpi()(i));

    BOOST_CHECK_CLOSE( p1.x()+1, p2.x()+1, 0.000001);
    BOOST_CHECK_CLOSE( p1.y()+1, p2.y()+1, 0.000001);
    BOOST_CHECK_CLOSE( p1.z()+1, p2.z()+1, 0.000001);
    BOOST_CHECK_CLOSE( p1.w()+1, p2.w()+1, 0.000001);
  }
}

BOOST_AUTO_TEST_SUITE_END()
