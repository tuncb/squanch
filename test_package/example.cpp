#include <squanch/Model.h>
#include <squanch/PointContainerTypes.h>
#include <squanch/Circle.h>

int main()
{
  using namespace squanch;

  Model<double> model;
  PointContainerTypes<double>::map_type cp;

  auto create_circle = [&](double r) ->squanch::Nurbs<double, 1>& {
    const double pi = 3.14159265359;
    Circle<double, Eigen::Vector3d> circle;
    circle.start_angle = 0.0;
    circle.end_angle = 2 * pi;
    circle.N1 = Eigen::Vector3d(1, 0, 0);
    circle.N2 = Eigen::Vector3d(0, 1, 0);
    circle.origin = Eigen::Vector3d(0, 0, 0);
    circle.r = r;
    return circle.to_nurbs(model);
  };

  auto&& nrb_c1 = create_circle(1.0);
  auto&& nrb_c2 = create_circle(0.9);
}