#pragma once
#include <squanch/Model.h>
#include <squanch/nurbs.h>


namespace squanch {

  template <typename T, typename PointType>
  class Line
  {
  public:
    Line() {}
    Line(PointType p1, PointType p2) : p1(p1), p2(p2) {}

    PointType p1;
    PointType p2;

    Nurbs<T, 1>& to_nurbs(squanch::Model<T>& model)
    {
      auto&& nurbs = model.new_nurbs<1>();
      nurbs.curves()[0] = this->to_curve();

      int id1, id2;

      std::tie(id1, std::ignore) = model.new_point(p1.x(), p1.y(), p1.z(), 1.0);
      std::tie(id2, std::ignore) = model.new_point(p2.x(), p2.y(), p2.z(), 1.0);

      nurbs.cpi().resize(2);
      nurbs.cpi()(0) = id1;
      nurbs.cpi()(1) = id2;

      return nurbs;
    }

    NurbsCurve<T> to_curve()
    {
      NurbsCurve<T> curve;

      curve.set_p(1);
      curve.knots().reserve(4);
      curve.knots().push_back((T)0);
      curve.knots().push_back((T)0);
      curve.knots().push_back((T)1);
      curve.knots().push_back((T)1);

      return curve;
    }
  };

}