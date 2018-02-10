#pragma once
#include <squanch/Model.h>
#include <squanch/nurbs.h>


namespace squanch {

  template <typename T>
  class RectangleXy
  {
  public:
    RectangleXy() {}
    RectangleXy(T width, T height) : w(width), h(height) {}

    T w;
    T h;

    Nurbs<T, 1>& to_nurbs(squanch::Model<T>& model)
    {
      auto&& nurbs = model.new_nurbs<1>();
      nurbs.curves()[0] = this->to_curve();

      int id1, id2, id3, id4;
      std::tie(id1, std::ignore) = model.new_point((T)0, (T)0, (T)0, (T)1);
      std::tie(id2, std::ignore) = model.new_point(   w, (T)0, (T)0, (T)1);
      std::tie(id3, std::ignore) = model.new_point(   w,    h, (T)0, (T)1);
      std::tie(id4, std::ignore) = model.new_point((T)0,    h, (T)0, (T)1);

      nurbs.cpi().resize(5);
      nurbs.cpi()(0) = id1;
      nurbs.cpi()(1) = id2;
      nurbs.cpi()(2) = id3;
      nurbs.cpi()(3) = id4;
      nurbs.cpi()(4) = id1;

      return nurbs;
    }

    NurbsCurve<T> to_curve()
    {
      NurbsCurve<T> curve;

      curve.set_p(1);
      curve.knots().reserve(7);
      curve.knots().push_back((T)0);
      curve.knots().push_back((T)0);
      curve.knots().push_back((T)0.25);
      curve.knots().push_back((T)0.5);
      curve.knots().push_back((T)0.75);
      curve.knots().push_back((T)1);
      curve.knots().push_back((T)1);

      curve.set_closed(true);

      return curve;
    }
  };

}