#pragma once
#include <Eigen/geometry>
#include <squanch/Model.h>
#include <squanch/nurbs.h>
#include <squanch/Line.h>


namespace squanch {

  namespace detail {

    template <typename T> NurbsCurve<T> create_line_curve()
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
  }


  /// <summary>
  /// Extrudes the specified in.
  /// </summary>
  /// <param name="in">The in.</param>
  /// <param name="vec">The vec.</param>
  /// <param name="cp">The cp.</param>
  /// <param name="id_s">The id_s.</param>
  /// <returns></returns>
  template <typename T> Nurbs<T, 2>& extrude(squanch::Model<T>& model, const Nurbs<T, 1>& in, const typename PointContainerTypes<T>::point_type& vec)
  {
    PointContainerTypes<T>::point_type ext = vec;
    ext.w() = 0;
    auto&& out = model.new_nurbs<2>();
    out.curves()[0] = in.curves()[0];
    out.curves()[1] = detail::create_line_curve<T>();
    auto n = in.curves()[0].n();
    out.cpi().resize(n, 2);

    auto end = n;
    if (in.curves()[0].is_closed()) --end;

    for (int i = 0; i < end; ++i) {
      auto id = in.cpi()(i);

      auto pt_tuple = model.new_point();
      auto id_pt = std::get<0>(pt_tuple);
      auto&& pt = std::get<1>(pt_tuple);
      
			auto&& p = model.get_point(id);

      pt = p + ext * p.w();
      out.cpi()(i,0) = id;
      out.cpi()(i,1) = id_pt;
    }

    if (in.curves()[0].is_closed()) {
      out.cpi()((int)(n-1),0) = in.cpi()(0);
      out.cpi()(int(n-1),1) = out.cpi()(0,1);
    }

    return out;
  }

  /// <summary>
  /// Extrudes the specified in.
  /// </summary>
  /// <param name="in">The in.</param>
  /// <param name="vec">The vec.</param>
  /// <param name="cp">The cp.</param>
  /// <param name="id_s">The id_s.</param>
  /// <returns></returns>
  template <typename T> Nurbs<T, 3>& extrude(squanch::Model<T>& model, const Nurbs<T, 2>& in, const typename PointContainerTypes<T>::point_type& vec)
  {
    PointContainerTypes<T>::point_type ext = vec;
    ext.w() = 0;
    auto&& out = model.new_nurbs<3>();
    out.curves()[0] = in.curves()[0];
    out.curves()[1] = in.curves()[1];
    out.curves()[2] = detail::create_line_curve<T>();
    out.cpi().resize(in.curves()[0].n(), in.curves()[1].n(), 2);

    std::array<size_t, 2> ends;
    ends[0] = in.curves()[0].n();
    ends[1] = in.curves()[1].n();
    for (unsigned int i = 0; i < 2; ++i) {
      if (in.curves()[i].is_closed()) ends[i] -= 1;
    }

    blitz::TinyVector<int, 2> index;

    for (index[0] = 0; index[0] < ends[0]; ++index[0]) {
      for (index[1] = 0; index[1] < ends[1]; ++index[1]) {
        auto id = in.cpi()(index);

        auto pt_tuple = model.new_point();
        auto id_pt = std::get<0>(pt_tuple);
        auto&& pt = std::get<1>(pt_tuple);
        
				auto&& p = model.get_point(id);

        pt = p + ext * p.w();
        out.cpi()(index[0], index[1], 0) = id;
        out.cpi()(index[0], index[1], 1) = id_pt;
      }
    }

    if (in.curves()[0].is_closed()) {
      for (int i = 0; i < in.curves()[1].n(); ++i) {
        out.cpi()((int)in.curves()[0].n()-1,i,0) = out.cpi()(0,i,0);
        out.cpi()((int)in.curves()[0].n()-1,i,1) = out.cpi()(0,i,1);
      }
    }

    if (in.curves()[1].is_closed()) {
      for (int i = 0; i < in.curves()[0].n(); ++i) {
        out.cpi()(i,(int)in.curves()[1].n()-1,0) = out.cpi()(i,0,0);
        out.cpi()(i,(int)in.curves()[1].n()-1,1) = out.cpi()(i,0,1);
      }
    }


    return out;
  }

}