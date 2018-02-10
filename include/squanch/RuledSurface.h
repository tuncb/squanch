#pragma once
#include <squanch/Model.h>
#include <squanch/nurbs.h>
#include <squanch/refinement.h>


namespace squanch {

  namespace detail {

    template <typename T> std::vector<T> find_merging_knots(const std::vector<T>& base_knots, 
                                                            const std::vector<T>& target_knots)
    {
      std::vector<T> merging_knots;

      auto base_iter = std::begin(base_knots);
      for (auto knot : target_knots) {
        auto iter = std::find(base_iter, std::end(base_knots), knot);
        if (iter == std::end(base_knots)) {
          merging_knots.push_back(knot);
        }
        else base_iter = iter;
      }
      return merging_knots;
    }


    template <typename T, unsigned int NDIM> void p_refine_helper(squanch::Model<T>& model, Nurbs<T, NDIM>& nurbs, unsigned int curve_index, unsigned int p_delta)
    {
      p_refine(model, nurbs, curve_index, p_delta);
    }

    template <typename T> void p_refine_helper(squanch::Model<T>& model, Nurbs<T, 1>& nurbs, unsigned int, unsigned int p_delta)
    {
      p_refine(model, nurbs, p_delta);
    }


    template <typename T, unsigned int NDIM> void make_same_degree(squanch::Model<T>& model, Nurbs<T, NDIM>& nurbs1, Nurbs<T, NDIM>& nurbs2)
    {
      for (auto i = 0u; i < NDIM; ++i) {
        auto p1 = nurbs1.curves()[i].p();
        auto p2 = nurbs2.curves()[i].p();
        auto pmax = std::max(p1, p2);
        if (p1 < pmax) p_refine_helper(model, nurbs1, i, pmax - p1);
        if (p2 < pmax) p_refine_helper(model, nurbs2, i, pmax - p2);
      }
    }


    template <typename T, unsigned int NDIM> void h_refine_helper(squanch::Model<T>& model, Nurbs<T, NDIM>& nurbs, unsigned int curve_index, const std::vector<T>& knots)
    {
      h_refine(model, nurbs, curve_index, knots);
    }

    template <typename T> void h_refine_helper(squanch::Model<T>& model, Nurbs<T, 1>& nurbs, unsigned int, const std::vector<T>& knots)
    {
      h_refine(model, nurbs, knots);
    }

    template <typename T, unsigned int NDIM> void merge_knots(squanch::Model<T>& model, Nurbs<T, NDIM>& nurbs1, Nurbs<T, NDIM>& nurbs2)
    {
      for (auto i = 0u; i < NDIM; ++i) {
        auto&& u1 = nurbs1.curves()[i].knots();
        auto&& u2 = nurbs2.curves()[i].knots();

        auto merging_u1 = find_merging_knots(u1, u2);
        if (!merging_u1.empty()) h_refine_helper(model, nurbs1, i, merging_u1);

        auto merging_u2 = find_merging_knots(u2, u1);
        if (!merging_u2.empty()) h_refine_helper(model, nurbs2, i, merging_u2);
      }
    }

    template <typename T, unsigned int NDIM> void make_line_dimension(Nurbs<T, NDIM>& nurbs, unsigned int curve_index)
    {
      nurbs.curves()[curve_index].set_p(1);
      nurbs.curves()[curve_index].knots().push_back(0);
      nurbs.curves()[curve_index].knots().push_back(0);
      nurbs.curves()[curve_index].knots().push_back(1);
      nurbs.curves()[curve_index].knots().push_back(1);
    }
  }




  template <typename T> Nurbs<T, 2>& create_ruled_surface(squanch::Model<T>& model, Nurbs<T, 1>& c1, Nurbs<T, 1>& c2)
  {
    using namespace squanch;

    detail::make_same_degree(model, c1, c2);
    detail::merge_knots(model, c1, c2);
    // both curves should have the same degree
    //auto p1 = c1.curves()[0].p();
    //auto p2 = c2.curves()[0].p();
    //auto pmax = std::max(p1,p2);
    //if (p1 < pmax) p_refine(model, c1, pmax-p1);
    //if (p2 < pmax) p_refine(model, c2, pmax-p2);

    // merge knot vectors
    //auto&& u1 = c1.curves()[0].knots();
    //auto&& u2 = c2.curves()[0].knots();

    //auto merging_u1 = detail::find_merging_knots(u1, u2);
    //if (! merging_u1.empty()) h_refine(model, c1, merging_u1);
    //
    //auto merging_u2 = detail::find_merging_knots(u2, u1);
    //if (! merging_u2.empty()) h_refine(model, c2, merging_u2);

    // finalize product:
    auto&& nurbs = model.new_nurbs<2>();
    nurbs.curves()[0] = c1.curves()[0];
    detail::make_line_dimension(nurbs, 1);
    //nurbs.curves()[1].set_p(1);
    //nurbs.curves()[1].knots().push_back(0);
    //nurbs.curves()[1].knots().push_back(0);
    //nurbs.curves()[1].knots().push_back(1);
    //nurbs.curves()[1].knots().push_back(1);

    auto n = c1.cpi().size();
    nurbs.cpi().resize(n, 2);

    for (unsigned int i = 0; i < n; ++i) {
      nurbs.cpi()(blitz::TinyVector<int, 2>(i,0)) = c1.cpi()(i);
      nurbs.cpi()(blitz::TinyVector<int, 2>(i,1)) = c2.cpi()(i);
    }
    
    return nurbs;
  }
  
  template <typename T> Nurbs<T, 3>& create_ruled_volume(squanch::Model<T>& model, Nurbs<T, 2>& surface1, Nurbs<T, 2>& surface2)
  {
    using namespace squanch;

    detail::make_same_degree(model, surface1, surface2);
    detail::merge_knots(model, surface1, surface2);

    auto&& nurbs = model.new_nurbs<3>();
    nurbs.curves()[0] = surface1.curves()[0];
    nurbs.curves()[1] = surface1.curves()[1];

    detail::make_line_dimension(nurbs, 2);

    auto n = surface2.cpi().size();
    nurbs.cpi().resize(surface1.cpi().rows(), surface1.cpi().cols(), 2);

    for (unsigned int i = 0; i < surface1.cpi().rows(); ++i) {
      for (unsigned int j = 0; j < surface1.cpi().cols(); ++j) {
        nurbs.cpi()(blitz::TinyVector<int, 3>(i, j, 0)) = surface1.cpi()(blitz::TinyVector<int, 2>(i, j));
        nurbs.cpi()(blitz::TinyVector<int, 3>(i, j, 1)) = surface2.cpi()(blitz::TinyVector<int, 2>(i, j));
      }
    }

    return nurbs;
  }

}