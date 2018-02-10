#pragma once
#include <vector>
#include <array>
#include <Eigen/core>
#include <blitz/array.h>
#include <squanch/NurbsCurve.h>
#include <squanch/PatchBase.h>

namespace squanch {
  template <typename T> class Model;

 namespace detail {

   /// <summary>
   /// Represents a NURBS geometry composed of Dim NurbsCurves. T is basic value type (generally float or double)
   /// </summary>
  template <typename T, unsigned int Dim> class NurbsBase : public squanch::PatchBase
  {
  public:
    /// <summary>
    /// Initializes a new instance of the <see cref="NurbsBase{T, Dim}"/> class.
    /// </summary>
    NurbsBase(int id, squanch::PatchType pt) : PatchBase(id, pt) {}
    /// <summary>
    /// Copy Constructor
    /// </summary>
    NurbsBase(const NurbsBase<T, Dim>& other) : PatchBase(other.id(), other.type())
    {
      auto&& shape = other._cpi.shape();
      _cpi.resize(other._cpi.shape());
      _cpi = other._cpi;
      _curves = other._curves;
    }
    /// <summary>
    /// Move Constructor
    /// </summary>
    NurbsBase(NurbsBase<T, Dim>&& other) : PatchBase(other.id(), other.type())
    {
      _curves = std::move(other._curves);
      _cpi.reference(other._cpi); // reference to other cpi's memory space
    }

    /// <summary>
    /// Copy assignment operator.
    /// </summary>
    NurbsBase<T, Dim>& operator=(NurbsBase<T, Dim>& other)
    {
      if (this != &other) {
        _curves = other._curves;
        _cpi.resize(other._cpi.shape());
        _cpi = other._cpi;
        PatchBase::operator=(other);
      }
      return *this;
    }
    /// <summary>
    /// Move assignment operator.
    /// </summary>
    NurbsBase<T, Dim>& operator=(NurbsBase<T, Dim>&& other)
    {
      if (this != &other) {
        _curves = std::move(other._curves);
        _cpi.reference(other._cpi); // reference to other cpi's memory space
        PatchBase::operator=(other);
      }
      return *this;
    }

    typedef T value_type;

    /// <summary>
    /// Gives the curves of the NURBS geometry.
    /// </summary>
    const std::array<squanch::NurbsCurve<T>, Dim>& curves() const {return _curves;}
    /// <summary>
    /// Gives the curves of the NURBS geometry.
    /// </summary>
    std::array<squanch::NurbsCurve<T>, Dim>& curves()       {return _curves;}

    /// <summary>
    /// Returns the control point index of the geometry.
    /// </summary>
    /// <returns>The control point index of the geometry</returns>
    const blitz::Array<int, Dim> &  cpi() const { return _cpi; }
    /// <summary>
    /// Returns the control point index of the geometry.
    /// </summary>
    /// <returns>The control point index of the geometry</returns>
    blitz::Array<int, Dim> &  cpi()       { return _cpi; }
    
    /// <summary>
    /// Returns the number of dimensions of the NURBS geometry.
    /// </summary>
    /// <returns>the number of dimensions of the NURBS geometry</returns>
    unsigned int ndim() const {return Dim;}

    size_t num_cp() const
    {
      unsigned int num_cp = 1;
      for (unsigned i = 0; i < Dim; ++i) {
        num_cp *= _curves[i].p() + 1;
      }
      return num_cp;
    }
  private:
    std::array<squanch::NurbsCurve<T>, Dim> _curves;
    blitz::Array<int, Dim> _cpi;
  };
}


  template <template<class, unsigned int> class NurbsType, typename T, unsigned int NDIM> 
  std::vector<int> get_node_ids(const NurbsType<T, NDIM>& spline)
  {
    return std::vector<int>(spline.cpi().begin(), spline.cpi().end());
  }

  template <typename CurveArray> unsigned int num_nodes_per_span(const CurveArray& curves)
  {
    auto num = 1u;
    const auto dim = curves.size();
    for (auto i = 0u; i < dim; ++i) {
      num *= (curves[i].p() + 1);
    }
    return num;
  }

  inline unsigned int find_other_dimension(unsigned int dim1)
  {
    if (dim1 == 0) return 1; else return 0;
    throw;
  }

  inline unsigned int find_other_dimension(unsigned int dim1, unsigned int dim2)
  {
    if (dim1 == 0 && dim2 == 1) return 2;
    if (dim1 == 1 && dim2 == 0) return 2;
    if (dim1 == 0 && dim2 == 2) return 1;
    if (dim1 == 2 && dim2 == 0) return 1;
    if (dim1 == 1 && dim2 == 2) return 0;
    if (dim1 == 2 && dim2 == 1) return 0;

    throw;
  }

  inline std::array<unsigned int, 2> find_other_dimensions(unsigned int dim)
  {
    if (dim == 0) return{ { 1, 2 } };
    if (dim == 1) return{ { 0, 2 } };
    if (dim == 2) return{ { 0, 1 } };
    throw;
  }

  template <typename T, unsigned int Dim> class Nurbs : public detail::NurbsBase<T, Dim> {};

  template <typename T> class Nurbs<T, 1> : public detail::NurbsBase<T, 1>
  {
  public:
    typedef T ParametricPointType;
    typedef unsigned int SpanType;

    Nurbs<T, 1>& operator=(Nurbs<T, 1>& other)  { NurbsBase<T, 1>::operator=(other); return *this; }
    
    Nurbs(const Nurbs<T, 1>& nurbs) = delete;
    Nurbs(Nurbs<T, 1>&& nurbs)      = delete;
    Nurbs<T, 1>& operator=(Nurbs<T, 1>&& other) = delete;
  private:
    friend class squanch::Model<T>;
    Nurbs(int id) : detail::NurbsBase<T, 1>(id, squanch::PatchType::nurbs1d) {}
  };

  template <typename T> class Nurbs<T, 2> : public detail::NurbsBase<T, 2>
  {
  public:
    typedef std::array<T, 2> ParametricPointType;
    typedef std::array<unsigned int, 2> SpanType;

    Nurbs<T, 2>& operator=(Nurbs<T, 2>& other)  { NurbsBase<T, 2>::operator=(other); return *this; }
    
    Nurbs(const Nurbs<T, 2>& nurbs) = delete;
    Nurbs(Nurbs<T, 2>&& nurbs)      = delete;
    Nurbs<T, 2>& operator=(Nurbs<T, 2>&& other) = delete;
  private:
    friend class squanch::Model<T>;
    Nurbs(int id) : detail::NurbsBase<T, 2>(id, squanch::PatchType::nurbs2d) {}
  };

  template <typename T> class Nurbs<T, 3> : public detail::NurbsBase<T, 3>
  {
  public:
    typedef std::array<T, 3> ParametricPointType;
    typedef std::array<unsigned int, 3> SpanType;

    Nurbs<T, 3>& operator=(Nurbs<T, 3>& other)  { NurbsBase<T, 3>::operator=(other); return *this; }
    
    Nurbs(const Nurbs<T, 3>& nurbs) = delete;
    Nurbs(Nurbs<T, 3>&& nurbs)      = delete;
    Nurbs<T, 3>& operator=(Nurbs<T, 3>&& other) = delete;
  private:
    friend class squanch::Model<T>;
    Nurbs(int id) : detail::NurbsBase<T, 3>(id, squanch::PatchType::nurbs3d) {}
  };

  enum class LineEnds : unsigned int
  {
    Start = 0, End = 1
  };

  template <typename T, unsigned int Dim> class SkeletonCurvesContainer
  {
  public:
    template <unsigned int ParentDim> SkeletonCurvesContainer(std::array<squanch::NurbsCurve<T>, ParentDim>& parent_curves, std::array<unsigned int, Dim> curve_positions)
    {
      auto c = 0u;
      for (auto pos : curve_positions) {
        _curves[c++] = &(parent_curves[pos]);
      }
    }

    size_t size() const { return Dim; }

    const squanch::NurbsCurve<T>& operator[](size_t pos) const { return *(_curves[pos]); }
    squanch::NurbsCurve<T>& operator[](size_t pos) { return *(_curves[pos]); }
  private:
    std::array<squanch::NurbsCurve<T>*, Dim> _curves;
  };


  template <typename T, unsigned int Dim> class NurbsSkeletonBase
  {
  public:
    NurbsSkeletonBase() {}
    NurbsSkeletonBase(SkeletonCurvesContainer<T, Dim>& curves) : _curves(curves) {}
    NurbsSkeletonBase(SkeletonCurvesContainer<T, Dim>& curves, blitz::Array<int, Dim>& parent_cpi) : _curves(curves) 
    {
      _cpi = parent_cpi; // use assignment to share memory
    }


    const blitz::Array<int, Dim> &  cpi() const { return _cpi; }
    const SkeletonCurvesContainer<T, Dim>& curves() const { return _curves; }

    blitz::Array<int, Dim> &  cpi()       { return _cpi; }
    SkeletonCurvesContainer<T, Dim>& curves()       { return _curves; }

    unsigned int ndim() const { return Dim; }

  protected:
    SkeletonCurvesContainer<T, Dim> _curves;
    blitz::Array<int, Dim> _cpi; // this will be a reference array
  };


  template <typename T, unsigned int Dim> class NurbsSkeleton : public NurbsSkeletonBase<T, Dim> {};

  template <typename T> class NurbsSkeleton<T, 3> : public NurbsSkeletonBase<T, 3>
  {
  public:
    typedef std::array<T, 3> ParametricPointType;
    typedef std::array<unsigned int, 3> SpanType;
    template <template<class, unsigned int> class NurbsType>
    NurbsSkeleton(NurbsType<T, 3>& parent) : NurbsSkeletonBase<T, 3>(SkeletonCurvesContainer<T, 3>(parent.curves(), { { 0, 1, 2 } }), parent.cpi()) { }
  };

  template <typename T> class NurbsSkeleton<T, 2> : public NurbsSkeletonBase<T, 2>
  {
  public:
    typedef std::array<T, 2> ParametricPointType;
    typedef std::array<unsigned int, 2> SpanType;
    
    template <template<class, unsigned int> class NurbsType>
    NurbsSkeleton(NurbsType<T, 2>& parent) : NurbsSkeletonBase<T, 2>(SkeletonCurvesContainer<T, 2>(parent.curves(), { { 0, 1 } }), parent.cpi()) { }

    template <template<class, unsigned int> class NurbsType>
    NurbsSkeleton(NurbsType<T, 3>& parent, std::array<unsigned, 2> subdimensions, LineEnds otherdim_position) : NurbsSkeletonBase<T, 2>(SkeletonCurvesContainer<T, 2>(parent.curves(), subdimensions))
    {
      auto&& shape = parent.cpi().shape();
      auto all = blitz::Range::all();

      auto getCollapsedDimCpi = [otherdim_position, &shape](int dim){
        if (otherdim_position == LineEnds::Start) return 0;
        else return shape[dim] - 1;
      };
      
      auto other_dim = find_other_dimension(subdimensions[0], subdimensions[1]);
      if (other_dim == 0) {
        _cpi.reference(parent.cpi()(getCollapsedDimCpi(0), all, all));
      }
      else if (other_dim == 1) {
        _cpi.reference(parent.cpi()(all, getCollapsedDimCpi(1), all));
      }
      else {
        _cpi.reference(parent.cpi()(all, all, getCollapsedDimCpi(2)));
      }
    }
  };

  template <typename T> class NurbsSkeleton<T, 1> : public NurbsSkeletonBase<T, 1>
  {
  public:
    typedef T ParametricPointType;
    typedef unsigned int SpanType;

    template <template<class, unsigned int> class NurbsType>
    NurbsSkeleton(NurbsType<T, 1>& parent) : NurbsSkeletonBase<T, 1>(SkeletonCurvesContainer<T, 1>(parent.curves(), { { 0 } } ), parent.cpi()) { }

    template <template<class, unsigned int> class NurbsType>
    NurbsSkeleton(NurbsType<T, 2>& parent, unsigned int dim, LineEnds otherdim_position) : NurbsSkeletonBase<T, 1>(SkeletonCurvesContainer<T, 1>(parent.curves(), { { dim } }))
    {
      auto&& shape = parent.cpi().shape();
      auto all = blitz::Range::all();

      auto getCollapsedDimCpi = [otherdim_position, &shape](int dim){
        if (otherdim_position == LineEnds::Start) return 0;
        else return shape[dim] - 1;
      };

      if (dim == 0) {
        _cpi.reference(parent.cpi()(all, getCollapsedDimCpi(1)));
      }
      else {
        _cpi.reference(parent.cpi()(getCollapsedDimCpi(0), all));
      }
    }

    template <template<class, unsigned int> class NurbsType>
    NurbsSkeleton(NurbsType<T, 3>& parent, unsigned int dim, std::array<LineEnds, 2> otherdim_positions) : NurbsSkeletonBase<T, 1>(SkeletonCurvesContainer<T, 1>(parent.curves(), { { dim } }))
    {
      // there must be a better way of doing this, but I couldn't find it!
      auto&& shape = parent.cpi().shape();
      auto all = blitz::Range::all();
      
      auto getCollapsedDimCpi = [otherdim_positions, &shape](unsigned int od, int dim){
        if (otherdim_positions[od] == LineEnds::Start) return 0;
        else return shape[dim] - 1;
      };

      if (dim == 0) {
        _cpi.reference(parent.cpi()(all, getCollapsedDimCpi(0, 1), getCollapsedDimCpi(1, 2)));
      }
      else if (dim == 1) {
        _cpi.reference(parent.cpi()(getCollapsedDimCpi(0, 0), all, getCollapsedDimCpi(1, 2)));
      }
      else {
        _cpi.reference(parent.cpi()(getCollapsedDimCpi(0, 0), getCollapsedDimCpi(1, 1), all));
      }
    }
  };



}