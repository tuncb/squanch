#pragma once
#include <vector>

namespace squanch {

  /// <summary>
  /// Represents a NURBS curve.
  /// </summary>
  template <typename T> class NurbsCurve
  {
  public:
    NurbsCurve() : _is_closed(false), _p(0) {}
    NurbsCurve (unsigned p) : _p(p), _is_closed(false) {}

    unsigned int p() const {return _p;}

    void set_p(unsigned int p){_p = p;}

    const std::vector<T>& knots() const {return _knots;}

    std::vector<T>& knots()       {return _knots;}
 
    size_t n() const { return _knots.size() - _p - 1; }

    void set_closed(bool b){_is_closed = b;}
    bool is_closed() const { return _is_closed;}

  private:
    unsigned int _p;
    std::vector<T> _knots;
    bool _is_closed;
  };
}