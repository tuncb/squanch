#pragma once
#include <vector>

namespace squanch {

  /// <summary>
  /// Represents a NURBS curve.
  /// </summary>
  template <typename T> class NurbsCurve
  {
  public:
    /// <summary>
    /// Initializes a new instance of the <see cref="NurbsCurve{T}"/> class.
    /// </summary>
    NurbsCurve() : _is_closed(false), _p(0) {}
    /// <summary>
    /// Initializes a new instance of the <see cref="NurbsCurve{T}"/> class.
    /// </summary>
    /// <param name="p">The degree of the curve</param>
    NurbsCurve (unsigned p) : _p(p), _is_closed(false) {}

    /// <summary>
    /// Initializes a new instance of the <see cref="NurbsCurve{T}"/> class.
    /// </summary>
    /// <param name="other">The other NURBS curve to be moved.</param>
    NurbsCurve(NurbsCurve&& other)
    {
      _p = other._p;
      _knots = std::move(other._knots);
      _is_closed = other._is_closed;
    }

    /// <summary>
    /// Move assignment operator.
    /// </summary>
    /// <param name="other">The other NURBS curve to be moved.</param>
    /// <returns></returns>
    NurbsCurve& operator=(NurbsCurve&& other)
    {
      if (this != &other) {
        _p = other._p;
        _knots = std::move(other._knots);
        _is_closed = other._is_closed;
      }
      return *this;
    }

    /// <summary>
    /// Gives the degree of the curve.
    /// </summary>
    /// <returns>Degree of the curve</returns>
    unsigned int p() const {return _p;}

    /// <summary>
    /// Sets the degree of the curve the specified.
    /// </summary>
    /// <param name="p">The degree of the curve.</param>
    void set_p(unsigned int p){_p = p;}

    /// <summary>
    /// Gives the knots of the curve.
    /// </summary>
    /// <returns>The knots of the curve</returns>
    const std::vector<T>& knots() const {return _knots;}
    /// <summary>
    /// Gives the knots of the curve.
    /// </summary>
    /// <returns>The knots of the curve</returns>>
    std::vector<T>& knots()       {return _knots;}
    
    /// <summary>
    /// Computes the total number of basis functions of the curve.
    /// </summary>
    /// <returns>Total number of basis functions of the curve</returns>
    size_t n() const { return _knots.size() - _p - 1; }

    void set_closed(bool b){_is_closed = b;}
    bool is_closed() const { return _is_closed;}

  private:
    unsigned int _p;
    std::vector<T> _knots;
    bool _is_closed;
  };

}