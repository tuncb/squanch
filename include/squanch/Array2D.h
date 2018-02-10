#pragma once
#include <ostream>
#include <squanch/PointContainerTypes.h>

namespace squanch {

template <typename T> class PointArray2D
  {
  public:
    PointArray2D() {}
    PointArray2D(size_t nrows, size_t ncols) { this->resize(nrows, ncols); }

    size_t rows() const { return m_nrows; }
    size_t cols() const { return m_ncols; }
    const typename PointContainerTypes<T>::point_type& operator()(size_t r, size_t c) const { return m_vec[r * m_ncols + c]; }
    std::array<size_t, 2> shape() const { return {{ m_nrows, m_ncols }}; }


    typename PointContainerTypes<T>::point_type& operator()(size_t r, size_t c) { return m_vec[r * m_ncols + c]; }
    void resize(size_t nrows, size_t ncols) 
    { 
      m_nrows = nrows;
      m_ncols = ncols;
      m_vec.resize(m_nrows * m_ncols); 
    }
    template <typename ArrType> void resize(const ArrType& arr)
    {
      this->resize(arr[0], arr[1]);
    }

    friend std::ostream& operator<< (std::ostream& stream, const PointArray2D& arr) {
      for (auto i = 0u; i < arr.rows(); ++i) {
        for (auto j = 0u; j < arr.cols(); ++j) {
          auto&& p = arr(i, j);
          stream << "(" << i << "," << j << "): " << p.transpose() << "\n";
        }
      }
      return stream;
    }

  private:
    size_t m_nrows;
    size_t m_ncols;
    typename PointContainerTypes<T>::vector_type m_vec;
  };



}