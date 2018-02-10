#include <squanch\Nurbs.h>


namespace squanch {

  template <typename T> void knit_closed_shape(squanch::Nurbs<T, 1>& nrb)
  {
    auto&& curve = nrb.curves()[0];
    if (! curve.is_closed()) return;
    auto n = curve.n();
    nrb.cpi()[n-1] = nrb.cpi()[0];
  }

  template <typename T> void knit_closed_shape(squanch::Nurbs<T, 2>& nrb)
  {
    int n1 = nrb.curves()[0].n();
    int n2 = nrb.curves()[1].n();

    if (nrb.curves()[0].is_closed()) {
      for (int i = 0; i < n2; ++i) {
        nrb.cpi()(n1-1,i) = nrb.cpi()(0,i);
      }
    }

    if (nrb.curves()[1].is_closed()) {
      for (int i = 0; i < n1; ++i) {
        nrb.cpi()(i,n2-1) = nrb.cpi()(i,0);
      }
    }
  }

  template <typename T> void knit_closed_shape(squanch::Nurbs<T, 3>& nrb)
  {
    int n1 = nrb.curves()[0].n();
    int n2 = nrb.curves()[1].n();
    int n3 = nrb.curves()[2].n();

    if (nrb.curves()[0].is_closed()) {
      for (int i = 0; i < n2; ++i) {
        for (int j = 0; j < n3; ++j) {
          nrb.cpi()(n1-1,i,j) = nrb.cpi()(0,i,j);
        }
      }
    }

    if (nrb.curves()[1].is_closed()) {
      for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n3; ++j) {
          nrb.cpi()(i,n2-1,j) = nrb.cpi()(i,0,j);
        }
      }
    }

    if (nrb.curves()[2].is_closed()) {
      for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
          nrb.cpi()(i,j,n3-1) = nrb.cpi()(i,j,0);
        }
      }
    }

  }

}