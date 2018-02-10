#pragma once
#include <boost/multi_array.hpp>
#include <vector>
#include <squanch/Nurbs.h>
#include <squanch/Model.h>
#include <squanch/IgaAlgorithms.h>
#include <squanch/knit.h>
#include <squanch/UtilityDetail.h>

namespace squanch {

namespace detail {



  template <typename T> std::vector<T> find_new_knots(const std::vector<T>& knots, unsigned int factor)
  {
    if (factor <= 1) return std::vector<T>();

    std::vector<T> new_knots;
    const auto s = knots.size();
    for (size_t i = 1; i < s; ++i) {
      auto delta = knots[i] - knots[i - 1];
      if ( delta > (T)0.0000001) { // if knot span > 0
        float du = delta / factor;
        for (unsigned int j = 1; j < factor; ++j) {
          new_knots.push_back(knots[i - 1] + j*du);
        }
      }
    }
    return new_knots;
  }

}

  template <typename T> void h_refine(squanch::Model<T>& model, Nurbs<T, 1>& nurbs, const std::vector<T>& x)
  {
    /* algorithm A5.4 from the NURBS book*/
    if (x.empty()) return;

    auto&& curve = nurbs.curves()[0];
    auto&& u = curve.knots();
    auto p = curve.p();
    auto n = curve.n()-1;
    auto m = n + p + 1;
    auto r = x.size()-1;

    blitz::TinyVector<int, 1> arr; arr[0] = n + r + 2;
    nurbs.cpi().resizeAndPreserve(arr);

    // find start and end spans effected for knot insertion.
    auto a = squanch::find_span(curve.knots(), curve.p(), x.front());
    auto b = squanch::find_span(curve.knots(), curve.p(), x.back());
    b += 1;

    // ubar is the new knot vector.
    // copy unaffected knots from the beginning.
    std::vector<T> ubar(u.size() + x.size());
    for (size_t j = 0; j <= a; ++j) ubar[j] = u[j];
    // copy unaffected knots from the end.
    for (size_t j = b+p; j <= m; ++j) { ubar[j+r+1] = u[j]; }
    
    // copy unaffected control points from the end 
    for (size_t j = b - 1; j <= n; ++j) { nurbs.cpi()(j + r + 1) = nurbs.cpi()(j); }

    // insert new knots and control points
    auto i = b+p-1;
    auto k = b+p+r;
    for (int j = r; j>=0; --j) {
      while(i > a && x[j] <= u[i]) {
        auto&& old_point = model.get_point(nurbs.cpi()(i - p - 1));
        auto&& tup = model.new_point(old_point.x(), old_point.y(), old_point.z(), old_point.w());
        nurbs.cpi()(k - p - 1) = std::get<0>(tup);
        ubar[k] = u[i];
        k--; i--;
      }
      nurbs.cpi()(k-p-1) = nurbs.cpi()(k-p);
      for (unsigned int l = 1; l <= p; ++l) {
        auto ind = k-p+l;
        auto alpha = ubar[k+l]-x[j];
        if (abs(alpha) == 0.0) 
          nurbs.cpi()(ind-1) = nurbs.cpi()(ind);
        else {
          alpha = alpha / (ubar[k+l]-u[i-p+l]);
          auto&& tup = model.new_point();
          auto&& point = std::get<1>(tup);
          point = PointContainerTypes<T>::point_type( alpha * (model.get_point(nurbs.cpi()(ind-1))) + ((T)1 - alpha) * (model.get_point(nurbs.cpi()(ind))) ) ;
          nurbs.cpi()(ind-1) = std::get<0>(tup);
        }
      }
      ubar[k] = x[j];
      --k;
    }
    u = ubar;
  }

  template <typename T> void h_refine(squanch::Model<T>& model, Nurbs<T, 2>& nurbs, unsigned int curve_i, const std::vector<T>& x)
  {
    /* algorithm A5.4 from the NURBS book*/
    if (x.empty()) return;

    // same stuff from the one dimensional version
    auto&& curve = nurbs.curves()[curve_i];
    auto&& u = curve.knots();
    int p = curve.p();
    auto n = curve.n()-1;
    auto m = n + p + 1;
    int r = x.size()-1;

    //auto id_s_s = id_s;

    auto arr = nurbs.cpi().shape();
    arr[curve_i] = n + r + 2;
    nurbs.cpi().resizeAndPreserve(arr);

    auto a = squanch::find_span(curve.knots(), curve.p(), x.front());
    auto b = squanch::find_span(curve.knots(), curve.p(), x.back());
    b += 1;

    std::vector<T> ubar(u.size() + x.size());
    for (size_t j = 0; j <= a; ++j) ubar[j] = u[j];
    for (size_t j = b+p; j <= m; ++j) { ubar[j+r+1] = u[j]; }
    // 2D set point function
    auto set_point = [&model, &nurbs](size_t pos, size_t nr, typename PointContainerTypes<T>::point_type& p){
      int id;
      std::tie(id, std::ignore) = model.new_point(p.x(), p.y(), p.z(), p.w());
      nurbs.cpi()((int)pos,(int)nr) = id;
    };
		auto set_point_w_id = [&model, &nurbs](size_t pos, size_t nr, int id){
			auto&& pair = model.new_point();
			pair.second = model.get_point(id);
			nurbs.cpi()((int)pos, (int)nr) = pair.first;
		};

    if (curve_i == 0) {
      auto num_row = arr[1];
      for (int i = 0; i < num_row; ++i) {
        for (int j = b-1; j <= n; ++j) nurbs.cpi()((int)(j+r+1),i) = nurbs.cpi()(j,i);
      }

      auto i = b+p-1;
      int k = b+p+r;
      for (int j = r; j>=0; --j) {
        while(i > a && x[j] <= u[i]) {
          for (int nr = 0; nr < num_row; ++nr) set_point_w_id(k-p-1, nr, (nurbs.cpi()((int)i-p-1,nr)));
          ubar[k] = u[i];
          k--; i--;
        }
        for (int nr = 0; nr < num_row; ++nr) nurbs.cpi()(k-p-1,nr) = nurbs.cpi()(k-p,nr);
        for (unsigned int l = 1; l <= p; ++l) {
          int ind = k-p+l;
          auto alpha = ubar[k+l]-x[j];
          if (abs(alpha) == 0.0) 
            for (int nr = 0; nr < num_row; ++nr) nurbs.cpi()(ind-1,nr) = nurbs.cpi()(ind,nr);
          else {
            alpha = alpha / (ubar[k+l]-u[i-p+l]);
            for (int nr = 0; nr < num_row; ++nr){
              PointContainerTypes<T>::point_type point = PointContainerTypes<T>::point_type(alpha*(model.get_point(nurbs.cpi()(ind - 1, nr))) + ((T)1 - alpha)*(model.get_point(nurbs.cpi()(ind, nr))));
              set_point(ind-1, nr, point);
            }
          }
        }
        ubar[k] = x[j];
        --k;
      }
    } else {
      auto num_col = arr[0];
      for (int i = 0; i < num_col; ++i) {
        for (int j = b-1; j <= n; ++j) nurbs.cpi()(i,j+r+1) = nurbs.cpi()(i,j);
      }

      int i = b+p-1;
      int k = b+p+r;
      for (int j = r; j>=0; --j) {
        while(i > a && x[j] <= u[i]) {
          for (int nr = 0; nr < num_col; ++nr) set_point_w_id(nr, k - p - 1, (nurbs.cpi()(nr, i - p - 1)));
          ubar[k] = u[i];
          k--; i--;
        }
        for (int nr = 0; nr < num_col; ++nr) nurbs.cpi()(nr,k-p-1) = nurbs.cpi()(nr,k-p);
        for (unsigned int l = 1; l <= p; ++l) {
          int ind = k-p+l;
          auto alpha = ubar[k+l]-x[j];
          if (abs(alpha) == 0.0) 
            for (int nr = 0; nr < num_col; ++nr) nurbs.cpi()(nr,ind-1) = nurbs.cpi()(nr,ind);
          else {
            alpha = alpha / (ubar[k+l]-u[i-p+l]);
            for (int nr = 0; nr < num_col; ++nr) {
              PointContainerTypes<T>::point_type point = PointContainerTypes<T>::point_type(alpha*(model.get_point(nurbs.cpi()(nr, ind - 1))) + ((T)1 - alpha)*(model.get_point(nurbs.cpi()(nr, ind))));
              set_point(nr, ind-1, point);
            }
          }
        }
        ubar[k] = x[j];
        --k;
      }
    }

    u = ubar;
  }

  template <typename T> void h_refine(squanch::Model<T>& model, Nurbs<T, 3>& nurbs, unsigned int curve_i, const std::vector<T>& x)
  {
    /* algorithm A5.4 from the NURBS book*/
    if (x.empty()) return;

    // same stuff from the one dimensional version
    auto&& curve = nurbs.curves()[curve_i];
    auto&& u = curve.knots();
    int p = curve.p();
    auto n = curve.n()-1;
    auto m = n + p + 1;
    auto r = x.size()-1;

    //auto id_s_s = id_s;

    auto arr = nurbs.cpi().shape();
    arr[curve_i] = n + r + 2;
    nurbs.cpi().resizeAndPreserve(arr);

    auto a = squanch::find_span(curve.knots(), curve.p(), x.front());
    auto b = squanch::find_span(curve.knots(), curve.p(), x.back());
    b += 1;

    std::vector<T> ubar(u.size() + x.size());
    for (size_t j = 0; j <= a; ++j) ubar[j] = u[j];
    for (size_t j = b+p; j <= m; ++j) { ubar[j+r+1] = u[j]; }
    
    // 3D set point function
    auto set_point = [&model, &nurbs](size_t p1, size_t p2, size_t p3, typename PointContainerTypes<T>::point_type& p){
      auto&& tup = model.new_point(p.x(), p.y(), p.z(), p.w());
      nurbs.cpi()((int)p1, (int)p2, (int)p3) = std::get<0>(tup);
    };
		auto set_point_w_id = [&model, &nurbs](size_t p1, size_t p2, size_t p3, int id){
			auto&& tup = model.new_point();
			tup.second = model.get_point(id);
			nurbs.cpi()((int)p1, (int)p2, (int)p3) = std::get<0>(tup);
		};

    if (curve_i == 0) {
      auto num_row = arr[1];
      auto num_col = arr[2];
      
      for (int i = 0; i < num_row; ++i) { 
        for (int i2 = 0; i2 < num_col; ++i2) {
          for (int j = b-1; j <= n; ++j) nurbs.cpi()((int)(j+r+1),i,i2) = nurbs.cpi()(j,i,i2);
        }
      }

      auto i = b+p-1;
      auto k = b+p+r;
      for (int j = r; j>=0; --j) {
        while(i > a && x[j] <= u[i]) {
          for (int nr = 0; nr < num_row; ++nr) for (int nc = 0; nc < num_col; ++nc) set_point_w_id(k - p - 1, nr, nc, (nurbs.cpi()((int)(i - p - 1), nr, nc)));
          ubar[k] = u[i];
          k--; i--;
        }

        for (int nr = 0; nr < num_row; ++nr) for(int nc = 0; nc < num_col; ++nc) nurbs.cpi()((int)(k-p-1),nr,nc) = nurbs.cpi()((int)(k-p),nr,nc);

        for (unsigned int l = 1; l <= p; ++l) {
          int ind = k-p+l;
          auto alpha = ubar[k+l]-x[j];
          if (abs(alpha) == 0.0) 
            for (int nr = 0; nr < num_row; ++nr) for (int nc = 0; nc < num_col; ++nc) nurbs.cpi()(ind - 1, nr, nc) = nurbs.cpi()(ind, nr, nc);
          else {
            alpha = alpha / (ubar[k+l]-u[i-p+l]);
            for (int nr = 0; nr < num_row; ++nr) for(int nc = 0; nc < num_col; ++nc)  {{
              PointContainerTypes<T>::point_type point = PointContainerTypes<T>::point_type(alpha*(model.get_point(nurbs.cpi()(ind-1,nr,nc))) + ((T)1 - alpha)*(model.get_point(nurbs.cpi()(ind,nr,nc))));
              set_point(ind-1, nr, nc, point);
            }}
          }
        }
        ubar[k] = x[j];
        --k;
      }
    } else if (curve_i == 1) {
      auto num_row = arr[0];
      auto num_col = arr[2];

      for (int i = 0; i < num_row; ++i) { 
        for (int i2 = 0; i2 < num_col; ++i2) {
          for (int j = b-1; j <= n; ++j) nurbs.cpi()(i, (int)(j+r+1),i2) = nurbs.cpi()(i,j,i2);
        }
      }

      int i = b+p-1;
      int k = b+p+r;
      for (int j = r; j>=0; --j) {
        while(i > a && x[j] <= u[i]) {
          for (int nr = 0; nr < num_row; ++nr) for(int nc = 0; nc < num_col; ++nc) set_point_w_id(nr, k-p-1, nc, (nurbs.cpi()(nr,(int)(i-p-1),nc)));
          ubar[k] = u[i];
          k--; i--;
        }

        for (int nr = 0; nr < num_row; ++nr) for(int nc = 0; nc < num_col; ++nc) nurbs.cpi()(nr,(int)(k-p-1),nc) = nurbs.cpi()(nr,(int)(k-p),nc);

        for (unsigned int l = 1; l <= p; ++l) {
          int ind = k-p+l;
          auto alpha = ubar[k+l]-x[j];
          if (abs(alpha) == 0.0) 
            for (int nr = 0; nr < num_row; ++nr) for(int nc = 0; nc < num_col; ++nc) nurbs.cpi()(nr,ind-1,nc) = nurbs.cpi()(nr,ind,nc);
          else {
            alpha = alpha / (ubar[k+l]-u[i-p+l]);
            for (int nr = 0; nr < num_row; ++nr) for(int nc = 0; nc < num_col; ++nc)  {{
              PointContainerTypes<T>::point_type point = PointContainerTypes<T>::point_type(alpha*(model.get_point(nurbs.cpi()(nr,ind-1,nc))) + ((T)1 - alpha)*(model.get_point(nurbs.cpi()(nr,ind,nc))));
              set_point(nr, ind-1, nc, point);
            }}
          }
        }
        ubar[k] = x[j];
        --k;
      }
    } else if (curve_i == 2) {
      auto num_row = arr[0];
      auto num_col = arr[1];

      for (int i = 0; i < num_row; ++i) { 
        for (int i2 = 0; i2 < num_col; ++i2) {
          for (int j = b-1; j <= n; ++j) nurbs.cpi()(i,i2,(int)(j+r+1)) = nurbs.cpi()(i,i2,j);
        }
      }

      auto i = b+p-1;
      int k = b+p+r;
      for (int j = r; j>=0; --j) {
        while(i > a && x[j] <= u[i]) {
          for (int nr = 0; nr < num_row; ++nr) for(int nc = 0; nc < num_col; ++nc) set_point_w_id(nr, nc, k-p-1, (nurbs.cpi()(nr,nc,(int)(i-p-1))));
          ubar[k] = u[i];
          k--; i--;
        }

        for (int nr = 0; nr < num_row; ++nr) for(int nc = 0; nc < num_col; ++nc) nurbs.cpi()(nr,nc,k-p-1) = nurbs.cpi()(nr,nc,k-p);

        for (unsigned int l = 1; l <= p; ++l) {
          int ind = k-p+l;
          auto alpha = ubar[k+l]-x[j];
          if (abs(alpha) == 0.0) 
            for (int nr = 0; nr < num_row; ++nr) for(int nc = 0; nc < num_col; ++nc) nurbs.cpi()(nr,nc,ind-1) = nurbs.cpi()(nr,nc,ind);
          else {
            alpha = alpha / (ubar[k+l]-u[i-p+l]);
            for (int nr = 0; nr < num_row; ++nr) for(int nc = 0; nc < num_col; ++nc)  {{
              //PointContainerTypes<T>::point_type point = PointContainerTypes<T>::point_type(alpha*cp[nurbs.cpi()[nr][nc][ind-1]] + ((T)1 - alpha)*cp[nurbs.cpi()[nr][nc][ind]]);
              PointContainerTypes<T>::point_type point_1 = model.get_point(nurbs.cpi()(nr,nc,ind-1));
              PointContainerTypes<T>::point_type point_2 = model.get_point(nurbs.cpi()(nr,nc,ind));
              PointContainerTypes<T>::point_type point = alpha*point_1 + ((T)1-alpha)*point_2;

              set_point(nr, nc, ind-1, point);
            }}
          }
        }
        ubar[k] = x[j];
        --k;
      }    
    }

    u = ubar;
  }


  template <typename T> void h_refine(squanch::Model<T>& model, Nurbs<T, 2>& nurbs, std::vector<T> x1, std::vector<T> x2)
  {
    h_refine(model, nurbs, 0, x1);
    h_refine(model, nurbs, 1, x2);

    knit_closed_shape(nurbs);
  }

  template <typename T> void h_refine(squanch::Model<T>& model, Nurbs<T, 3>& nurbs, std::vector<T> x1, std::vector<T> x2, std::vector<T> x3)
  {
    h_refine(model, nurbs, 0, x1);
    h_refine(model, nurbs, 1, x2);
    h_refine(model, nurbs, 2, x3);
    knit_closed_shape(nurbs);
  }

  template <typename T> void h_refine_elements(squanch::Model<T>& model, Nurbs<T, 1>& nurbs, unsigned int f1)
  {
    h_refine(model, nurbs, detail::find_new_knots(nurbs.curves()[0].knots(), f1));
  }

  template <typename T> void h_refine_elements(squanch::Model<T>& model, Nurbs<T, 2>& nurbs, unsigned int f1, unsigned int f2)
  {
    h_refine(model, nurbs, detail::find_new_knots(nurbs.curves()[0].knots(), f1), detail::find_new_knots(nurbs.curves()[1].knots(), f2));
    knit_closed_shape(nurbs);
  }

  template <typename T> void h_refine_elements(squanch::Model<T>& model, Nurbs<T, 3>& nurbs, unsigned int f1, unsigned int f2, unsigned int f3)
  {
    h_refine(model, nurbs, detail::find_new_knots(nurbs.curves()[0].knots(), f1), detail::find_new_knots(nurbs.curves()[1].knots(), f2), detail::find_new_knots(nurbs.curves()[2].knots(), f3));
    knit_closed_shape(nurbs);
  }


  template <typename T> void p_refine(squanch::Model<T>& model, Nurbs<T, 1>& nurbs, int t)
  {
		if (t == 0) return;

    auto&& curve = nurbs.curves()[0];
    auto&& cpi = nurbs.cpi();

    int i,j,k ;
    int n = curve.n()-1;
    int p = curve.p();
    int m = n+p+1;
    int ph = p+t ;
    int ph2 = ph/2 ;

    //auto Uh = curve.knots();

    //auto  id_s_s = id_s;

    PointContainerTypes<T>::map_type Qw;

    Eigen::Matrix<T, -1, -1> bezalfs(p+t+1,p+1) ; // coefficients for degree elevating the Bezier segment
    PointContainerTypes<T>::vector_type bpts(p+1); // pth-degree Bezier control points of the current segment
    PointContainerTypes<T>::vector_type ebpts(p+t+1) ; // (p+t)th-degree Bezier control points of the  current segment
    PointContainerTypes<T>::vector_type Nextbpts(p-1) ; // leftmost control points of the next Bezier segment
    Eigen::Matrix<T,-1,1> alphas(p-1) ; // knot insertion alphas.

    // Compute the binomial coefficients
    Eigen::Matrix<T,-1,-1> Bin(ph+1,ph2+1) ;
    detail::binomialCoef(Bin) ;

    // Compute Bezier degree elevation coefficients
    T inv,mpi ;
    bezalfs(0,0) = bezalfs(ph,p) = 1.0 ;
    for(i=1;i<=ph2;i++){
      inv= 1.0/Bin(ph,i) ;
      mpi = std::min(p,i) ;
      for(j=std::max<int>(0,i-t); j<=mpi; j++){
        bezalfs(i,j) = inv*Bin(p,j)*Bin(t,i-j) ;
      }
    }

    for(i=ph2+1;i<ph ; i++){
      mpi = std::min(p,i) ;
      for(j=std::max<int>(0,i-t); j<=mpi ; j++)
        bezalfs(i,j) = bezalfs(ph-i,p-j) ;
    }

    auto size_arr = nurbs.cpi().shape()[0];
    size_arr = size_arr + size_arr*t;
    nurbs.cpi().resizeAndPreserve(size_arr);
    std::vector<T> Uh(size_arr + ph + 1);
    //Uh.resize(size_arr + ph + 1);
    //resize(c.P.n()+c.P.n()*t,ph) ; // Allocate more control points than necessary

    int mh = ph ;
    int kind = ph+1 ;
    T ua = curve.knots()[0];
    Qw[0] = model.get_point(cpi(0));
    T ub = 0.0 ;
    int r=-1 ; 
    int oldr ;
    int a = p ;
    int b = p+1 ; 
    int cind = 1 ;
    int rbz,lbz = 1 ; 
    int mul,save,s;
    T alf;
    int first, last, kj ;
    T den,bet,gam,numer ;

    //P[0] = c.P[0] ;
    for (size_t i = 0; i <= ph; ++i) Uh[i] = ua;


    // Initialize the first Bezier segment

    for(i=0;i<=p ;i++) 
      bpts[i] = model.get_point(nurbs.cpi()(i));

    while(b<m){ // Big loop over knot vector
      i=b ;
      while(b<m && curve.knots()[b] >= curve.knots()[b+1]) // for some odd reasons... == doesn't work
        b++ ;
      mul = b-i+1 ; 
      mh += mul+t ;
      ub = curve.knots()[b];
      oldr = r ;
      r = p-mul ;
      if(oldr>0)
        lbz = (oldr+2)/2 ;
      else
        lbz = 1 ;
      if(r>0) 
        rbz = ph-(r+1)/2 ;
      else
        rbz = ph ;
      if(r>0){ // Insert knot to get Bezier segment
        numer = ub-ua ;
        for(k=p;k>mul;k--){
          alphas[k-mul-1] = numer/(curve.knots()[a+k]-ua) ;
        }
        for(j=1;j<=r;j++){
          save = r-j ; s = mul+j ;
          for(k=p;k>=s;k--){
            bpts[k] = alphas[k-s] * bpts[k]+(1.0-alphas[k-s])*bpts[k-1] ;
          }
          Nextbpts[save] = bpts[p] ;
        }
      }

      for(i=lbz;i<=ph;i++){ // Degree elevate Bezier,  only the points lbz,...,ph are used
        ebpts[i] = PointContainerTypes<T>::point_type::Zero();
        mpi = std::min(p,i) ;
        for(j=std::max<int>(0,i-t); j<=mpi ; j++)
          ebpts[i] += bezalfs(i,j)*bpts[j] ;
      }

      if(oldr>1){ // Must remove knot u=c.U[a] oldr times
        // if(oldr>2) // Alphas on the right do not change
        //	alfj = (ua-U[kind-1])/(ub-U[kind-1]) ;
        first = kind-2 ; last = kind ;
        den = ub-ua ;
        bet = (ub-Uh[kind-1])/den ;  
        for(int tr=1; tr<oldr; tr++){ // Knot removal loop
          i = first ; j = last ;
          kj = j-kind+1 ;
          while(j-i>tr){ // Loop and compute the new control points for one removal step
            if(i<cind){
              alf=(ub-Uh[i])/(ua-Uh[i]) ;
              PointContainerTypes<T>::point_type po = PointContainerTypes<T>::point_type(alf*Qw[i] + (1.0-alf)*Qw[i-1]);
              Qw[i] = po;
            }
            if( j>= lbz){
              if(j-tr <= kind-ph+oldr){
                gam = (ub-Uh[j-tr])/den ;
                ebpts[kj] = gam*ebpts[kj] + (1.0-gam)*ebpts[kj+1] ;
              }
              else{
                ebpts[kj] = bet*ebpts[kj]+(1.0-bet)*ebpts[kj+1] ;
              }
            }
            ++i ; --j; --kj ;
          }
          --first ; ++last ;
        }
      }

      if(a!=p) // load the knot u=c.U[a]
        for(i=0;i<ph-oldr; i++){
          Uh[kind++] = ua ; 
        }
        for(j=lbz; j<=rbz ; j++) { // load control points onto the curve
          Qw[cind] = ebpts[j];
          cind++;
        }

        if(b<m){ // Set up for next pass through loop
          for(j=0;j<r;j++)
            bpts[j] = Nextbpts[j] ;
          for (j = r; j <= p; j++) 
            bpts[j] = model.get_point(nurbs.cpi()(b - p + j));
          a=b ; 
          b++ ;
          ua = ub ;
        }
        else{
          for(i=0;i<=ph;i++)
            Uh[kind+i] = ub ;
        }
    }

    //resize(mh-ph,ph) ; // Resize to the proper number of control points
    nurbs.cpi().resizeAndPreserve(mh - ph);
    Uh.resize(mh + 1);
    
    curve.set_p(curve.p() + t);
    curve.knots() = Uh;
    //cpi.resize(boost::extents[curve.n()]);
    
    for (auto&& pair : Qw) {
      int id;
      std::tie(id, std::ignore) = model.new_point(pair.second.x(), pair.second.y(), pair.second.z(), pair.second.w());
      cpi(pair.first) = id;
    }
  }



  template <typename T> void p_refine(squanch::Model<T>& model, Nurbs<T, 2>& nurbs, unsigned int curve_i, int t)
  {
		if (t == 0) return;

    auto&& curve = nurbs.curves()[curve_i];
    auto&& cpi = nurbs.cpi();
    auto cpisize = cpi.shape();

    int i,j,k ;
    int n = curve.n()-1;
    int p = curve.p();
    int m = n+p+1;
    int ph = p+t ;
    int ph2 = ph/2 ;

    //auto  id_s_s = id_s;

    PointContainerTypes<T>::map_type Qw;
    blitz::Array<int, 2> qcpi;

    size_t dim2 = 0;
    if (curve_i == 0) dim2 = cpisize[1]; else dim2 = cpisize[0];

    Eigen::Matrix<T, -1, -1> bezalfs(p+t+1,p+1) ; // coefficients for degree elevating the Bezier segment
    boost::multi_array<PointContainerTypes<T>::vector_type, 1> bpts(boost::extents[dim2]); for (auto mi = 0; mi < dim2; ++mi) bpts[mi].resize(p+1); // pth-degree Bezier control points of the current segment
    boost::multi_array<PointContainerTypes<T>::vector_type, 1> ebpts(boost::extents[dim2]); for (auto mi = 0; mi < dim2; ++mi) ebpts[mi].resize(p+t+1) ; // (p+t)th-degree Bezier control points of the  current segment
    boost::multi_array<PointContainerTypes<T>::vector_type, 1> Nextbpts(boost::extents[dim2]); for (auto mi = 0; mi < dim2; ++mi) Nextbpts[mi].resize(p-1) ; // leftmost control points of the next Bezier segment
    Eigen::Matrix<T,-1,1> alphas(p-1); // knot insertion alphas.

    // Compute the binomial coefficients
    Eigen::Matrix<T,-1,-1> Bin(ph+1,ph2+1) ;
    detail::binomialCoef(Bin) ;

    // Compute Bezier degree elevation coefficients
    T inv,mpi ;
    bezalfs(0,0) = bezalfs(ph,p) = 1.0 ;
    for(i=1;i<=ph2;i++){
      inv= 1.0/Bin(ph,i) ;
      mpi = std::min(p,i) ;
      for(j=std::max<int>(0,i-t); j<=mpi; j++){
        bezalfs(i,j) = inv*Bin(p,j)*Bin(t,i-j) ;
      }
    }
    for(i=ph2+1;i<ph ; i++){
      mpi = std::min(p,i) ;
      for(j=std::max<int>(0,i-t); j<=mpi ; j++)
        bezalfs(i,j) = bezalfs(ph-i,p-j) ;
    }

    // resize cpi
    auto arrp = nurbs.cpi().shape();
    std::array<size_t, 2> arr; arr[0] = arrp[0]; arr[1] = arrp[1];
    arr[curve_i] += arr[curve_i]*t;
    qcpi.resizeAndPreserve(dim2, arr[curve_i]); // in order to provide similar notation between u and v directions. 
    std::fill(qcpi.data(), qcpi.data() + qcpi.size(), -1); 

    auto set_point = [&Qw, &qcpi, &model](size_t p1, size_t p2, PointContainerTypes<T>::point_type& ptn)
    {
      auto&& tup = model.new_point();
      auto id = std::get<0>(tup);
      qcpi((int)p1,(int)p2) = id;
      Qw[id] = ptn;
    };
		auto set_point_w_id = [&Qw, &qcpi, &model](size_t p1, size_t p2, int id)
		{
			auto&& tup = model.new_point();
			auto&& ptn = model.get_point(id);
			qcpi((int)p1, (int)p2) = tup.first;
			Qw[tup.first] = ptn;
		};

    std::vector<T> Uh(arr[curve_i] + ph + 1);

    int mh = ph ;
    int kind = ph+1 ;
    T ua = curve.knots()[0];
    T ub = 0.0 ;
    int r=-1 ; 
    int oldr ;
    int a = p ;
    int b = p+1 ; 
    int cind = 1 ;
    int rbz,lbz = 1 ; 
    int mul,save,s;
    T alf;
    int first, last, kj ;
    T den,bet,gam,numer ;
    

    for (size_t i = 0; i <= ph; ++i) Uh[i] = ua;
    // Initialize the first Bezier segment

    if (curve_i == 0)
      for (auto di = 0; di < dim2; ++di) set_point_w_id(di, 0, (cpi(0,di)));
    else
			for (auto di = 0; di < dim2; ++di) set_point_w_id(di, 0, (cpi(di, 0)));
    
    if (curve_i == 0)
      for (auto di = 0; di < dim2; ++di) for(i=0;i<=p;i++) bpts[di][i] = model.get_point(nurbs.cpi()(i,di));
    else
      for (auto di = 0; di < dim2; ++di) for(i=0;i<=p;i++) bpts[di][i] = model.get_point(nurbs.cpi()(di,i));

    while(b<m){ // Big loop through knot vector
      i=b ;
      while(b<m && curve.knots()[b] >= curve.knots()[b+1]) // for some odd reasons... == doesn't work
        b++ ;
      mul = b-i+1 ; 
      mh += mul+t ;
      ub = curve.knots()[b];
      oldr = r ;
      r = p-mul ;
      if(oldr>0)
        lbz = (oldr+2)/2 ;
      else
        lbz = 1 ;
      if(r>0) 
        rbz = ph-(r+1)/2 ;
      else
        rbz = ph ;
      if(r>0){ // Insert knot to get Bezier segment
        numer = ub-ua ;
        for(k=p;k>mul;k--){
          alphas[k-mul-1] = numer/(curve.knots()[a+k]-ua) ;
        }
        for(j=1;j<=r;j++){
          save = r-j ; s = mul+j ;
          for(k=p;k>=s;k--){
            for (auto di = 0; di < dim2; ++di){
              bpts[di][k] = alphas[k-s] * bpts[di][k]+(1.0-alphas[k-s])*bpts[di][k-1] ;
            }
          }
          for (auto di = 0; di < dim2; ++di) Nextbpts[di][save] = bpts[di][p] ;
        }
      }

      for(i=lbz;i<=ph;i++) { // Degree elevate Bezier,  only the points lbz,...,ph are used
        for (auto di = 0; di < dim2; ++di) ebpts[di][i] = PointContainerTypes<T>::point_type::Zero();
        mpi = std::min(p,i) ;
        for(j=std::max<int>(0,i-t); j<=mpi ; j++)
          for (auto di = 0; di < dim2; ++di) ebpts[di][i] += bezalfs(i,j)*bpts[di][j] ;
      }

      if(oldr>1){ // Must remove knot u=c.U[a] oldr times
        // if(oldr>2) // Alphas on the right do not change
        //	alfj = (ua-U[kind-1])/(ub-U[kind-1]) ;
        first = kind-2 ; last = kind ;
        den = ub-ua ;
        bet = (ub-Uh[kind-1])/den ;  
        for(int tr=1; tr<oldr; tr++){ // Knot removal loop
          i = first ; j = last ;
          kj = j-kind+1 ;
          while(j-i>tr){ // Loop and compute the new control points for one removal step
            if(i<cind){
              alf=(ub-Uh[i])/(ua-Uh[i]) ;
              PointContainerTypes<T>::point_type po = PointContainerTypes<T>::point_type(alf*Qw[i] + (1.0-alf)*Qw[i-1]);
              for (auto di = 0; di < dim2; ++di) set_point(di, i, po);
            }
            if( j>= lbz){
              if(j-tr <= kind-ph+oldr){
                gam = (ub-Uh[j-tr])/den ;
                for (auto di = 0; di < dim2; ++di) ebpts[di][kj] = gam*ebpts[di][kj] + (1.0-gam)*ebpts[di][kj+1] ;
              }
              else{
                for (auto di = 0; di < dim2; ++di) ebpts[di][kj] = bet*ebpts[di][kj]+(1.0-bet)*ebpts[di][kj+1] ;
              }
            }
            ++i ; --j; --kj ;
          }
          --first ; ++last ;
        }
      }

      if(a!=p) // load the knot u=c.U[a]
        for(i=0;i<ph-oldr; i++){
          Uh[kind++] = ua ; 
        }
        for(j=lbz; j<=rbz ; j++) { // load control points onto the curve
          for (auto di = 0; di < dim2; ++di) set_point(di, cind, ebpts[di][j]);
          cind++;
        }

        if(b<m){ // Set up for next pass through loop
          for (auto di = 0; di < dim2; ++di) for (j=0;j<r;j++)  bpts[di][j] = Nextbpts[di][j];
          if (curve_i == 0)
            for (auto di = 0; di < dim2; ++di) for(j=r;j<=p;j++) bpts[di][j] = model.get_point(nurbs.cpi()((int)(b-p+j),di));
          else
            for (auto di = 0; di < dim2; ++di) for(j=r;j<=p;j++) bpts[di][j] = model.get_point(nurbs.cpi()(di,(int)(b-p+j)));

          a=b ; 
          b++ ;
          ua = ub ;
        }
        else{
          for(i=0;i<=ph;i++) Uh[kind+i] = ub;
        }
    }

    // Resize to the proper number of control points
    Uh.resize(mh + 1);
    curve.set_p(curve.p() + t);
    curve.knots() = Uh;
    if (curve_i == 0)
      cpi.resizeAndPreserve(mh - ph, dim2);
    else
      cpi.resizeAndPreserve(dim2, mh - ph);

    if (curve_i == 0) {
      for (auto di = 0; di < dim2; ++di) {
        for (int i = 0; i < arr[curve_i]; ++i) {
          auto id = qcpi(di,i);
          if (id != -1) {
            cpi((int)i,di) = id;
            model.get_point(id) = Qw[id];
          }
        }
      }
    } else {
      for (auto di = 0; di < dim2; ++di) {
        for (int i = 0; i < arr[curve_i]; ++i) {
          auto id = qcpi(di,i);
          if (id != -1) {
            cpi(di,(int)i) = id;
            model.get_point(id) = Qw[id];
          }
        }
      }
    }
  }


  template <typename T> void p_refine(squanch::Model<T>& model, Nurbs<T, 3>& nurbs, unsigned int curve_i, int t)
  {
		if (t == 0) return;

    auto&& curve = nurbs.curves()[curve_i];
    auto&& cpi = nurbs.cpi();

    int i,j,k ;
    int n = curve.n()-1;
    int p = curve.p();
    int m = n+p+1;
    int ph = p+t ;
    int ph2 = ph/2 ;

    //auto  id_s_s = id_s;

    PointContainerTypes<T>::map_type Qw;
    blitz::Array<int, 3> qcpi;

    auto arr = nurbs.cpi().shape();
    size_t dim2 = 0;
    size_t dim3 = 0;
    switch (curve_i)
    {
    case 0:
      dim2 = arr[1];
      dim3 = arr[2];
      break;
    case 1:
      dim2 = arr[0];
      dim3 = arr[2];
      break;
    case 2:
      dim2 = arr[0];
      dim3 = arr[1];
      break;
    }

    Eigen::Matrix<T, -1, -1> bezalfs(p+t+1,p+1) ; // coefficients for degree elevating the Bezier segment
    boost::multi_array<PointContainerTypes<T>::vector_type, 2> bpts    (boost::extents[dim2][dim3]); 
    for (auto mi = 0; mi < dim2; ++mi) for (auto mj = 0; mj < dim3; ++mj) bpts[mi][mj].resize(p+1); // pth-degree Bezier control points of the current segment
    boost::multi_array<PointContainerTypes<T>::vector_type, 2> ebpts   (boost::extents[dim2][dim3]); 
    for (auto mi = 0; mi < dim2; ++mi)  for (auto mj = 0; mj < dim3; ++mj) ebpts[mi][mj].resize(p+t+1) ; // (p+t)th-degree Bezier control points of the  current segment
    boost::multi_array<PointContainerTypes<T>::vector_type, 2> Nextbpts(boost::extents[dim2][dim3]); 
    for (auto mi = 0; mi < dim2; ++mi)  for (auto mj = 0; mj < dim3; ++mj) Nextbpts[mi][mj].resize(p-1) ; // leftmost control points of the next Bezier segment
    Eigen::Matrix<T,-1,1> alphas(p-1); // knot insertion alphas.

    // Compute the binomial coefficients
    Eigen::Matrix<T,-1,-1> Bin(ph+1,ph2+1) ;
    detail::binomialCoef(Bin) ;

    // Compute Bezier degree elevation coefficients
    T inv,mpi ;
    bezalfs(0,0) = bezalfs(ph,p) = 1.0 ;
    for(i=1;i<=ph2;i++){
      inv= 1.0/Bin(ph,i) ;
      mpi = std::min(p,i) ;
      for(j=std::max<int>(0,i-t); j<=mpi; j++){
        bezalfs(i,j) = inv*Bin(p,j)*Bin(t,i-j) ;
      }
    }
    for(i=ph2+1;i<ph ; i++){
      mpi = std::min(p,i) ;
      for(j=std::max<int>(0,i-t); j<=mpi ; j++)
        bezalfs(i,j) = bezalfs(ph-i,p-j) ;
    }

    // resize cpi
    arr[curve_i] += arr[curve_i]*t;
    qcpi.resize(dim2,dim3,arr[curve_i]); // in order to provide similar notation between u and v directions. 
    std::fill(qcpi.data(), qcpi.data() + qcpi.size(), -1); 

    auto set_point = [&Qw, &qcpi, &model](size_t p1, size_t p2, size_t p3, PointContainerTypes<T>::point_type& ptn)
    {
      auto&& tup = model.new_point();
      auto id = std::get<0>(tup);
      qcpi((int)p1,(int)p2,(int)p3) = id;
      Qw[id] = ptn;
    };

		auto set_point_w_id = [&Qw, &qcpi, &model](size_t p1, size_t p2, size_t p3, int ptn_id)
		{
			auto&& tup = model.new_point();
			auto id = std::get<0>(tup);
			auto&& ptn = model.get_point(ptn_id);
			qcpi((int)p1, (int)p2, (int)p3) = id;
			Qw[id] = ptn;
		};

    std::vector<T> Uh(arr[curve_i] + ph + 1);

    int mh = ph ;
    int kind = ph+1 ;
    T ua = curve.knots()[0];
    T ub = 0.0 ;
    int r=-1 ; 
    int oldr ;
    int a = p ;
    int b = p+1 ; 
    int cind = 1 ;
    int rbz,lbz = 1 ; 
    int mul,save,s;
    T alf;
    int first, last, kj ;
    T den,bet,gam,numer ;
    

    for (size_t i = 0; i <= ph; ++i) Uh[i] = ua;
    // Initialize the first Bezier segment

    switch(curve_i)
    {
    case 0:
      for (auto di = 0; di < dim2; ++di) for (auto di3 = 0; di3 < dim3; ++di3) set_point_w_id(di, di3, 0, (cpi(0,di,di3)));
      for (auto di = 0; di < dim2; ++di) for (auto di3 = 0; di3 < dim3; ++di3) for(i=0;i<=p;i++) bpts[di][di3][i] = model.get_point(nurbs.cpi()(i,di,di3));
      break;
    case 1:
      for (auto di = 0; di < dim2; ++di) for (auto di3 = 0; di3 < dim3; ++di3) set_point_w_id(di, di3, 0, (cpi(di, 0, di3)));
      for (auto di = 0; di < dim2; ++di) for (auto di3 = 0; di3 < dim3; ++di3) for (i = 0; i <= p; i++) bpts[di][di3][i] = model.get_point(nurbs.cpi()(di, i, di3));
      break;
    case 2:
      for (auto di = 0; di < dim2; ++di) for (auto di3 = 0; di3 < dim3; ++di3) set_point_w_id(di, di3, 0, (cpi(di, di3, 0)));
      for (auto di = 0; di < dim2; ++di) for (auto di3 = 0; di3 < dim3; ++di3) for (i = 0; i <= p; i++) bpts[di][di3][i] = model.get_point(nurbs.cpi()(di, di3, i));
      break;
    }

    while(b<m){ // Big loop through knot vector
      i=b ;
      while(b<m && curve.knots()[b] >= curve.knots()[b+1]) // for some odd reasons... == doesn't work
        b++ ;
      mul = b-i+1 ; 
      mh += mul+t ;
      ub = curve.knots()[b];
      oldr = r ;
      r = p-mul ;
      if(oldr>0)
        lbz = (oldr+2)/2 ;
      else
        lbz = 1 ;
      if(r>0) 
        rbz = ph-(r+1)/2 ;
      else
        rbz = ph ;
      if(r>0){ // Insert knot to get Bezier segment
        numer = ub-ua ;
        for(k=p;k>mul;k--){
          alphas[k-mul-1] = numer/(curve.knots()[a+k]-ua) ;
        }
        for(j=1;j<=r;j++){
          save = r-j ; s = mul+j ;
          for(k=p;k>=s;k--){
            for (auto di = 0; di < dim2; ++di)  for (auto di3 = 0; di3 < dim3; ++di3) {
              bpts[di][di3][k] = alphas[k-s] * bpts[di][di3][k]+(1.0-alphas[k-s])*bpts[di][di3][k-1];
            }
          }
          for (auto di = 0; di < dim2; ++di)  for (auto di3 = 0; di3 < dim3; ++di3) Nextbpts[di][di3][save] = bpts[di][di3][p];
        }
      }

      for(i=lbz;i<=ph;i++){ // Degree elevate Bezier,  only the points lbz,...,ph are used
        for (auto di = 0; di < dim2; ++di)  for (auto di3 = 0; di3 < dim3; ++di3) ebpts[di][di3][i] = PointContainerTypes<T>::point_type::Zero();
        mpi = std::min(p,i) ;
        for(j=std::max<int>(0,i-t); j<=mpi ; j++)
          for (auto di = 0; di < dim2; ++di)  for (auto di3 = 0; di3 < dim3; ++di3) ebpts[di][di3][i] += bezalfs(i,j)*bpts[di][di3][j] ;
      }

      if(oldr>1){ // Must remove knot u=c.U[a] oldr times
        // if(oldr>2) // Alphas on the right do not change
        //	alfj = (ua-U[kind-1])/(ub-U[kind-1]) ;
        first = kind-2 ; last = kind ;
        den = ub-ua ;
        bet = (ub-Uh[kind-1])/den ;  
        for(int tr=1; tr<oldr; tr++){ // Knot removal loop
          i = first ; j = last ;
          kj = j-kind+1 ;
          while(j-i>tr){ // Loop and compute the new control points for one removal step
            if(i<cind){
              alf=(ub-Uh[i])/(ua-Uh[i]) ;
              PointContainerTypes<T>::point_type po = PointContainerTypes<T>::point_type(alf*Qw[i] + (1.0-alf)*Qw[i-1]);
              for (auto di = 0; di < dim2; ++di)  for (auto di3 = 0; di3 < dim3; ++di3) set_point(di, di3, i, po);
            }
            if( j>= lbz){
              if(j-tr <= kind-ph+oldr){
                gam = (ub-Uh[j-tr])/den ;
                for (auto di = 0; di < dim2; ++di) for(auto di3 = 0; di3 < dim3; ++di3) ebpts[di][di3][kj] = gam*ebpts[di][di3][kj] + (1.0-gam)*ebpts[di][di3][kj+1];
              }
              else{
                for (auto di = 0; di < dim2; ++di) for(auto di3 = 0; di3 < dim3; ++di3) ebpts[di][di3][kj] = bet*ebpts[di][di3][kj]+(1.0-bet)*ebpts[di][di3][kj+1];
              }
            }
            ++i ; --j; --kj ;
          }
          --first ; ++last ;
        }
      }

      if(a!=p) // load the knot u=c.U[a]
        for(i=0;i<ph-oldr; i++){
          Uh[kind++] = ua ; 
        }
        for(j=lbz; j<=rbz ; j++) { // load control points onto the curve
          for (auto di = 0; di < dim2; ++di)  for (auto di3 = 0; di3 < dim3; ++di3) set_point(di, di3, cind, ebpts[di][di3][j]);
          cind++;
        }

        if(b<m){ // Set up for next pass through loop
          for (auto di = 0; di < dim2; ++di)  for (auto di3 = 0; di3 < dim3; ++di3) for (j=0;j<r;j++)  bpts[di][di3][j] = Nextbpts[di][di3][j];
          switch(curve_i)
          {
          case 0:
            for (auto di = 0; di < dim2; ++di) for (auto di3 = 0; di3 < dim3; ++di3) for (j = r; j <= p; j++) bpts[di][di3][j] = model.get_point(nurbs.cpi()(b - p + j, di, di3));
            break;
          case 1:
            for (auto di = 0; di < dim2; ++di) for (auto di3 = 0; di3 < dim3; ++di3) for (j = r; j <= p; j++) bpts[di][di3][j] = model.get_point(nurbs.cpi()(di, b - p + j, di3));
            break;
          case 2:
            for (auto di = 0; di < dim2; ++di) for (auto di3 = 0; di3 < dim3; ++di3) for (j = r; j <= p; j++) bpts[di][di3][j] = model.get_point(nurbs.cpi()(di, di3, b - p + j));
            break;
          }
          a=b ; 
          b++ ;
          ua = ub ;
        }
        else{
          for(i=0;i<=ph;i++) Uh[kind+i] = ub;
        }
    }

    Uh.resize(mh + 1);
    
    curve.set_p(curve.p() + t);
    curve.knots() = Uh;

    switch(curve_i)
    {
    case 0:
      cpi.resize(mh-ph,dim2,dim3);
      for (auto di = 0; di < dim2; ++di) for(auto di3 = 0; di3 < dim3; ++di3) {
        for (int mi = 0; mi < arr[curve_i]; ++mi) {
          auto id = qcpi(di,di3,mi);
          if (id != -1) {
            cpi(mi,di,di3) = id;
            model.get_point(id) = Qw[id];
          }
        }
      }
      break;
    case 1:
      cpi.resize(dim2,mh-ph,dim3);
      for (auto di = 0; di < dim2; ++di) for(auto di3 = 0; di3 < dim3; ++di3) {
        for (int mi = 0; mi < arr[curve_i]; ++mi) {
          auto id = qcpi(di,di3,mi);
          if (id != -1) {
            cpi(di,mi,di3) = id;
            model.get_point(id) = Qw[id];
          }
        }
      }
      break;
    case 2:
      cpi.resize(dim2,dim3,mh-ph);
      for (auto di = 0; di < dim2; ++di) for(auto di3 = 0; di3 < dim3; ++di3) {
        for (int mi = 0; mi < arr[curve_i]; ++mi) {
          auto id = qcpi(di,di3,mi);
          if (id != -1) {
            cpi(di,di3,mi) = id;
            model.get_point(id) = Qw[id];
          }
        }
      }
      break;
    }
  }

  template <typename T> void p_refine(squanch::Model<T>& model, Nurbs<T, 2>& nurbs, const std::array<unsigned int, 2>& t)
  {
    p_refine(model, nurbs, 0, t[0]);
    p_refine(model, nurbs,  1, t[1]);
    knit_closed_shape(nurbs);
  }

  template <typename T> void p_refine(squanch::Model<T>& model, Nurbs<T, 3>& nurbs, const std::array<unsigned int, 3>& t)
  {
    p_refine(model, nurbs, 0, t[0]);
    p_refine(model, nurbs, 1, t[1]);
    p_refine(model, nurbs, 2, t[2]);
    knit_closed_shape(nurbs);
  }
}