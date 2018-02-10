#pragma once
#include <squanch/Model.h>
#include <squanch/PointContainerTypes.h>
#include <squanch/nurbs.h>
#include <Eigen/core>
#include <Eigen/geometry>

namespace squanch {

  namespace detail {

    template <class T>
    int intersectLine(const typename PointContainerTypes<T>::point_type& p1, const typename PointContainerTypes<T>::point_type& t1, 
                      const typename PointContainerTypes<T>::point_type& p2, const typename PointContainerTypes<T>::point_type& t2, 
                      typename PointContainerTypes<T>::point_type& p){
      // a line is written like 
      // L(t) = Q + u*t
      // u is the parametric value P is a point in the line and T is the tangent
      // a plane is of the form
      // (X-P).v = 0 
      // where P is a point where the plane goes through and v is the normal to it
      // and X is (x,y,z)
      // solving for u

      T u = p1.dot(p2);
      if (abs((T)1 - u) < 1e-7) return 0; 

      PointContainerTypes<T>::point_type v,px ;

      //px = crossProduct(t1,p1-p2) ;
      //v = crossProduct(px,t1) ;
      px = t1.cross3(t2) ; 
      v = px.cross3(t1) ; 

      T t = (p1-p2).dot(v) ;
      T vw = v.dot(t2) ;
      if(vw * vw < 1e-7) return 0 ;
      //vw /= t;
      t /= vw;
      p = p2+(((p1-p2).dot(v))/vw)*t2 ;
      return 1 ;
    }

  }

  template <typename T, typename PointType>
  class Circle
  {
  public:
    T start_angle;
    T end_angle;
    PointType origin;
    PointType N1;
    PointType N2;
    T r;

    Nurbs<T, 1>& to_nurbs(squanch::Model<T>& model)
    {
      const T pi = (T)3.14159265359;

      auto&& nurbs = model.new_nurbs<1>();
      auto&& curve = nurbs.curves()[0];

      // check if the circle is closed or not 
      if (abs(2*pi - (end_angle - start_angle)) < 0.00001f) nurbs.curves()[0].set_closed(true);
      else nurbs.curves()[0].set_closed(false); 

      auto set_point = [&nurbs, &model](int pos, typename PointContainerTypes<T>::point_type& p, T w){
        if (nurbs.curves()[0].is_closed() && pos == 8) {
          nurbs.cpi()(8) = nurbs.cpi()(0);
          return;
        }

        auto&& pr = model.new_point(p.x(), p.y(), p.z(), (T)1);
				nurbs.cpi()(pos) = pr.first;
        pr.second *= w;
      };

      PointContainerTypes<T>::point_type O(this->origin.x(), this->origin.y(), this->origin.z(), (T)1);
      PointContainerTypes<T>::point_type X(this->N1.x(), this->N1.y(), this->N1.z(), (T)1);
      PointContainerTypes<T>::point_type Y(this->N2.x(), this->N2.y(), this->N2.z(), (T)1);

      T theta,angle,dtheta ;
      int narcs;

      auto as = this->start_angle;
      auto ae = this->end_angle;


      while(ae<as) ae += 2*pi;

      theta = ae-as ;
      if(theta<=pi/2.0)	
        narcs = 1 ;
      else{
        if(theta<=pi)
          narcs = 2 ;
        else{
          if(theta<=1.5*pi)
            narcs = 3 ;
          else
            narcs = 4 ;
        }
      }
      dtheta = theta/(T)narcs ;
      int n = 2*narcs+1 ; // n control points ;
      T w1 = cos(dtheta/2.0) ; // dtheta/2.0 is base angle

      PointContainerTypes<T>::point_type T0,T2, P0, P1, P2;
      //auto&& P0_t = model.new_point();
      //auto&& P1_t = model.new_point();
      //auto&& P2_t = model.new_point();

      //auto&& P0 = std::get<1>(P0_t);
      //auto&& P1 = std::get<1>(P1_t);
      //auto&& P2 = std::get<1>(P2_t);

      //std::cout << P0.transpose() << std::endl;
      //std::cout << O.transpose() << std::endl;
      //std::cout << X.transpose() << std::endl;
      //std::cout << Y.transpose() << std::endl;
      //std::cout << as << std::endl;
      P0 = O + r*cos(as)*X + r*sin(as)*Y ; 
/*      std::cout << P0.transpose() << std::endl*/;
      T0 = -sin(as)*X + cos(as)*Y ; // initialize start values

      curve.set_p(2);
      nurbs.cpi().resize(n);
      auto&& U = curve.knots();
      U.resize(n+3);

      set_point(0, P0, (T)1);
      //P[0] = P0;
      int i ;
      int index = 0 ;
      angle = as;
      for(i=1;i<=narcs;++i){
        angle += dtheta ;
        P2 = O + r*cos(angle)*X + r*sin(angle)*Y ;  
        set_point(index+2, P2, (T)1);
        //P[index+2] = P2 ;
        T2 = -sin(angle)*X + cos(angle)*Y ;
        detail::intersectLine<T>(P0,T0,P2,T2,P1);
        set_point(index+1, P1, w1);
        //P[index+1] = P1 ;
        //P[index+1] *= w1 ;
        index += 2 ;
        if(i<narcs){
          P0 = P2 ;
          T0 = T2 ;
        }
      }
      int j = 2*narcs+1 ; // load the knot vector
      for(i=0;i<3;++i){
        U[i] = 0.0 ;
        U[i+j] = 1.0 ;
      }
      switch(narcs){
      case 1: break ;
      case 2: 
        U[3] = U[4] = 0.5 ;
        break ;
      case 3:
        U[3] = U[4] = 1.0/3.0 ;
        U[5] = U[6] = 2.0/3.0 ;
        break ;
      case 4:
        U[3] = U[4] = 0.25 ;
        U[5] = U[6] = 0.50 ;  
        U[7] = U[8] = 0.75 ;
        break ;    
      }

      return nurbs;
    }
  };

}