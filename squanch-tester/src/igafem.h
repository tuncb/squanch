#pragma once
#include <vector>
#include <Eigen/core>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>

namespace igafem {

double **init2DArray(int x, int y);
void free2Darray(double **array, int x);

size_t findspan(size_t n, size_t p, double u, const std::vector<double> &U);
void BasisFuns(size_t i, double u, unsigned int p, const std::vector<double> &U,
               Eigen::VectorXd &N);
void dersBasisFuns(int i, double u, int p, int order, double knot[],
                   double** ders);

void nurbs1DBasisDers(double u1, unsigned int p, const std::vector<double> &knotU,
                      const std::vector<double> &weight, Eigen::VectorXd &R,
                      Eigen::VectorXd &dRdxi);

void nurbs2DBasisDers(double u1, double u2, unsigned int p, unsigned int q,
                      const std::vector<double> &knotU,
                      const std::vector<double> &knotV,
                      const std::vector<double> &weight, Eigen::VectorXd &R,
                      Eigen::VectorXd &dRdxi, Eigen::VectorXd &dRdet);

void nurbs3DBasisDers(double u1, double u2, double u3, 
                      unsigned int p, unsigned int q, unsigned int r,
                      const std::vector<double> &knotU,
                      const std::vector<double> &knotV,
                      const std::vector<double> &knotW,
                      const std::vector<double> &weight, Eigen::VectorXd &R,
                      Eigen::VectorXd &dRdxi, Eigen::VectorXd &dRdet, Eigen::VectorXd &dRdze);
}