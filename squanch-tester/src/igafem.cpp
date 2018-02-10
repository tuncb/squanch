#include "igafem.h"

#include <iostream>

double** igafem::init2DArray(int x, int y)
{
    double **array = (double **)malloc(x * sizeof(double *));

    int c;
    for(c = 0; c < x; c++)
    {
        array[c] = (double*)malloc(y * sizeof(double));
    }

    int d;
    return array;
}

void igafem::free2Darray(double **array, int x)
{
    int c;
    for(c = 0; c < x; c++)
    {
        free(array[c]);
    }
    free(array);
}

size_t igafem::findspan(size_t n, size_t p, double u, const std::vector<double>& U)
{
    int low, high, mid;                   
    // special case 
    if (u == U[n+1]) return(n); 
    // do binary search 
    low = p; 
    high = n + 1; 
    mid = (low + high) / 2; 
    while (u < U[mid] || u >= U[mid+1])  {
      if (u < U[mid]) 
        high = mid; 
      else 
        low = mid;                   

      mid = (low + high) / 2; 
    } 
    return(mid); 
}

void igafem::BasisFuns(size_t i, double u, unsigned int p, const std::vector<double>& U, Eigen::VectorXd& N)
{
   /* INPUT: 
  i - knot span  ( from FindSpan() ) 
  u - parametric point 
  p - spline degree 
  U - knot sequence 
            OUTPUT:    
  N - Basis functions vector[p+1] 
            Algorithm A2.2 from 'The NURBS BOOK' pg70. */
                                 
    int j,r; 
    double saved, temp; 
    std::vector<double> left(p+1) ;
    std::vector<double> right(p+1) ;

    N[0] = 1.0; 
    for (j = 1; j <= p; j++)
    {  
        left[j]  = u - U[i+1-j]; 
        right[j] = U[i+j] - u; 
        saved = 0.0; 
        for (r = 0; r < j; r++) 
        { 
            temp = N[r] / (right[r+1] + left[j-r]); 
            N[r] = saved + right[r+1] * temp; 
            saved = left[j-r] * temp; 
        } 
        N[j] = saved; 
    }
}


void igafem::dersBasisFuns(int i, double u, int p, int order, double knot[], double **ders)
{
	/*
	//	Calculate the non-zero derivatives of the b-spline functions
	*/
	
   double saved,temp;
   int j,k,j1,j2,r;
  
   double *left  = (double *)malloc(sizeof(double)*(p+1));
   double *right = (double *)malloc(sizeof(double)*(p+1));
	
   double **ndu  = init2DArray(p+1, p+1);
   double **a    = init2DArray(p+1, p+1);
   
   ndu[0][0]=1.;
   for( j=1; j<=p; j++ )
   {
     left[j]=u-knot[i+1-j];
     right[j]=knot[i+j]-u;
     saved=0.0;
     for( r=0; r<j; r++ )
     {
       ndu[j][r]=right[r+1]+left[j-r];
       temp=ndu[r][j-1]/ndu[j][r];
      
       ndu[r][j]=saved+right[r+1]*temp;
       saved=left[j-r]*temp;
     }
     ndu[j][j]=saved;
   }
   for( j=0; j<=p; j++ )
     ders[0][j]=ndu[j][p];  
    
   if( order==0 )
     return;
  

   for( r=0; r<=p; r++ )
   {
     int s1=0, s2=1;   
     a[0][0]=1.0;

     for( k=1; k<=order; k++ )
     {
       double d=0.;
       int rk=r-k, pk=p-k;
       if( r>=k )
       {
 			a[s2][0]=a[s1][0]/ndu[pk+1][rk];
 			d=a[s2][0]*ndu[rk][pk];
       }
       j1 = rk >= -1 ? 1 : -rk;
       j2 = (r-1<=pk) ? k-1 : p-r;
       for( j=j1; j<=j2; j++ )
       {
           a[s2][j]=(a[s1][j]-a[s1][j-1])/ndu[pk+1][rk+j];
           d+=a[s2][j]*ndu[rk+j][pk];
       }
       if( r<=pk )
       {
           a[s2][k]= -a[s1][k-1]/ndu[pk+1][r];
           d+=a[s2][k]*ndu[r][pk];
       }
       ders[k][r]=d;
       j=s1; s1=s2; s2=j;  
     }
   }
   r=p;
   for( k=1; k<=order; k++ )
   {
     for( j=0; j<=p; j++ ) 
       ders[k][j]*=r;
     r*=(p-k);
   }
   
   free(left); 
   free(right);
   
   free2Darray(ndu, p+1);
   free2Darray(a, p+1);
}

void igafem::nurbs2DBasisDers(double u1, double u2, unsigned int p, unsigned int q, 
                              const std::vector<double>& knotU, const std::vector<double>& knotV, const std::vector<double>& weight, 
                              Eigen::VectorXd& R, Eigen::VectorXd& dRdxi, Eigen::VectorXd& dRdet)
//void nurbs2DBasisDers(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Return the NURBS basis functions and first derivatives to matlab
    // All non-zero basis functions and derivatives at point [xi,eta] are
    // computed.
    //
    // We expect the function to be called as 
    // [R dRdxi dRdeta] = NURBSinterpolation(xi, p, q, knotU, ...
    //                        knotV, weights)
    //
    //	xi           = point, [xi eta], where we want to interpolate
    //	knotU, knotV = knot vectors
    //	weights      = vector of weights 
    //
    // Vinh Phu Nguyen, nvinhphu@gmail.com
    */
    
    /* First get the inputs */
    
    int    numKnotU = knotU.size();
    int    numKnotV = knotV.size();
    
    int    nU       = numKnotU - 1 - p - 1;
    int    nV       = numKnotV - 1 - q - 1;
    int    noFuncs  = (p+1)*(q+1); 
    
    //int    numWeights = weight.size();
        
    //mexPrintf("\n\np=%d\n", numWeights);
  
    double tol    = 100*DBL_EPSILON;
    
    if(fabs(u1 - knotU[numKnotU-1]) < tol) 
        u1 = (knotU[numKnotU-1] - tol);
    
    if(fabs(u2 - knotV[numKnotV-1]) < tol) 
        u2 = (knotV[numKnotV-1] - tol); 
    
    /* and evaluate the non-zero univariate B-spline basis functions
     * and first derivatives 
     */

    Eigen::VectorXd N(p+1);
    Eigen::VectorXd M(q+1);
    //Eigen::MatrixXd dersN(nU+1, p+1);
    //Eigen::MatrixXd dersM(nV+1, q+1);
    
    //double *N      = (double *)malloc(sizeof(double)*(p+1));
    //double *M      = (double *)malloc(sizeof(double)*(q+1));
    double **dersN = init2DArray(nU+1, p+1);
    double **dersM = init2DArray(nV+1, q+1);
    
    int spanU = findspan(nU, p, u1, knotU); 
    int spanV = findspan(nV, q, u2, knotV); 
    
    BasisFuns     (spanU, u1, p, knotU, N);
    BasisFuns     (spanV, u2, q, knotV, M);
    //std::cout << "----" << std::endl;
    //std::cout << spanU << "," << u1 << "," << p << std::endl;
    //std::cout << N.transpose() << std::endl;
    //std::cout << M.transpose() << std::endl;
    //std::cout << "----" << std::endl;

    dersBasisFuns (spanU, u1, p, nU, const_cast<double*>(knotU.data()), dersN);	
    dersBasisFuns (spanV, u2, q, nV, const_cast<double*>(knotV.data()), dersM);
    //dersBasisFuns(int i, double u, int p, int order, double knot[], Eigen::MatrixXd& ders)
    //std::cout << N << std::endl;
    //std::cout << M << std::endl;
    //std::cout << dersN << std::endl;
    //std::cout << dersM << std::endl;
    
    //for(size_t i=0;i<=nU;i++){
    //    for (size_t j = 0; j<=p; j++) {
    //        printf("dersN= %f\n", dersN[i][j]);
    //    //printf("dNdet= %f\n", dersM[1][i]);
    //    }
    //}
    
    //printf("vind = %d\n", numWeights);printf("vind = %f\n", xi[0]);
    
    /* and create NURBS approximation */
    
    int i, j, k, c;
    
    /*
    for(i=0;i<=p;i++){
        printf("dNdxi= %f\n", dersN[1][i]);
        printf("dNdet= %f\n", dersM[1][i]);
    }*/
    
    int uind = spanU - p;
    int vind;
    
    double w     = 0.0;
    double dwdxi = 0.0;
    double dwdet = 0.0;
    double wgt;
    
    c = 0;
    for(j = 0; j <= q; j++) {
        vind = spanV - q + j;  
        for(i = 0; i <= p; i++)
        {               
            //c   = uind + i + vind * (nU+1);            
            wgt = weight[c++];
            
            w     += N[i]        * M[j] * wgt;
            dwdxi += dersN[1][i] * M[j] * wgt;
            dwdet += dersM[1][j] * N[i] * wgt;
        }
    }

    /* create output */
    
    //plhs[0] = mxCreateDoubleMatrix(1,noFuncs,mxREAL); 
    //plhs[1] = mxCreateDoubleMatrix(1,noFuncs,mxREAL); 
    //plhs[2] = mxCreateDoubleMatrix(1,noFuncs,mxREAL); 
    //
    //double *R      = mxGetPr(plhs[0]);
    //double *dRdxi  = mxGetPr(plhs[1]);
    //double *dRdet  = mxGetPr(plhs[2]);
    
    uind = spanU - p;
    k    = 0;
    
    double fac;
    
    /*printf("uind= %d\n", uind);  
    printf("vind= %d\n", vind);  
    printf("nU+1= %d\n", nU+1);  */
   
    c = 0;
    for(j = 0; j <= q; j++) {
        vind = spanV - q + j; 
        for(i = 0; i <= p; i++) {               
            //c        = uind + i + vind*(nU+1);
            fac      = weight[c++]/(w*w);
            
            R[k]     = N[i]*M[j]*fac*w;
            dRdxi[k] = (dersN[1][i]*M[j]*w - N[i]*M[j]*dwdxi) * fac;
            dRdet[k] = (dersM[1][j]*N[i]*w - N[i]*M[j]*dwdet) * fac;
            
            k += 1;
        }
    }
    
        /*mexPrintf("\nWe have a knot vector with %d components\n", numWeights);
                for(k = 0; k < numWeights; k++) mexPrintf("%2.2f\t", weight[k]);
                mexPrintf("\n");*/
                
    //free(N);
    //free(M);
    free2Darray(dersN, (nU+1));
    free2Darray(dersM, (nV+1));	
}


void igafem::nurbs3DBasisDers(double u1, double u2, double u3, 
                      unsigned int p, unsigned int q, unsigned int r,
                      const std::vector<double> &knotU,
                      const std::vector<double> &knotV,
                      const std::vector<double> &knotW,
                      const std::vector<double> &weight, Eigen::VectorXd &R,
                      Eigen::VectorXd &dRdxi, Eigen::VectorXd &dRdet, Eigen::VectorXd &dRdze)
{
    /* Return the 3D NURBS basis functions and first derivatives to matlab
     * // All non-zero basis functions and derivatives at point [xi,eta,zeta] are
     * // computed.
     * //
     * // We expect the function to be called as
     * // [R dRdxi dRdeta dRdzeta] = NURBS3DBasisDers(xi,p,q,r, knotU, ...
     * //                                               knotV,knotZ, weights)
     * //
     * //	xi                  = point, [xi eta zeta], where we want to interpolate
     * //	knotU, knotV, knotZ = knot vectors
     * //	weights             = vector of weights
     * //
     * // Vinh Phu Nguyen, nvinhphu@gmail.com
     */
    
    std::array<double, 3> xi; 
    xi[0] = u1; xi[1] = u2; xi[2] = u3;

    int    numKnotU = knotU.size();
    int    numKnotV = knotV.size();
    int    numKnotW = knotW.size();
    
    int    nU       = numKnotU - 1 - p - 1;
    int    nV       = numKnotV - 1 - q - 1;
    int    nW       = numKnotW - 1 - r - 1;
    int    noFuncs  = (p+1)*(q+1)*(r+1);
    
    //int    numWeights = mxGetN(prhs[7]);
    
    //mexPrintf("\n\np=%d\n", numWeights);
    
    double tol    = 100*DBL_EPSILON;
    
    if(fabs(xi[0]-knotU[numKnotU-1]) < tol)
        xi[0] = knotU[numKnotU-1] - tol;
    
    if(fabs(xi[1]-knotV[numKnotV-1]) < tol)
        xi[1] = knotV[numKnotV-1] - tol;
    
    if(fabs(xi[2]-knotW[numKnotW-1]) < tol)
        xi[2] = knotW[numKnotW-1] - tol;
    
    /* and evaluate the non-zero univariate B-spline basis functions
     * and first derivatives
     */
    
    Eigen::VectorXd N(p+1);
    Eigen::VectorXd M(q+1);
    Eigen::VectorXd P(r+1);
    
    double **dersN = init2DArray(nU+1, p+1);
    double **dersM = init2DArray(nV+1, q+1);
    double **dersP = init2DArray(nW+1, r+1);
    
    int spanU = igafem::findspan(nU, p, xi[0], knotU);
    int spanV = igafem::findspan(nV, q, xi[1], knotV);
    int spanW = igafem::findspan(nW, r, xi[2], knotW);
    
    BasisFuns     (spanU, xi[0], p, knotU, N);
    BasisFuns     (spanV, xi[1], q, knotV, M);
    BasisFuns     (spanW, xi[2], r, knotW, P);
    
    dersBasisFuns (spanU, xi[0], p, nU, const_cast<double*>(knotU.data()), dersN);
    dersBasisFuns (spanV, xi[1], q, nV, const_cast<double*>(knotV.data()), dersM);
    dersBasisFuns (spanW, xi[2], r, nW, const_cast<double*>(knotW.data()), dersP);
    
    //printf("vind = %d\n", p);
    //printf("vind = %f\n", xi[0]);
    
    /* and create NURBS approximation */
    
    int i, j, k, c, kk;
    
    /*
     * for(i=0;i<=p;i++){
     * printf("dNdxi= %f\n", dersN[1][i]);
     * printf("dNdet= %f\n", dersM[1][i]);
     * }*/
    
    int uind = spanU - p;
    int vind, wind;
    
    double w     = 0.0;
    double dwdxi = 0.0;
    double dwdet = 0.0;
    double dwdze = 0.0;
    double wgt;
    
		c = 0;
		for (k = 0; k <= r; k++)
		{
			wind = spanW - r + k;

			for (j = 0; j <= q; j++)
			{
				vind = spanV - q + j;

				for (i = 0; i <= p; i++)
				{
					//c   = uind + i + (nU+1) * ( (nV+1)*wind + vind);
					wgt = weight[c++];

					w += N[i] * M[j] * P[k] * wgt;
					dwdxi += dersN[1][i] * M[j] * P[k] * wgt;
					dwdet += dersM[1][j] * N[i] * P[k] * wgt;
					dwdze += dersP[1][k] * N[i] * M[j] * wgt;
				}
			}
		}

    //printf("vind after = %d\n", numWeights);printf("vind = %f\n", xi[0]);

    /* create output */
    
    uind = spanU - p;
    kk   = 0;
    
    double fac;
    double nmp;
    
    /*printf("uind= %d\n", uind);
     * printf("vind= %d\n", vind);
     * printf("nU+1= %d\n", nU+1);  */
    
		c = 0;
		for (k = 0; k <= r; k++)
		{
			wind = spanW - r + k;

			for (j = 0; j <= q; j++)
			{
				vind = spanV - q + j;

				for (i = 0; i <= p; i++)
				{
					//c        = uind + i + (nU+1) * ( (nV+1)*wind + vind);
					fac = weight[c++] / (w*w);
					nmp = N[i] * M[j] * P[k];

					R[kk] = nmp * fac * w;
					dRdxi[kk] = (dersN[1][i] * M[j] * P[k] * w - nmp*dwdxi) * fac;
					dRdet[kk] = (dersM[1][j] * N[i] * P[k] * w - nmp*dwdet) * fac;
					dRdze[kk] = (dersP[1][k] * N[i] * M[j] * w - nmp*dwdze) * fac;

					kk += 1;
				}
			}
		}

    /*mexPrintf("\nWe have a knot vector with %d components\n", numWeights);
     * for(k = 0; k < numWeights; k++) mexPrintf("%2.2f\t", weight[k]);
     * mexPrintf("\n");*/
    
    free2Darray(dersN, (nU+1));
    free2Darray(dersM, (nV+1));
    free2Darray(dersP, (nW+1));
}


void igafem::nurbs1DBasisDers(double u1, unsigned int p, const std::vector<double> &knotU,
                      const std::vector<double> &weight, Eigen::VectorXd &R,
                      Eigen::VectorXd &dRdxi)
{
    /* Return the univariate NURBS basis functions and first derivatives to matlab
     * // All non-zero basis functions and derivatives at point [xi] are
     * // computed.
     * //
     * // We expect the function to be called as
     * // [R dRdxi]      = NURBS1DBasisDers(xi,p,knotU,weights)
     * //
     * //	xi           = where we want to interpolate
     * //	knotU        = knot vector
     * //	weights      = vector of weights
     * //
     * // Vinh Phu Nguyen, nvinhphu@gmail.com
     */
    
    //if(nrhs != 4) mexErrMsgTxt("You fool! You haven't passed in 6 arguments to the function."
    //        "We expect it to be in the form [dRdxi dRdeta] = NURBSinterpolation(xi, p, q, knotU, knotV, weights)\n");
    
    /* First get the inputs */
    
    int    numKnotU = knotU.size();
    int    nU       = numKnotU - 1 - p - 1;
    int    noFuncs  = (p+1);
    
    //mexPrintf("\n\np=%d\n", numWeights);
    
    double tol    = 100*DBL_EPSILON;
    
    if(fabs(u1-knotU[numKnotU-1]) < tol)
        u1 = knotU[numKnotU-1] - tol;
    
    /* and evaluate the non-zero univariate B-spline basis functions
     * and first derivatives
     */
    
    Eigen::VectorXd N(p+1);
    double **dersN = init2DArray(nU+1, p+1);
    
    int    spanU   = igafem::findspan(nU, p, u1, knotU);
    
    BasisFuns     (spanU, u1, p, knotU, N);
    dersBasisFuns (spanU, u1, p, nU, const_cast<double*>(knotU.data()), dersN);
    
    //printf("vind = %d\n", numWeights);printf("vind = %f\n", xi[0]);
    
    /* and create NURBS approximation */
    
    int i;
    
    /*
     * for(i=0;i<=p;i++){
     * printf("dNdxi= %f\n", dersN[1][i]);
     * printf("dNdet= %f\n", dersM[1][i]);
     * }*/
    
    int    uind  = spanU - p;
    double w     = 0.0;
    double dwdxi = 0.0;
    double wgt;
    
    int c = 0;
    for(i = 0; i <= p; i++)
    {       
        wgt    = weight[c++];        
        w     += N[i]        * wgt;
        dwdxi += dersN[1][i] * wgt;        
    }
        
    /* create output */
    
    double fac;
    
    /*printf("uind= %d\n", uind);
     * printf("vind= %d\n", vind);
     * printf("nU+1= %d\n", nU+1);  */
    c=0;
    for(i = 0; i <= p; i++)
    {
        fac      = weight[c++]/(w*w);
        
        R[i]     = N[i]*fac*w;
        dRdxi[i] = (dersN[1][i]*w - N[i]*dwdxi) * fac;                       
    }
       
    free2Darray(dersN, (nU+1));
}