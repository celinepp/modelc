class MatrixUtil
{
    public:
    /**********************
     genCorrelatedDeviates: comoute 4 correlated deviates for Monte Carlo simulation
    [in] const symmetricMatrix R: correlation matrix
         double dt: time step size
         double z[]: array to store correlated deviates
    [out] doubel z[]
    ************************/
        double * genCorrelatedDeviates (const SymmetricMatri R, double dt, double z[])
        {
            int i, j;
            double sum[4] = {0,0};
            double deviate = 0.0        //standard normal deviate
            int m = R.Nrows();          //number of rows in correlation matrix

            std::vector<double> dz;     //vector of correlated deviates
            std::vector<double> eigenValue; //vector of eigenvalues
            std::vector<double> eigenVector[4] //array of vector of eigenvector

            std::vector<double>::iterator eigenvecIter //vector iterator
            double lambda[4] = {0.0};                  //store eigenvalues of correlation matrix R

            double dw [4] = {0.0};                     //stores correlated deviates
            DiagonalMatrix D(m);                       //diagonal matrix
            Matrix V(m,m);                             //m x n matrix
            D = genEigenValue (R);                     //get eigenvalues
            V = genEigenVectors(R);                    //get eigenvectors

            //store eigenvalues
            for (i =0; i<m; i++)
            {
                eigenValue.push_back(D.element(i,i));
                lambda[i] = D.element(i,i);
            }
            //store rows of eigenvectors so that we can compute
            //dz [i] = v[i][1] * sqrt(eigenvalue [1] * dw1 + v[i][2] * sqrt(eigenvalue[2]) * dw2)
            //+.....
            for (i = 0; i< m; i ++)
            {
                for (j = 0; j<m, j++)
                {
                    eigenVector[i].push_back(V.element(i,j));
                }
            }
            srand (0);
            long seed = (long) rand() %100
            long * idum = &seed;

            //generate uncorrelated deviates
            for (i = 0; i < m; i++)
            {
                deviate = util.NormalDeviate (idum);
                dw[i] = deviate * sqrt(dt);
            }
            //generate correlated deviates
            for (i = 0; i < m; i++)
            {
                eigenVecIter = eigenVector[i].begin();
                for (j = 0; j < m; j++)
                {
                    sum[i] += (*eigenVecIter) * sqrt(lambda[j]) * dw[j];
                    eigenVecIter ++;
                }
                z =[i] = sum[i];
            }
            return z
        }
        //other definde methods...
}

#include "newmatap.h"
#include <vector>
#include <math.h>
#include "Constants.h"
#include "StatUtility.h"
class MatrixUtil
{
    public:
    /******************************
    genCorrelatedDeviatesCholesky: compute correlated deviates from
                                    Cholesky decomposition
    [in]: SymmetricMatrix& R     : correlation matrix
          double dt              : step size
          double z[]             : correlated deviates array to be returned
    [out] double z[]             : array of correlated deviates
    *******************************/
        double* genCorrelatedDeviatesCholesky (const symmetricMatrix& R, double dt, double z[])
        {
            int m = R.Nrows()
            int n = R.Ncols()
            Matrix lb(m,n);           //lower banded matrix
            StatUtil util;            
            double deviate = 0.0;      //standard normal deviate
            double dw[4] = {0.0};      //stores deviate * sqrt(dt)
            double sum = 0.0;          
            long seed = 0;
            long* idum = 0;
            int i,j;
            
            lb = Cholesky (R);          //calls Cholesky routine in NEWMAT library

            srand(time(0));             //initialize RNG
            seed = (long) rand() %100;  //generate seed
            idum = &seed;               //store address of seed
            //generate uncorrelated deviates
            for (i = 0; i < m; i++)
            {
                deviate = util.gasdev(idum); //generate normal deviate
                dw[i] = deviate * sqrt(dt)
            }
            //generate correlated deviates
            for (i = 0; i < m; i++)
            {
                sum = 0;
                for (j = 0; j < m; j++)
                {
                    sum += lb.element (i,j) * dw[j]
                }
                z[i] = sum;
            }
            return z;
        }
}
