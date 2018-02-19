#ifndef _MATRIX
#define _MATRIX

#include <iostream>
#include <vector>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics_double.h>
#define THRESHOLD 0.000001
#define MAXITER 2000

using namespace std;

class Matrix
{
	public:
		Matrix();
		Matrix(int ,int);
		~Matrix();
		int init(int,int);
		int initAsIdentity();

		Matrix* addMatrix(Matrix*);
		Matrix* subtractMatrix(Matrix*);
		Matrix* multiplyMatrix(Matrix*);

		int addWithMatrix(Matrix*);
		int subtractWithMatrix(Matrix*);
		int multiplyWithMatrix(Matrix*);
		
		int addScalar(double);
		int subtractScalar(double);
		int multiplyScalar(double);
		int divideScalar(double);

		int setValue(double,int,int);
		int setAllValues(double);
		double getValue(int,int);

		Matrix* invMatrix(); 		
		Matrix* invMatrix(gsl_matrix*,gsl_permutation*); 		
		Matrix* transMatrix();

		bool dimequal(Matrix*);
		int getRowCnt();
		int getColCnt();

		double detMatrix();
		double detMatrix(gsl_matrix*,gsl_permutation*);
		int showMatrix(ostream& o=cout);
		int normalize();
		int normalizeVector();
		double getMax();
		int makeUncorrelated();
		bool rowZero();
		bool colZero();
		int makePositive();
		Matrix* copyMe();
		Matrix* findClosest();
		Matrix* getNextClosest(Matrix*,double);
		double getDistance(Matrix*);

		int readFromFile(char* inname, int r, int c);
		int writeToFile(char* outname);

		double sumRow(int i);
		double sumCol(int j);

		int removeRows(gsl_vector* rows);
		int removeCols(gsl_vector* cols);
		int removeRows(vector<int>& rows);
		int removeCols(vector<int>& cols);

		gsl_vector_view getRowView(int i);
		gsl_vector_view getColView(int i);
		
		int mldivide(gsl_vector* y, gsl_vector* c);
		int setCol(gsl_vector* y, int i);
		int setRow(gsl_vector* y, int i);

		double getFNorm();

		int rowStandardize();
		int colStandardize();

		int copyTo(Matrix* M);
		int copyFrom(Matrix* M);
		int setMultiplyMatrix(Matrix* a, Matrix* b);
		Matrix* dotMultiplyMatrix(Matrix* b);
		int dotMultiplyWithMatrix(Matrix* b);

		Matrix* getColMatrix(int c);

		double getColMean(int c);
		double getRowMean(int r);
		int getColMeanSTD(int c, double& mean, double& std);
		int getRowMeanSTD(int r, double& mean, double& std);

		int getZeroCols(gsl_vector* cvec);
		int getZeroRows(gsl_vector* rvec);
		int countColEq(int c, double v);
		int countRowEq(int r, double v);
	private:
		static void throwExp(const char * reason, const char * file, int line, int gsl_errno);
		/* // not used 
		int convertToFloat(gsl_matrix_float*,const gsl_matrix*,int,int);
		int convertFromFloat(const gsl_matrix_float*, gsl_matrix*,int,int);
		*/
		gsl_matrix* matrix;
		int row;
		int col;
};

#endif //_MATRIX
