#include"Matrix.H"
#include <math.h>
#include <map>
#include <exception>
#include <stdexcept>

using namespace std;

Matrix::Matrix()
{
	gsl_set_error_handler((void (*)(const char*, const char*, int, int))&throwExp);
	matrix=NULL;
	row=0;
	col=0;
}

Matrix::Matrix(int r,int c) : 
	row(r),
	col(c)
{
	gsl_set_error_handler((void (*)(const char*, const char*, int, int))&throwExp);
	matrix=gsl_matrix_alloc(row,col);
	gsl_matrix_set_zero(matrix);
}

Matrix::~Matrix()
{
	if(matrix!=NULL)
	{
		gsl_matrix_free(matrix);
		matrix=NULL;
	}
}

int 
Matrix::init(int r,int c)
{
	row=r;
	col=c;
	if(matrix!=NULL)
	{
		gsl_matrix_free(matrix);
	}
	matrix=gsl_matrix_alloc(r,c);
	gsl_matrix_set_zero(matrix);
	return 0;
}

int
Matrix::initAsIdentity()
{
	if((row==0)||(col==0))
	{
		cout << "Warning!! Row or col = 0" << endl;
	}
	gsl_matrix_set_identity (matrix);
	return 0;
}
	
//Add matrix b to this matrix
//Return NULL on failure
//Caller must free this matrix
Matrix* 
Matrix::addMatrix(Matrix* b)
{
	if(!dimequal(b))
	{
		cout << "addMatrix, Dimensions do not match" << endl;
		return NULL;	
	}
	Matrix* res=new Matrix(row,col);
	gsl_matrix_memcpy (res->matrix, matrix);
	gsl_matrix_add(res->matrix,b->matrix);
	return res;	
}

Matrix* 
Matrix::subtractMatrix(Matrix* b)
{
	if(!dimequal(b))
	{
		cout << "subtractMatrix, Dimensions do not match" << endl;
		return NULL;	
	}
	Matrix* res=new Matrix(row,col);
	gsl_matrix_memcpy (res->matrix, matrix);
	gsl_matrix_sub(res->matrix,b->matrix);
	return res;	
}

//Multiply this with b
Matrix* 
Matrix::multiplyMatrix(Matrix* b)
{
	if(col!=b->getRowCnt())
	{
		return NULL;
	}
	Matrix* res=new Matrix(row,b->getColCnt());
	gsl_matrix_set_zero (res->matrix);

/* //original
	gsl_matrix_float* resmatrix=gsl_matrix_float_alloc(row,b->getColCnt());
	convertToFloat(resmatrix,res->matrix,row,b->getColCnt());

	gsl_matrix_float* cmatrix=gsl_matrix_float_alloc(row,col);
	convertToFloat(cmatrix,matrix,row,col);

	gsl_matrix_float* bmatrix=gsl_matrix_float_alloc(b->getRowCnt(),b->getColCnt());
	convertToFloat(bmatrix,b->matrix,b->getRowCnt(),b->getColCnt());
	
	gsl_blas_sgemm (CblasNoTrans, CblasNoTrans, 1, cmatrix, bmatrix, 0, resmatrix);
	convertFromFloat(resmatrix,res->matrix,row,b->getColCnt());
	gsl_matrix_float_free(resmatrix);
	gsl_matrix_float_free(cmatrix);
	gsl_matrix_float_free(bmatrix);
	return res;	
*/

	// changed by soyoun: to multiply without converting to/from float matrix
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, matrix, b->matrix, 0, res->matrix);
	return res;
}

int 
Matrix::addWithMatrix(Matrix* b)
{
	if(!dimequal(b))
	{
		cout << "addWithMatrix, Dimensions do not match" << endl;
		return -1;	
	}
	gsl_matrix_add(matrix,b->matrix);
	return 0;	
}
	
int 
Matrix::subtractWithMatrix(Matrix* b)
{
	if(!dimequal(b))
	{
		cout << "subtractWithMatrix, Dimensions do not match" << endl;
		return -1;	
	}
	gsl_matrix_sub(matrix,b->matrix);
	return 0;	
}

int
Matrix::setMultiplyMatrix(Matrix* a, Matrix* b)
{
	if (row != a->row || col != b->col)
	{
		init(a->row,b->col);
	}
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, a->matrix, b->matrix, 0, matrix);
	return 0;
}

int 
Matrix::multiplyWithMatrix(Matrix* b)
{
	if(col!=b->getRowCnt())
	{
		return -1;
	}

	gsl_matrix* res=gsl_matrix_alloc(row,b->getColCnt());
	gsl_matrix_set_zero(res);
	
/* // original	
	gsl_matrix_float* resmatrix=gsl_matrix_float_alloc(row,b->getColCnt());
	convertToFloat(resmatrix,res,row,b->getColCnt());

	gsl_matrix_float* cmatrix=gsl_matrix_float_alloc(row,col);
	convertToFloat(cmatrix,matrix,row,col);

	gsl_matrix_float* bmatrix=gsl_matrix_float_alloc(b->getRowCnt(),b->getColCnt());
	convertToFloat(bmatrix,b->matrix,b->getRowCnt(),b->getColCnt());
	
	gsl_blas_sgemm(CblasNoTrans, CblasNoTrans, 1, cmatrix, bmatrix, 0, resmatrix);
	convertFromFloat(resmatrix,res,row,b->getColCnt());

	gsl_matrix_float_free(resmatrix);
	gsl_matrix_float_free(cmatrix);
	gsl_matrix_float_free(bmatrix);
	
	gsl_matrix_memcpy(matrix, res);
	gsl_matrix_free(res);
	return 0;
*/
	// changed by soyoun : 
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, matrix, b->matrix, 0, res);
	gsl_matrix_memcpy(matrix, res);
	gsl_matrix_free(res);
	return 0;
}
	
int 
Matrix::addScalar(double aVal)
{
	gsl_matrix_add_constant(matrix,aVal);
	return 0;
}

int 
Matrix::subtractScalar(double aVal)
{
	gsl_matrix_add_constant(matrix,-1*aVal);
	return 0;
}

int 
Matrix::multiplyScalar(double aVal)
{
	gsl_matrix_scale(matrix,aVal);
	return 0;
}

int 
Matrix::divideScalar(double aVal)
{
	if(aVal==0)
	{
		return -1;
	}
	gsl_matrix_scale(matrix,1/aVal);
	return 0;
}

int 
Matrix::setValue(double val,int i,int j)
{
	gsl_matrix_set(matrix, i, j, val);
	return 0;
}

int 
Matrix::setAllValues(double val)
{
	gsl_matrix_set_all (matrix,val);
	return 0;
}

double
Matrix::getValue(int i,int j)
{
	double val=gsl_matrix_get(matrix, i, j);
	return val;
}

Matrix* 
Matrix::invMatrix() 		
{
	Matrix* minv=new Matrix(row,col);

	gsl_matrix* ludecomp=gsl_matrix_alloc(row,col);
	gsl_matrix_memcpy(ludecomp, matrix);
	//cout << "Old Value : " << gsl_matrix_get(matrix,0,0) << endl;
	//cout << "New Value : " << gsl_matrix_get(ludecomp,0,0) << endl;
	gsl_permutation* p=gsl_permutation_alloc(row);
	int signum=0;

	gsl_linalg_LU_decomp(ludecomp, p, &signum);
	gsl_linalg_LU_invert(ludecomp, p, minv->matrix);
	gsl_matrix_free(ludecomp);
	gsl_permutation_free(p);
	return minv;
}

Matrix*
Matrix::invMatrix(gsl_matrix* ludecomp, gsl_permutation* p) 		
{
	//cout << "Old Value : " << gsl_matrix_get(matrix,0,0) << endl;
	//cout << "New Value : " << gsl_matrix_get(ludecomp,0,0) << endl;
	Matrix* minv=new Matrix(row,col);
	ludecomp->size1=row;
	ludecomp->size2=col;
	p->size=row;
	gsl_matrix_memcpy(ludecomp, matrix);
	int signum=0;
	gsl_linalg_LU_decomp(ludecomp, p, &signum);
	gsl_linalg_LU_invert(ludecomp, p, minv->matrix);
	return minv;
}

Matrix*
Matrix::transMatrix()
{
	Matrix* transMatrix=new Matrix(col,row);
	gsl_matrix_transpose_memcpy(transMatrix->matrix,matrix);
	return transMatrix;
}

bool 
Matrix::dimequal(Matrix* aMatrix)
{
	if((aMatrix->getRowCnt()==row)&&
	   (aMatrix->getColCnt()==col))
	{
		return true;
	}	
	return false;
}

int 
Matrix::getRowCnt()
{
	return row;
}

int 
Matrix::getColCnt()
{
	return col;
}

double
Matrix::detMatrix()
{
	gsl_matrix* ludecomp=gsl_matrix_alloc(row,col);
	gsl_matrix_memcpy(ludecomp, matrix);
	gsl_permutation* p=gsl_permutation_alloc(row);
	int signum=0;

	gsl_linalg_LU_decomp(ludecomp, p, &signum);
	double det=gsl_linalg_LU_det(ludecomp, signum);
	gsl_matrix_free(ludecomp);
	gsl_permutation_free(p);
	return det;
}

double
Matrix::detMatrix(gsl_matrix* ludecomp, gsl_permutation* p)
{
	ludecomp->size1=row;
	ludecomp->size2=col;
	p->size=row;
	gsl_matrix_memcpy(ludecomp, matrix);
	int signum=0;

	gsl_linalg_LU_decomp(ludecomp, p, &signum);
	double det=gsl_linalg_LU_det(ludecomp, signum);
	return det;
}


/* // not used
int 
Matrix::convertToFloat(gsl_matrix_float* dest,const gsl_matrix* source, int r, int c)
{
	for(int i=0;i<r;i++)
	{
		for(int j=0;j<c;j++)
		{
			float val=gsl_matrix_get(source,i,j);	
			gsl_matrix_float_set(dest,i,j,val);
		}
	}
	return 0;
}

int 
Matrix::convertFromFloat(const gsl_matrix_float* source, gsl_matrix* dest,int r, int c)
{
	for(int i=0;i<r;i++)
	{
		for(int j=0;j<c;j++)
		{
			float val=gsl_matrix_float_get(source,i,j);
			gsl_matrix_set(dest,i,j,val);
		}
	}
	return 0;
}

*/

int
Matrix::showMatrix(ostream& o)
{
	for(int i=0;i<row;i++)
	{
		for(int j=0;j<col;j++)
		{
			o << getValue(i,j) << " ";
		}	
		o << endl;
	}
	return 0;
}

//Normalize column wise
int 
Matrix::normalize()
{
	double total=0;
	for(int i=0;i<col;i++)
	{
		for(int j=0;j<row;j++)
		{
			double val=getValue(j,i);
			total=total + val;
		}
	}	
	
	for(int i=0;i<col;i++)
	{
		for(int j=0;j<row;j++)
		{
			double val=getValue(j,i);
			val=val/total;	
			setValue(val,j,i);
		}	
	}
	return 0;
}

int
Matrix::normalizeVector()
{
	if(col>1)
	{
		return 0;
	}	
	double minValue=gsl_matrix_min(matrix);
	if(minValue<0)
	{
		minValue=-1*minValue;
		for(int i=0;i<row;i++)
		{
			double val=getValue(i,0);
			if(val<0)
			{
				val=val+minValue;
				setValue(val,i,0);	
			}
		}
	}

	double total=0;
	for(int i=0;i<row;i++)
	{
		total=total+getValue(i,0);
	}

	for(int i=0;i<row;i++)
	{
		double val=getValue(i,0)/total;
		setValue(val,i,0);
	}

	return 0;
}

double
Matrix::getMax()
{
	double val=gsl_matrix_max(matrix);
	return val;
}

int
Matrix::makeUncorrelated()
{
	for(int i=0;i<row;i++)
	{
		for(int j=0;j<col;j++)
		{
			if(i!=j)
			{
				setValue(0,i,j);	
				//setValue(getValue(i,j),i,j);	
				//setValue(getValue(i,j),j,i);	
			}
			else 
			{
				if(getValue(i,j)==0)
				{
					cout << "What the f**k!! Zero variance u duffer!! at " << i << " : " << j <<endl; 
					setValue(0.001,i,j);
				}
			}
		}
	}
	return 0;
}

bool 
Matrix::rowZero()
{
	bool zeroCheck=false;
	for(int i=0;i<row;i++)
	{
		int zeroCol=0;
		for(int j=0;j<col;j++)
		{
			if(getValue(i,j)==0)	
			{
				zeroCol++;
			}
		}
		if(zeroCol==col)
		{
			cout << "Row " << i << " is zero" << endl;
			zeroCheck=true;
		}
	}
	return zeroCheck;
}

bool 
Matrix::colZero()
{
	bool zeroCheck=false;
	for(int i=0;i<col;i++)
	{
		int zeroRow=0;
		for(int j=0;j<row;j++)
		{
			if(getValue(j,i)==0)	
			{
				zeroRow++;
			}
		}
		if(zeroRow==row)
		{
			cout << "Col " << i << " is zero" <<  endl;
			zeroCheck=true;
		}
	}
	return zeroCheck;
}

//Make matrix non negative by making all elements less than zero
//as zero
int 
Matrix::makePositive()
{
	double minValue=gsl_matrix_min(matrix);
	if(minValue>=0)
	{
		//This means matrix is positive
		return 0;
	}

	minValue=minValue*-1;

	for(int i=0;i<row;i++)
	{
		for(int j=0;j<col;j++)
		{
			double value=getValue(i,j);
			value=value+minValue;
			setValue(value,i,j);
		}
	}
	return 0;	
}

Matrix*
Matrix::copyMe()
{
	Matrix* aMatrix=new Matrix(row,col);
	//memcpy dest src
	gsl_matrix_memcpy(aMatrix->matrix,matrix);
	/*
	for(int i=0;i<row;i++)
	{	
		for(int j=0;j<col;j++)
		{
			double val=getValue(i,j);
			aMatrix->setValue(val,i,j);
		}
	}
	*/
	return aMatrix;
}

int
Matrix::copyTo(Matrix* M)
{
	if (M->row != row || M->col != col)
	{
		M->init(row,col);
	}
	//memcpy dest src
	gsl_matrix_memcpy(M->matrix,matrix);
	return 0;
}

int
Matrix::copyFrom(Matrix* M)
{
	if (M->row != row || M->col != col)
	{
		init(M->row,M->col);
	}
	//memcpy dest src
	gsl_matrix_memcpy(matrix,M->matrix);
	return 0;
}

//Finds the closest vector to this vector
//which satisfies the constraint that it is a probability
//vector and also belongs to the same distribution

Matrix*
Matrix::findClosest()
{
	if(col!=1)
	{
		cout << "This is valid only for a vector " << endl;
		return NULL;
	}
	
	bool foundClosest=false;
	int iter=0;
	//Matrix* randVector=copyMe();
	Matrix* randVector=new Matrix(row,1);
	for(int i=0;i<row;i++)
	{
		double val=(double)rand()/RAND_MAX;
		randVector->setValue(val,i,0);
	}
	
	randVector->normalizeVector();
	double dist=0;

	while(iter < MAXITER)
	{
		//getDistance
		dist=getDistance(randVector);
		if(iter==0)
		{
			cout << "Initial Distance : " << dist << endl;
		}
		Matrix* nextClosest=getNextClosest(randVector,dist);
		if(nextClosest!=NULL)
		{
			delete randVector;
			randVector=nextClosest;
		}
		iter++;
	}	
	cout << "Final Distance : " << dist << endl;
	return randVector;
}

Matrix*
Matrix::getNextClosest(Matrix* current,double currDist)
{
	//Generate a vector for a list of possible closest
	//vectors. This is done by adding or subtracting a 
	//small fraction from the current 
	map<int,Matrix*> candidates;
	int id=0;
	double corr=(double)rand()/RAND_MAX;
	if(corr>0.1)
	{
		corr=corr*0.1;
	}
	
	for(int i=0;i<row;i++)
	{
		Matrix* cand1=current->copyMe();
		double val=cand1->getValue(i,0);
		double newVal=val+(val*corr);
		cand1->setValue(newVal,i,0);
		cand1->normalizeVector();
		Matrix* cand2=current->copyMe();
		newVal=val-(val*corr);
		cand2->setValue(newVal,i,0);
		cand2->normalizeVector();
		candidates[id]=cand1;
		id++;
		candidates[id]=cand2;
		id++;
	}

	double minDist=10000;
	int minVec=-1;
	for(int i=0;i<id;i++)
	{
		double dist=getDistance(candidates[i]);	
		if(dist<minDist)
		{
			minDist=dist;
			minVec=i;
		}
	}
	Matrix* bestFound=NULL;
	if(minDist<currDist)
	{
		bestFound=candidates[minVec];
	}
	
	map<int,Matrix*>::iterator anIter;
	for(anIter=candidates.begin();anIter!=candidates.end();anIter++)
	{
		Matrix* temp=anIter->second;
		candidates.erase(anIter);
		if(temp!=bestFound)
		{
			delete temp;
		}
	}
	
	return bestFound;
	
}

double
Matrix::getDistance(Matrix* a)
{
	double dist=0;
	for(int i=0;i<row;i++)
	{
		double val=a->getValue(i,0)-getValue(i,0);
		dist=dist + (val*val);
	}

	return sqrt(dist);
}

int
Matrix::readFromFile(char* inname, int r, int c)
{
	init(r,c);
	FILE* f;
	f = fopen(inname,"r");
	gsl_matrix_fscanf(f,matrix);
	fclose(f);
	return 0;
}

int
Matrix::writeToFile(char* outname)
{
	FILE* f;
	f = fopen(outname,"w");
	//gsl_matrix_fprintf(f,X,"%f");
	for (int i=0;i<row;i++)
	{
		for (int j=0;j<col;j++)
		{
			double v = gsl_matrix_get(matrix, i, j);
			fprintf(f,"%f\t",v);
		}
		fprintf(f,"\n");
	}
	fclose(f);
	return 0;
}

double
Matrix::sumRow(int i)
{
	double res = 0;
	if (i<0 || i>=row)
		return res;
	for (int j=0;j<col;j++)
		res = res+getValue(i,j);
	return res;
}

double
Matrix::sumCol(int j)
{
	double res = 0;
	if (j<0 || j>=col)
		return res;
	for (int i=0;i<row;i++)
		res = res+getValue(i,j);
	return res;
}

int
Matrix::removeRows(vector<int>& rows)
{
	gsl_vector* rr = gsl_vector_alloc(row);
	gsl_vector_set_zero(rr);
	for (int i=0;i<rows.size();i++)
	{
		gsl_vector_set(rr,rows[i],1);
	}
	removeRows(rr);
	gsl_vector_free(rr);
	return 0;
}

int
Matrix::removeCols(vector<int>& cols)
{
	gsl_vector* rr = gsl_vector_alloc(col);
	gsl_vector_set_zero(rr);
	for (int i=0;i<cols.size();i++)
	{
		gsl_vector_set(rr,cols[i],1);
	}
	removeCols(rr);
	gsl_vector_free(rr);
	return 0;
}

int
Matrix::removeRows(gsl_vector* rows)
{
	double rr = gsl_blas_dasum(rows);
	int m = row-(int)rr;
	gsl_matrix* temp = gsl_matrix_alloc(m,col);
	int k=0;
	for (int i=0;i<row;i++)
	{
		double v = gsl_vector_get(rows,i);
		if (v==1)
			continue;
		gsl_vector_view v1 = gsl_matrix_row(matrix,i);
		gsl_vector_view v2 = gsl_matrix_row(temp,k);
		gsl_vector_memcpy(&v2.vector,&v1.vector);

		k++;
	}
	row = m;
	gsl_matrix_free(matrix);
	matrix = temp;
	return 0;
}

int
Matrix::removeCols(gsl_vector* cols)
{
	double rc = gsl_blas_dasum(cols);
	int n = col-(int)rc;
	gsl_matrix* temp = gsl_matrix_alloc(row,n);
	int k=0;
	for (int i=0;i<col;i++)
	{
		double v = gsl_vector_get(cols,i);
		if (v==1)
			continue;
		gsl_vector_view v1 = gsl_matrix_column(matrix,i);
		gsl_vector_view v2 = gsl_matrix_column(temp,k);
		gsl_vector_memcpy(&v2.vector,&v1.vector);

		k++;
	}
	col = n;
	gsl_matrix_free(matrix);
	matrix = temp;
	return 0;
}

gsl_vector_view
Matrix::getRowView(int i)
{
	return gsl_matrix_row(matrix,i);
}

gsl_vector_view
Matrix::getColView(int i)
{
	return gsl_matrix_column(matrix,i);
}

int
Matrix::mldivide(gsl_vector* y, gsl_vector* c)
{
	gsl_matrix* cov;
	gsl_multifit_linear_workspace * work;
	double chisq;
	work = gsl_multifit_linear_alloc (row,col);
	cov = gsl_matrix_alloc (col,col);
	gsl_multifit_linear (matrix, y, c, cov, &chisq, work);
	gsl_matrix_free(cov);
	gsl_multifit_linear_free (work);
	return 0;
}

int
Matrix::setCol(gsl_vector* y, int i)
{
	gsl_vector_view d = gsl_matrix_column(matrix,i);
	gsl_vector_memcpy(&d.vector,y);
	return 0;
}

int
Matrix::setRow(gsl_vector* y, int i)
{
	gsl_vector_view d = gsl_matrix_row(matrix,i);
	gsl_vector_memcpy(&d.vector,y);
	return 0;
}

double
Matrix::getRowMean(int r)
{
	if (r<0 || r>=row)
	{
		cerr << "ERROR: row ID out of range!" << endl;
		return 0;
	}
	double m=0;
	for (int j=0;j<col;j++)
	{
		m=m+getValue(r,j);
	}
	m=m/col;
	return m;
}

int
Matrix::getRowMeanSTD(int r, double& mean, double& std)
{
	if (r<0 || r>=row)
	{
		cerr << "ERROR: row ID out of range!" << endl;
		return 0;
	}
	double m=getRowMean(r);
	double s=0;
	double v=0;
	if (col > 1)
	{
		for (int j=0;j<col;j++)
		{
			v = getValue(r,j);
			s += (v-m)*(v-m);
		}
		s = s/(col-1);
		s = sqrt(s);
	}
	if (s==0)
		s=1;
	mean = m;
	std  = s;
	return 0;
}

int
Matrix::rowStandardize()
{
	for (int i=0;i<row;i++)
	{
		double m=0;
		double s=0;
		double v=0;
		getRowMeanSTD(i,m,s);
		for (int j=0;j<col;j++)
		{
			v = getValue(i,j);
			setValue((v-m)/s,i,j);
		}
	}
	return 0;
}

double
Matrix::getColMean(int c)
{
	if (c<0 || c>=col)
	{
		cerr << "ERROR: column ID out of range!" << endl;
		return 0;
	}
	double m=0;
	for (int j=0;j<row;j++)
	{
		m=m+getValue(j,c);
	}
	m=m/row;
	return m;
}

int
Matrix::getColMeanSTD(int c, double& mean, double& std)
{
	if (c<0 || c>=col)
	{
		cerr << "ERROR: column ID out of range!" << endl;
		return 0;
	}
	double m=getColMean(c);
	double s=0;
	double v=0;
	if (row > 1)
	{
		for (int j=0;j<row;j++)
		{
			v = getValue(j,c);
			s += (v-m)*(v-m);
		}
		s = s/(row-1);
		s = sqrt(s);
	}
	if (s==0)
		s=1;
	mean = m;
	std  = s;
	return 0;
}

int
Matrix::colStandardize()
{
	for (int i=0;i<col;i++)
	{
		double m=0;
		double s=0;
		double v=0;
		getColMeanSTD(i,m,s);
		for (int j=0;j<row;j++)
		{
			v = getValue(j,i);
			setValue((v-m)/s,j,i);
		}
	}
	return 0;
}

double
Matrix::getFNorm()
{
	double res = 0;
	double v   = 0;
	for (int i=0;i<row;i++)
	{
		for (int j=0;j<col;j++)
		{
			v = getValue(i,j);
			res = res+v*v;
		}
	}
	res = sqrt(res);
	return res;
}

Matrix* 
Matrix::dotMultiplyMatrix(Matrix* b)
{
	if(!dimequal(b))
	{
		cout << "dotMultiplyMatrix, Dimensions do not match" << endl;
		return NULL;
	}
	Matrix* a = new Matrix;
	a->init(row,col);
	a->copyFrom(b);
	gsl_matrix_mul_elements(a->matrix,matrix);
	return a;
}

int 
Matrix::dotMultiplyWithMatrix(Matrix* b)
{
	if(!dimequal(b))
	{
		cout << "dotMultiplyWithMatrix, Dimensions do not match" << endl;
		return 0;	
	}
	gsl_matrix_mul_elements(matrix,b->matrix);
	return 0;
}

Matrix*
Matrix::getColMatrix(int c)
{
	if (c>=col || c<0)
	{
		cout << "Invalid column id" << endl;
		return NULL;
	}
	Matrix* a = new Matrix;
	a->init(row,1);
	gsl_matrix_view temp = gsl_matrix_submatrix(matrix,0,c,row,1);
	gsl_matrix_memcpy(a->matrix,&temp.matrix);
	return a;
}

void Matrix::throwExp(const char * reason, const char * file, int line, int gsl_errno)
{
	char c_error[1024];
	sprintf(c_error,"In file %s, line %d\n%s\n%s\n",file,line,reason,gsl_strerror(gsl_errno));
	string err(c_error);
	throw runtime_error(err);
}
		
int 
Matrix::getZeroCols(gsl_vector* cvec)
{
	gsl_vector_set_zero(cvec);
	double v;
	for (int i=0;i<col;i++)
	{
		bool flag=true;
		for (int j=0;j<row;j++)
		{
			v = getValue(j,i);
			if (v!=0)
			{
				flag = false;
				break;
			}
		}
		if (flag)
		{
			gsl_vector_set(cvec,i,1);
		}
	}
	return 0;
}

int 
Matrix::getZeroRows(gsl_vector* rvec)
{
	gsl_vector_set_zero(rvec);
	double v;
	for (int i=0;i<row;i++)
	{
		bool flag=true;
		for (int j=0;j<col;j++)
		{
			v = getValue(i,j);
			if (v!=0)
			{
				flag = false;
				break;
			}
		}
		if (flag)
		{
			gsl_vector_set(rvec,i,1);
		}
	}
	return 0;
}

int 
Matrix::countColEq(int c, double v)
{
	int cnt = 0;
	double cv;
	for (int i=0;i<row;i++)
	{
		cv = getValue(i,c);
		if (cv == v)
		{
			cnt += 1;
		}
	}
	return cnt;
}

int 
Matrix::countRowEq(int r, double v)
{
	int cnt = 0;
	double rv;
	for (int i=0;i<col;i++)
	{
		rv = getValue(r,i);
		if (rv == v)
		{
			cnt += 1;
		}
	}
	return cnt;
}
