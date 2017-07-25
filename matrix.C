#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <essl.h>



using namespace std;
	template<typename T>
	class Matrix {
		
		private:
			T *MAT;
			unsigned row;
			unsigned col;
			
		public:
		// Constructor (default)
  		Matrix<T>() {};
  		
  		// Destructor
  		~Matrix<T>() {};
  		
  		// Copy constructor
  		Matrix<T>(const Matrix<T>& rhs); 

  		// Constructor (initialize square with diagonal _initial matrix)	
  		Matrix<T>(unsigned _row, unsigned _col, const T& _initial);
  		
  		// Constructor (initialize tri-diagonal matrix with values on main diagonal / sub super diagonal)
	   	Matrix<T>(unsigned _row, unsigned _col, const T& _maindiag, const T& _subsuperdiag);
		
		
		
		// OPERATIONS     		
        // Show matrix in row column form
	    void display()const;
		
	    // Conjugate transpose of the matrix (dagger) 
		Matrix<T> dagger();
		
		// Make matrix elements zero 
		Matrix<T> init_zero();
	    
	    	    
	    // OPERATORS OVERLOADED
	    
	    // Assignment 
        Matrix<T>& operator=(const Matrix<T>& rhs);
  		
	    // Access individual elements                                                                                                                                             
        T& operator()(const unsigned& _row, const unsigned& _col) const{
           return this->MAT[_row*this->col+_col];
        }
		
	    // Inverse matrix unitary operator                                                                                                                                            
        Matrix<T> inverse();
        
	    // Increment matrix
		Matrix<T>& operator+=(const Matrix<T>& rhs);
		
		// Decrement operator matrix
		Matrix<T>& operator-=(const Matrix<T>& rhs);
	    
	    // Matrix minus
	    Matrix<T> operator-(const Matrix<T>& rhs);
	    
		// Martrix product
	    Matrix<T>& operator*=(const Matrix<T>& rhs); 
	    
	    // Matrix times constant product
	    Matrix<T> operator*(const T& lhs); 
		
		// Ax=b solver returns x (double)
		vector<double> operator/(const vector<double>& rhs);
	             	
	};


// Copy constructor
	template<typename T>
  	Matrix<T>::Matrix(const Matrix<T>& rhs){
  			row = rhs.row;
			col = rhs.col;
  		    MAT = rhs.MAT;
  	    
	}
        
// Constructor (initialize diagonal matrix with value _initial)
	template<typename T>																																				
	Matrix<T>::Matrix(unsigned _row, unsigned _col, const T& _initial) {
		row = _row;
        col = _col;		
		unsigned k=0;
  		MAT = new T[row*col];
  			for (unsigned i=0; i<row; i++) {
  				for (unsigned j=0; j<col; j++) {
    			   if (i==j) {
				      	MAT[k] = _initial; 
				      	k++;
					}
				    else {
				      	MAT[k] = 0; 
				        k++;
				    }
    			   
    		    }   
    		}
      
	}
	
// Constructor (initialize tri-diagonal matrix with values on main diagonal / sub super diagonal)
	template<typename T>
	Matrix<T>::Matrix(unsigned _row, unsigned _col, const T& _maindiag, const T& _subsuperdiag){
        row = _row;
        col = _col;		
		unsigned k=0;
  		MAT = new T[row*col];
  
  			for(unsigned i=0;i<row;i++){
  				for(unsigned j=0;j<col;j++){
  					if (i==j){ 
						MAT[k]=_maindiag;
						k++;
					} 
                	else if ((i==j+1)||(i==j-1)){
                   		MAT[k]=_subsuperdiag;
                   		k++;
                	}
					else {
				      	MAT[k] = 0; 
				        k++;
				    } 
				} 
   			}
  
	}
	



// OPERATIONS
// Show matrix in row column form 
	template<typename T>
	void Matrix<T>::display() const{
    unsigned k=0;
     		for (unsigned i=0; i<row; i++) {
    			for (unsigned j=0; j<col; j++) {
      				cout << MAT[k] << " ";
      				k++;
    			} cout << endl;  
    		}
  	} 
 // Conjugate transpose of the square matrix (dagger)  
    
 
	template <>
	Matrix<complex<double> > Matrix<complex<double> >::dagger() {
  		complex<double> tmp;
        unsigned dim = this->row;
  		for (unsigned _row=0; _row<row; _row++) {
  			for (unsigned _col=_row; _col<col; _col++) {
				if (_row==_col){
					this->MAT[_row*dim+_col]=conj(this->MAT[_row*dim+_col]);
				}else{
					tmp = this->MAT[_row*dim+_col];
					this->MAT[_row*dim+_col]=conj(this->MAT[_col*dim+_row]);
					this->MAT[_col*dim+_row]=conj(tmp);
				}
			}
    	}
  		return *this;
	}

	// Make matrix elements zero 
	template<typename T>
	Matrix<T> Matrix<T>::init_zero() {
		unsigned dim = row*col;
		for (unsigned i=0; i<dim; i++){ 
    			this->MAT[i] = 0;
      	} 
		return *this;
	}


	
// OPERATORS OVERLOADED	
/////////////////////////////////////////////////
// Assignment operator
	template<typename T>
	Matrix<T>& Matrix<T>::operator=(const Matrix& rhs){
        row = rhs.row;
		row = rhs.col;
        	for (unsigned i=0; i<row*col; i++) {
  				this->MAT[i] = rhs.MAT[i];   
    		}
    return *this; 		
	}	

// Inverse matrix operator ~A
	template<>
	Matrix<complex<double> > Matrix<complex<double> >::inverse(){
  		complex<double> *work = new complex<double>[row*col];
  		int *index = new int[row];
		int info;
  		int s = row;
  		int LDA = row;
  		int lwork = row*col;

  			zgetrf(s,LDA,this->MAT,LDA,index,info);
  				//cout << info << endl;
  				
  			zgetri(s,this->MAT, LDA, index, work, lwork, info);
  				//cout << info << endl;
  				
  		if(info != 0) {// failure!
    			cout << "Inversion failed! (info = " << info << " ) " << endl;
  		}
   			delete[] work;
  			delete[] index;
	  return *this;
	}
	
// Increment operator matrix
    template<typename T>
	Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& rhs) {                           
    	row = rhs.row;
		col = rhs.col;
		for(unsigned i=0; i<row*col; i++){
			this->MAT[i]+=rhs.MAT[i];
		}
    return *this; // return the result by reference
	}
	
// Decrement operator matrix
	template<typename T>
	Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& rhs) {                           
    	row = rhs.row;
		col = rhs.col;
		for(unsigned i=0; i<row*col; i++){
			this->MAT[i]-=rhs.MAT[i];
		}
    return *this; // return the result by reference
	}


 // Matrix multiplication operator
   	template<>
	Matrix<complex<double> >& Matrix<complex<double> >::operator*=(const Matrix<complex<double> >& rhs) {
		int MatSize = rhs.row;	
		complex<double> *result = new complex<double>[MatSize*MatSize];
    	int LDA = MatSize, LDB = MatSize, LDC = MatSize; 
    	int l=MatSize, m=MatSize, n=MatSize;                    
    	char norm = 'N'; 
     
        zgemul(this->MAT, LDA, &norm, rhs.MAT, LDB, &norm, result, LDC, l, m, n);
              
			for (unsigned i=0; i<row*col; i++) {
    			this->MAT[i] = result[i];
			}
			  delete[] result;
  	return *this;
	}   
	

	
// Matrix const product
	template<typename T>
    Matrix<T> Matrix<T>::operator*(const T& rhs) {
			for (unsigned i=0; i<row*col; i++) {
    			this->MAT[i] = (this->MAT[i])*rhs;
    		}
  		return *this;
		
	}

// Matrix solver Ax=b (double)
	template<>
    vector<double> Matrix<double>::operator/(const vector<double>& rhs) {
		
		int MatSize = this->row;
		double res[MatSize];
		vector<double> result(MatSize,0);
		int LDA = MatSize, LDB = MatSize, iopt = 1; 
        int n = MatSize;                    
        int ipvt[MatSize];
		
			for(unsigned i=0; i<MatSize; ++i){
				res[i] = rhs[i];
				
			}
				dgef(this->MAT, LDA, n, ipvt);  
				
				dges(this->MAT, LDA, n, ipvt, res, iopt);
					
			for(unsigned i=0; i<MatSize; ++i){
				result[i] = res[i];
			}
			
	    return result;       
	}
  		
		

// Types specialization	
template class Matrix<complex<double> >;
template class Matrix<double>;

