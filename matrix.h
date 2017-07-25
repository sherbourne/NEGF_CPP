using namespace std;
	template<typename T>
	class Matrix {
		
		private:
			T *MAT;
			unsigned row;
			unsigned col;
			
		public:
		// Constructor (default)
  		Matrix() {};
  		
  		// Destructor
  		~Matrix() {};
  		
  		// Copy constructor
  		Matrix(const Matrix<T>& rhs); 

  		// Constructor (initialize matrix)	
  		Matrix<T>(unsigned _row, unsigned _col, const T& _initial);
  		
  		// Constructor (initialize tri-diagonal matrix with values on main diagonal / sub super diagonal)
	   	Matrix<T>(unsigned _row, unsigned _col, const T& _maindiag, const T& _subsuperdiag);
		
		
		
		// Operations     		
        // Show matrix in row column form
	    void display() const;
		
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
        //Matrix<T>& operator~();
		
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

