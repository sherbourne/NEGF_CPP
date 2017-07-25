CXX=mpicxx

SCINET_ESSL_INC  = /opt/ibmmath/essl/5.1/include
SCINET_GSL_INC   = /scinet/bgq/Libraries/gsl-1.15/include
SCINET_ESSL_LIB  = /opt/ibmmath/essl/5.1/lib64
SCINET_GSL_LIB   = /scinet/bgq/Libraries/gsl-1.15/lib
SCINET_XLF_LIB   = /opt/ibmcmp/xlf/bg/14.1/bglib64


CPPFLAGS = -I${SCINET_ESSL_INC} -I${SCINET_GSL_INC}
LDFLAGS  = -L${SCINET_ESSL_LIB} -L${SCINET_XLF_LIB} -L${SCINET_GSL_LIB}
LDLIBS   = -lesslbg -lgsl -lxlf90_r -lxlfmath

all:mainMPI
mainMPI.o: mainMPI.C
		 $(CXX) mainMPI.C -c -O3 -qarch=qp -qtune=qp -I${SCINET_ESSL_INC} -L${SCINET_ESSL_LIB} -lesslbg
matrix.o : matrix.C
		 $(CXX) matrix.C -c -O3 -qarch=qp -qtune=qp -I${SCINET_ESSL_INC} -L${SCINET_ESSL_LIB} -lesslbg
matlib.o : matlib.C
		 $(CXX) matlib.C -c -O3 -qarch=qp -qtune=qp -I${SCINET_GSL_INC} -L${SCINET_GSL_LIB} -lgsl		
mainMPI  : mainMPI.o matrix.o matlib.o
		 $(CXX) mainMPI.o matrix.o matlib.o -o mainMPI -O3 -qarch=qp -qtune=qp -v -qphsinfo $(CPPFLAGS) $(LDFLAGS) $(LDLIBS)
clean:
		rm mainMPI *.o *.dat core.* *.eps
  
      
