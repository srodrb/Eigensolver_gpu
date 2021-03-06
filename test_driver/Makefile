#CUDAROOT = /usr/local/cuda-9.2

# Flags for GPUs with Volta architecture.  Modify cc value as required.
FLAGS = -O3 -mp -pgf90libs -Mcuda=cc70,cuda10.0 -Mlarge_arrays

INCS = -I../lib_eigsolve 
LIBS = -Mcudalib=cublas,cusolver -L${CUDAHOME}/lib64 -lnvToolsExt
LIBS += ../lib_eigsolve/lib_eigsolve.a

# Populate this section for your LAPACK installation. Default here is for MKL.
INCS += -I${INTELHOME}/compiler/include -I${MKLROOT}/include
LIBS += -L${INTELHOME}/compiler/lib/intel64 -liomp5  -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -pgf90libs -lpthread -lm -ldl
#LIBS += -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -pgf90libs -lpthread -lm -ldl
# LIBS += -L${LAPACKROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_pgi_thread -pgf90libs -lpthread -lm -ldl

# To enable comparison with MAGMA, uncomment and/or populate the following lines:
#MAGMAROOT = 
#INCS += -I${MAGMAROOT}/include
#LIBS += -L${MAGMAROOT}/lib -lmagma
#OPTFLAGS += -DHAVE_MAGMA

all: test_zhegvdx test_dsygvdx test_dsyevd test_zheevd

test_zhegvdx: test_zhegvdx.F90 
	pgcc -c wallclock.c
	pgf90 -c ${FLAGS} toolbox.F90
	pgf90 -o test_zhegvdx test_zhegvdx.F90 toolbox.o wallclock.o ${LIBS} ${FLAGS} ${OPTFLAGS} -pgf90libs ${INCS}

test_dsygvdx: test_dsygvdx.F90 
	pgcc -c wallclock.c
	pgf90 -c ${FLAGS} toolbox.F90
	pgf90 -o test_dsygvdx test_dsygvdx.F90 toolbox.o wallclock.o ${LIBS} ${FLAGS} ${OPTFLAGS} -pgf90libs ${INCS}

test_dsyevd: test_dsyevd.F90 
	pgcc -c wallclock.c
	pgf90 -c ${FLAGS} toolbox.F90
	pgf90 -o test_dsyevd test_dsyevd.F90 toolbox.o wallclock.o ${LIBS} ${FLAGS} ${OPTFLAGS} -pgf90libs ${INCS}
	
test_zheevd: test_zheevd.F90 
	pgcc -c wallclock.c
	pgf90 -c ${FLAGS} toolbox.F90
	pgf90 -o test_zheevd test_zheevd.F90 toolbox.o wallclock.o ${LIBS} ${FLAGS} ${OPTFLAGS} -pgf90libs ${INCS}

clean:
	rm test_dsyevd test_zheevd test_zhegvdx test_dsygvdx *.mod *.o
