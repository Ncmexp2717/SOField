
FFTW = -I$(MKLROOT)/include/fftw
LIB= -L$(MKLROOT)/lib/intel64/ -mkl=parallel -lifcore
CC = mpicc -qopenmp -O3 -xCORE-AVX2 -ip -no-prec-div $(FFTW)

#CC = mpifccpx -Dxt3 -Dkcomp -Kfast,openmp,SPARC64IXfx -I/usr/local/fftw/include
#CC = mpifccpx -Dxt3 -Dkcomp -KSPARC64IXfx -O0 -g -I/usr/local/fftw/include
#LIB = -L/usr/local/fftw/lib64 -lfftw3_omp -lfftw3_mpi -lfftw3 -SSL2MPI -SSL2BLAMP

DEL=rm

sofield: sofield.o Circular_Search.o Band_Dispersion.o read_scfout.o Inputtools.o Tools_BandCalc.o GetOrbital.o EigenValue_Problem.o Eigen_HH.o Tools_Search.o
	$(CC) sofield.o Circular_Search.o Band_Dispersion.o read_scfout.o Inputtools.o  Tools_BandCalc.o GetOrbital.o EigenValue_Problem.o Eigen_HH.o Tools_Search.o  $(LIB) -lm -o sofield

sofield.o: sofield.c Band_Dispersion.h read_scfout.h Inputtools.h Tools_BandCalc.h GetOrbital.h
	$(CC) -c sofield.c

Inputtools.o: Inputtools.c
	$(CC) -c Inputtools.c

read_scfout.o: read_scfout.c read_scfout.h
	$(CC) -c read_scfout.c -DSOField

dtime.o: dtime.c
	$(CC) -c dtime.c

Tools_BandCalc.o: Tools_BandCalc.c  lapack_prototypes.h read_scfout.h Tools_BandCalc.h
	$(CC) -c Tools_BandCalc.c

GetOrbital.o: GetOrbital.c Inputtools.h read_scfout.h Tools_BandCalc.h
	$(CC) -c GetOrbital.c

EigenValue_Problem.o: EigenValue_Problem.c read_scfout.h Tools_BandCalc.h lapack_prototypes.h f77func.h Eigen_HH.h EigenValue_Problem.h
	$(CC) -c EigenValue_Problem.c

Eigen_HH.o: Eigen_HH.c read_scfout.h Tools_BandCalc.h lapack_prototypes.h f77func.h Eigen_HH.h
	$(CC) -c Eigen_HH.c

Tools_Search.o: Tools_Search.c Tools_BandCalc.h EigenValue_Problem.h
	$(CC) -c Tools_Search.c

Band_Dispersion.o: Band_Dispersion.c GetOrbital.h Inputtools.h read_scfout.h Tools_BandCalc.h EigenValue_Problem.h Eigen_HH.h EigenValue_Problem.h
	$(CC) -c Band_Dispersion.c

Circular_Search.o: Circular_Search.c GetOrbital.h Inputtools.h read_scfout.h Tools_BandCalc.h EigenValue_Problem.h Eigen_HH.h EigenValue_Problem.h
	$(CC) -c Circular_Search.c


AMulPop: AMulPop.o Tools_BandCalc.o Inputtools.o
	$(CC) AMulPop.o Tools_BandCalc.o Inputtools.o $(LIB) -lm -o AMulPop

AMulPop.o: AMulPop.c Tools_BandCalc.h Inputtools.h
	$(CC) -c AMulPop.c


MulPCalc: MulPCalc.o Tools_BandCalc.o Inputtools.o
	$(CC) MulPCalc.o Tools_BandCalc.o Inputtools.o $(LIB) -lm -o MulPCalc

MulPCalc.o: MulPCalc.c Tools_BandCalc.h Inputtools.h
	$(CC) -c MulPCalc.c


ADenBand: ADenBand.o Tools_BandCalc.o Inputtools.o
	$(CC) ADenBand.o Tools_BandCalc.o Inputtools.o $(LIB) -lm -o ADenBand

ADenBand.o: ADenBand.c Tools_BandCalc.h Inputtools.h read_scfout.h
	$(CC) -c ADenBand.c


FermiLoop: FermiLoop.o read_scfout.o Inputtools.o Tools_BandCalc.o GetOrbital.o EigenValue_Problem.o Eigen_HH.o 
	$(CC) FermiLoop.o read_scfout.o Inputtools.o  Tools_BandCalc.o GetOrbital.o EigenValue_Problem.o Eigen_HH.o  $(LIB) -lm -o FermiLoop

FermiLoop.o: FermiLoop.c read_scfout.h Inputtools.h EigenValue_Problem.h Eigen_HH.h Tools_BandCalc.h GetOrbital.h
	$(CC) -c FermiLoop.c 

GridCalc: GridCalc.o read_scfout.o Inputtools.o Tools_BandCalc.o GetOrbital.o EigenValue_Problem.o Eigen_HH.o
	$(CC) GridCalc.o read_scfout.o Inputtools.o  Tools_BandCalc.o GetOrbital.o EigenValue_Problem.o Eigen_HH.o  $(LIB) -lm -o GridCalc 

GridCalc.o: GridCalc.c read_scfout.h Inputtools.h EigenValue_Problem.h Eigen_HH.h Tools_BandCalc.h GetOrbital.h
	$(CC) -c GridCalc.c


