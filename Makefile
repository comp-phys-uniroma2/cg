
FC = gfortran
FCOPT = -O2 -fmax-errors=3 
LIBS = -L../../lapack/lapack-3.9.1 -llapack -lrefblas
SOURCES = $(wildcard *.f90)
OBJS = $(SOURCES:.f90=.o)

TARGET = poisson

all: $(TARGET)


clean:
	rm *.o *.mod $(TARGET)

poisson.o: cg.o matdef.o adj_map.o precision.o sparsealg.o gmres.o
matdef.o: precision.o
sparsealg.o: precision.o matdef.o
cg.o: precision.o matdef.o preconditioners.o sparsealg.o
gmres.o: precision.o matdef.o sparsealg.o preconditioners.o

$(TARGET): $(OBJS)
	$(FC) -o $(TARGET) $(OBJS) $(LIBS)


%.o: %.f90
	$(FC) -c $(FCOPT) $<




