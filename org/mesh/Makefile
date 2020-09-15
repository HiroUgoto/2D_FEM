FC = ifort
FFLAGS = 
LIBS = 

EXEC = Mesh

Main = mk_mesh.o

OBJS = parameter_set.o
MODS = parameter_set.mod

all : $(EXEC) 

$(EXEC) : $(MODS) $(OBJS) $(Main)
	$(FC) -o $(EXEC) $(OBJS) $(Main) $(LIBS)

.f90.mod :
	$(FC) $(FFLAGS) -c $<

.f90.o :
	$(FC) $(FFLAGS) -c $< -o $@

clean : 
	rm -f *.o *.mod
	rm -f *~
	rm -f *.in
	rm -f $(EXEC)

.SUFFIXES:.f90 .mod
