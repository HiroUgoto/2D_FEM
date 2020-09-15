FC = ifort
FFLAGS = -O3
#FFLAGS = -CB -traceback -g
LIBS = 

EXEC = FEM

Main = main.o get_mesh.o get_solver.o \
	set_fem_mesh.o set_initial_condition.o set_initial_matrix.o \
	mk_local_matrix.o mk_local_damping_matrix.o shape_functions.o gauss_wgst.o \
	check_inside.o set_source.o set_site.o \
	update_time.o 


OBJS = module_node.o module_fem.o 
MODS = module_node.mod module_fem.mod 

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
	rm -f $(EXEC)

.SUFFIXES:.f90 .mod
