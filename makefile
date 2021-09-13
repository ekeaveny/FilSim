FFTLIBS = -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw -lm
InputFiles = fftw_parallel.cpp
FilamentFiles = main.cpp MobilityConstraints.cpp CollisionBarrierFilament.cpp print_functions.cpp file_functions.cpp filament_initialisation_functions.cpp FCMfunctions.cpp c_array_functions.cpp profilers.cpp spring_link_functions.cpp Filament.cpp
FilamentFilesGR = main_gaitreset.cpp MobilityConstraints.cpp CollisionBarrierFilament.cpp print_functions.cpp file_functions.cpp filament_initialisation_functions.cpp FCMfunctions.cpp c_array_functions.cpp profilers.cpp spring_link_functions.cpp Filament.cpp
FilamentFilesGRR = main_gaitresetandremove.cpp MobilityConstraints.cpp CollisionBarrierFilament.cpp print_functions.cpp file_functions.cpp filament_initialisation_functions.cpp FCMfunctions.cpp c_array_functions.cpp profilers.cpp spring_link_functions.cpp Filament.cpp
FilamentFilesTP = main_twoparticle.cpp MobilityConstraints.cpp CollisionBarrierFilament.cpp print_functions.cpp file_functions.cpp filament_initialisation_functions.cpp FCMfunctions.cpp c_array_functions.cpp profilers.cpp spring_link_functions.cpp Filament.cpp
SingleFiles = singlebead.cpp MobilityConstraints.cpp CollisionBarrierFilament.cpp
export OPENBLAS_NUM_THREADS=1
HAM_ARMADILLO_HOME = ./armadillo-10.1.2
HAM_FFTW_HOME = /usr/local/Cluster-Apps/fftw/openmpi/gcc/64/2.1.5
HAM_OPENBLAS_HOME = /ddn/apps/Cluster-Apps/openblas/gcc-8.2.0/0.2.20
DIO_FFTW_HOME = /home/aplm/smdg73/Documents/fftw-2.1.5

mac:
	g++-8 -o MPIFilament $(SingleFiles) -O3 -larmadillo -std=c++14 -DARMA_DONT_USE_WRAPPER -lopenblas -llapack -fopenmp $(FFTLIBS)

MPI:
	mpic++ -o paralleltest parallel.cpp

FFTWMPI:
	mpic++ -o $@ $(InputFiles) $(FFTLIBS)

test:
	mpic++ -o fftw_test fftw_parallel.cpp $(FFTLIBS)

test_dio:
	mpic++ -o main-test2 main-test2.cpp \
	-O3 -larmadillo -std=c++14 -DARMA_DONT_USE_WRAPPER -lopenblas -llapack \
	-I$(DIO_FFTW_HOME)/rfftw/ -L$(DIO_FFTW_HOME)/rfftw/ \
	-I$(DIO_FFTW_HOME)/fftw/ -L$(DIO_FFTW_HOME)/fftw/ -fopenmp \
	-I$(DIO_FFTW_HOME)/mpi/
	# -I$(DIO_FFTW_HOME)/mpi/ -L$(DIO_FFTW_HOME)/mpi/ \
	# $(FFTLIBS)

	# MobilityConstraints.cpp CollisionBarrierFilament.cpp \
	# print_functions.cpp file_functions.cpp \
	# filament_initialisation_functions.cpp FCMfunctions.cpp \
	# c_array_functions.cpp profilers.cpp \

clean:
	rm FFTWMPI

# For gehrig/rizzuto
MPI_FIL:
	mpic++  -o MPIFilament $(FilamentFiles) -O3 -larmadillo -std=c++14 -DARMA_DONT_USE_WRAPPER -lopenblas -llapack -fopenmp $(FFTLIBS)

MPI_FIL_gaitreset:
	mpic++  -o MPIFilament $(FilamentFilesGR) -O3 -larmadillo -std=c++14 -DARMA_DONT_USE_WRAPPER -lopenblas -llapack -fopenmp $(FFTLIBS)

MPI_FIL_gaitresetandremove:
	mpic++  -o MPIFilament $(FilamentFilesGRR) -O3 -larmadillo -std=c++14 -DARMA_DONT_USE_WRAPPER -lopenblas -llapack -fopenmp $(FFTLIBS)

MPI_FIL_twoparticle:
	mpic++  -o MPIFilament $(FilamentFilesTP) -O3 -larmadillo -std=c++14 -DARMA_DONT_USE_WRAPPER -lopenblas -llapack -fopenmp $(FFTLIBS)

MPI_FIL_debug:
	mpic++  -o MPIFilament $(FilamentFiles) -larmadillo -std=c++14 -DARMA_DONT_USE_WRAPPER -lopenblas -llapack -fopenmp $(FFTLIBS)

# For gehrig/rizzuto
fcm_test:
	mpic++  -o fcm_test fcm_test.cpp MobilityConstraints.cpp CollisionBarrierFilament.cpp print_functions.cpp file_functions.cpp filament_initialisation_functions.cpp FCMfunctions.cpp c_array_functions.cpp profilers.cpp -O3 -larmadillo -std=c++14 -DARMA_DONT_USE_WRAPPER -lopenblas -llapack -fopenmp $(FFTLIBS)

# For cx2 Imperial cluster
# Use make_wrapper_cx2.sh instead of this.
MPI_FIL_cx2:
	mpiicpc  -o MPIFilament $(FilamentFiles) -O3 -std=c++14 -DARMA_DONT_USE_WRAPPER -I/rds/general/user/akt12/home/OpenBLAS/include/ -L/rds/general/user/akt12/home/OpenBLAS/lib -lopenblas -llapack -fopenmp -I$(ARMADILLO_HOME)/include -L$(ARMADILLO_HOME)/lib -larmadillo -I$(FFTW_HOME)/include -L$(FFTW_HOME)/lib $(FFTLIBS)

MPI_FIL_cx2_gaitreset:
	mpiicpc  -o MPIFilament $(FilamentFilesGR) -O3 -std=c++14 -DARMA_DONT_USE_WRAPPER -I/rds/general/user/akt12/home/OpenBLAS/include/ -L/rds/general/user/akt12/home/OpenBLAS/lib -lopenblas -llapack -fopenmp -I$(ARMADILLO_HOME)/include -L$(ARMADILLO_HOME)/lib -larmadillo -I$(FFTW_HOME)/include -L$(FFTW_HOME)/lib $(FFTLIBS)

MPI_FIL_cx2_gaitresetandremove:
	mpiicpc  -o MPIFilament $(FilamentFilesGRR) -O3 -std=c++14 -DARMA_DONT_USE_WRAPPER -I/rds/general/user/akt12/home/OpenBLAS/include/ -L/rds/general/user/akt12/home/OpenBLAS/lib -lopenblas -llapack -fopenmp -I$(ARMADILLO_HOME)/include -L$(ARMADILLO_HOME)/lib -larmadillo -I$(FFTW_HOME)/include -L$(FFTW_HOME)/lib $(FFTLIBS)

# For Hamilton Durham cluster
# Use make_wrapper_hamilton.sh instead of this.
MPI_FIL_hamilton:
	mpicxx  -o MPIFilament $(FilamentFiles) -O3 -std=c++14 -DARMA_DONT_USE_WRAPPER -I$(HAM_OPENBLAS_HOME)/include/ -L$(HAM_OPENBLAS_HOME)/lib -lopenblas -llapack -fopenmp -I$(HAM_ARMADILLO_HOME)/include -L$(HAM_ARMADILLO_HOME)/lib -larmadillo -I$(HAM_FFTW_HOME)/include -L$(HAM_FFTW_HOME)/lib $(FFTLIBS)

single:
	mpic++  -o MPIFilament $(SingleFiles) -O3 -larmadillo -std=c++14 -DARMA_DONT_USE_WRAPPER -lopenblas -llapack -fopenmp $(FFTLIBS)


# execute resulting code like this: mpiexec -np 4 ./paralleltest where 4 is the number of processes.


# new:
# 	-rm MultiFilament
# 	g++  -I . -std=c++11 -O3 -fopenmp -msse2 -o MultiFilament main.cpp MobilityConstraints.cpp RPYeigen.cpp
#
# noopenmp:
# 	-rm MultiFilament
# 	g++  -I . -std=c++11 -O3 -msse2 -fstack-protector-all -o MultiFilament main.cpp MobilityConstraints.cpp RPYeigen.cpp
#
# warnings:
# 	-rm MultiFilament
# 	g++  -I . -std=c++11 -O3 -fopenmp  -msse2 -Wall  -o test main.cpp MobilityConstraints.cpp RPYeigen.cpp
#
# debug:
# 	-rm MultiFilament
# 	g++  -g -I . -std=c++11 -O3 -msse2 -fstack-protector-all -o MultiFilament main.cpp MobilityConstraints.cpp RPYeigen.cpp
#
# # arma:
# # 	g++ -I . -std=c++11 -O2 -msse2 -larmadillo -o MultiFilament main.cpp MobilityConstraints.cpp
#
#
# arma2:
# 	g++ main.cpp MobilityConstraints.cpp -o MultiFilamentArma -O0 -g -larmadillo -std=c++14
#
# arma2mac:
# 	clang++  main.cpp MobilityConstraints.cpp -o MultiFilamentArma -O3 -larmadillo  -std=c++14 #-fopenmp=libiomp5
# # -Xpreprocessor -fopenmp -lomp -I"$(brew --prefix libomp)/include" -L"$(brew --prefix libomp)/lib"
# macdebug:
# 	clang++ -std=c++14  main.cpp MobilityConstraints.cpp -o MultiFilamentArma -O0 -fsanitize=address -larmadillo -g ; ASAN_OPTIONS=detect_leaks=1 #./MultiFilamentArma
#
# macGcc:
# 	g++-8  main.cpp MobilityConstraints.cpp CollisionBarrierFilament.cpp -o MultiFilamentArma -O3 -larmadillo  -std=c++14 #-fopenmp
#
# macGccDebug:
# 	g++-8  main.cpp MobilityConstraints.cpp CollisionBarrierFilament.cpp -o MultiFilamentArma -O0 -larmadillo  -std=c++14 -g -Wall   -Wuninitialized #-fopenmp
#
# gcc:
# 	g++ main.cpp MobilityConstraints.cpp CollisionBarrierFilament.cpp -o MultiFilament -O3 -larmadillo -std=c++14 -DARMA_DONT_USE_WRAPPER -lopenblas -llapack #-ffast-math



# Note: // to find memory leakes or un-initialised pointers, run
# // valgrind --leak-check=yes --track-origins=yes ./MultiFilament
# // valgrind --tool=memcheck --track-origins=yes (compile using c++ with option -g -O0)
