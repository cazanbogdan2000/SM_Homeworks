NUM_PROC=4
NUM_THREAD=4

build: canny_edge_detection_seq canny_edge_detection_openmp canny_edge_detection_mpi canny_edge_detection_pthreads canny_edge_detection_mpi_and_opnemp canny_edge_detection_mpi_and_pthreads

canny_edge_detection_seq: canny_edge_detection_seq.c
	gcc canny_edge_detection_seq.c -o canny_edge_detection_seq -lm

canny_edge_detection_openmp: canny_edge_detection_openmp.c
	gcc canny_edge_detection_openmp.c -o canny_edge_detection_openmp -lm -fopenmp

canny_edge_detection_mpi: canny_edge_detection_mpi.c
	mpicc canny_edge_detection_mpi.c -o canny_edge_detection_mpi -lm

canny_edge_detection_pthreads: canny_edge_detection_pthreads.c
	gcc -g canny_edge_detection_pthreads.c -o canny_edge_detection_pthreads -lm -lpthread
	
canny_edge_detection_mpi_and_opnemp: canny_edge_detection_mpi_and_opnemp.c
	mpicc canny_edge_detection_mpi_and_opnemp.c -o canny_edge_detection_mpi_and_opnemp -lm -fopenmp

canny_edge_detection_mpi_and_pthreads: canny_edge_detection_mpi_and_pthreads.c
	mpicc canny_edge_detection_mpi_and_pthreads.c -o canny_edge_detection_mpi_and_pthreads -lm -lpthread

run_seq: canny_edge_detection_seq
	./canny_edge_detection_seq

run_openmp: canny_edge_detection_openmp
	./canny_edge_detection_openmp

run_mpi: canny_edge_detection_mpi
	mpirun -np $(NUM_PROC) ./canny_edge_detection_mpi

run_pthreads:
	./canny_edge_detection_pthreads $(NUM_THREAD)

run_mpi_and_pthreads: canny_edge_detection_mpi_and_pthreads
	mpirun -np $(NUM_PROC) ./canny_edge_detection_mpi_and_pthreads $(NUM_THREAD)

run_mpi_and_opnemp: canny_edge_detection_mpi_and_opnemp
	mpirun -np $(NUM_PROC) ./canny_edge_detection_mpi_and_opnemp $(NUM_THREAD)

diff_mpi: run_mpi run_seq
	diff image_seq.bmp image_mpi.bmp

diff_pthreads: run_pthreads run_seq
	diff image_seq.bmp image_pthreads.bmp

diff_openmp: run_openmp run_seq
	diff image_seq.bmp image_openmp.bmp

diff_mpi_and_openmp: run_mpi_and_opnemp run_seq
	diff image_seq.bmp image_mpi_and_openmp.bmp

diff_mpi_and_pthreads: run_mpi_and_pthreads run_seq
	diff image_seq.bmp image_mpi_and_pthreads.bmp

clean: 
	rm -f canny_edge_detection
	rm -f canny_edge_detection_openmp
	rm -f canny_edge_detection_mpi
	rm -f canny_edge_detection_pthreads
	rm -f canny_edge_detection_mpi_and_opnemp
	rm -f canny_edge_detection_mpi_and_pthreads
	rm -f image_*

