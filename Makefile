build: canny_edge_detection canny_edge_detection_openmp canny_edge_detection_mpi canny_edge_detection_pthreads

canny_edge_detection: canny_edge_detection.c
	gcc canny_edge_detection.c -o canny_edge_detection -lm

canny_edge_detection_openmp: canny_edge_detection_openmp.c
	gcc canny_edge_detection_openmp.c -o canny_edge_detection_openmp -lm -fopenmp

canny_edge_detection_mpi: canny_edge_detection_mpi.c
	mpicc canny_edge_detection_mpi.c -o canny_edge_detection_mpi -lm

canny_edge_detection_pthreads: canny_edge_detection_pthreads.c
	gcc -g canny_edge_detection_pthreads.c -o canny_edge_detection_pthreads -lm -lpthread

run_seq: canny_edge_detection
	./canny_edge_detection

run_openmp: canny_edge_detection_openmp
	./canny_edge_detection_openmp

run_mpi_4: canny_edge_detection_mpi
	mpirun -np 4 ./canny_edge_detection_mpi

run_mpi_8: canny_edge_detection_mpi
	mpirun -np 8 ./canny_edge_detection_mpi
	
run_pthreads_4:
	./canny_edge_detection_pthreads 4
	
run_pthreads_8:
	./canny_edge_detection_pthreads 8

clean: 
	rm -f canny_edge_detection
	rm -f canny_edge_detection_openmp
	rm -f canny_edge_detection_mpi
	rm -f canny_edge_detection_pthreads

