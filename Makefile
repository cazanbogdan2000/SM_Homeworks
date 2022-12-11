build: tema1_seq.c tema1_openmp.c
	gcc -g tema1_seq.c -o tema1_seq -lm
	gcc -g tema1_openmp.c -o tema1_openmp -lm -fopenmp
	gcc -g tema1_pthreads.c -o tema1_pthreads -lm -lpthread

run_openmp:
	./tema1_openmp

run_pthreads:
	./tema1_pthreads

run_seq:
	./tema1_seq

clean:
	rm -rf tema1
