CC = g++-7
OUT = main
#LIBS = -lm -L/home/grads/linling2/.linuxbrew/Cellar/libomp/8.0.1/lib/ -lgomp
FLAGS = -Wall -Wno-sign-compare -Ofast -std=gnu++14 -march=native -fopenmp

default:
	$(CC) $(FLAGS) $(INCS) main.cpp -o $(OUT) $(LIBS)

test: test.cpp
	$(CC) $(FLAGS) $(INCS) test.cpp -o test $(LIBS)
