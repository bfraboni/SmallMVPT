all: 
	g++ -o smallmvpt ./src/smallmvpt.cxx -O3 -std=c++0x -fopenmp

clean:
	rm smallmvpt

unreport:
	rm *.bmp index.html
