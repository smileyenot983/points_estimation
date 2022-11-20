This package estimates coordinates of N points from intersection of M lines.
Number of points is given.

Main steps are following:
1. read line coordinates(defined by 2 points lying on this line) from txt file("1-lines.txt" for example)
2. calculate line coefficients(slope and intercept) for all lines
3. calculate intersection points between all lines(1-M, 2-M), i.e. intersection between first line and all other lines, second and all others etc.
4. build dict of intersections(key - point coordinates, value - number of line intersections at this point) and sort this dict by values.
5. cluster points and get N clusters with highest number of intersections, that would be our points
6. write all that N points to "out.txt"

Several assumptions were made:
1. number of points N is given
2. all point coordinates(x,y) are in range [-1000,1000]
3. rounding point coordinates up to closest integer number


Original implementation was slow enough, it was taking around 14 minutes for 15000 lines and 100 points. Several improvements were done:
0. original algorithm (866 seconds)
1. not considering intersections outside [-1000, 1000] region (517 seconds)
2. using substitution instead of linear algebra for intersection estimation (354 seconds)
3. precomputing intersections beforehand (268 seconds)
4. adding multithreaded(4 threads) calculations in for cycle using openMP (80 seconds)

Exact timings might vary depending on PC.


## Requirements

1. OpenMP(https://learn.microsoft.com/ru-ru/cpp/parallel/openmp/reference/openmp-directives?view=msvc-170)

## Building

1. mkdir build && cd build
2. cmake .. 
3. make

## Running

executable name and path to line coordinates

./estimate_points ../Task-Points/1-lines.txt