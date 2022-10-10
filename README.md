# Master's Thesis

## Abstract
The Art Gallery Problem is one of the central problems in computational geometry. By aiming to maximise the area seen by the guards in the Art Gallery Problem, we get a continuous cost function, which allows us to compute a gradient. Using the gradient and other heuristics inspired from neural networks, we will try to solve the Art Gallery Problem practically using gradient descent. Specifically, we aim to use the library CGAL. 
Implementations will show how feasible this approach is. We will visualize the results and discuss the performance of the algorithm on different input shapes and sizes.


## Repository Structure
- `code` contains the source code of the program, some input polygons
- `code/visualisations` contains the visualisation scripts. Check `drawing.py` for the Drawing class
- `doc` contains the thesis LaTeX source code
- `presentation` contains the LaTeX and ipe source code for all the presentations done during the thesis


## How to Run
1. Install dependencies
  - CGAL 5.5
  - GCC 9
  - Python 3
2. Navigate into the `code` directory
3. Compile the program using `cmake .`
4. The executable `main` reads input from the `stdin`. Thus you can pipe input files and pipe the output to a file using `./main < [input] > [output]`
5. If you want to see a visualisation of the output, you can pipe the output file into the `visualisations/draw_guard_path.py` file as `cat [output] | python3 visualisations/draw_guard_path.py`

## Input File Structure
```
[gradient descent step size]
[pull gamma hyperparameter]
[number of edges]
[segments in anti-clockwise order: x1 y1 x2 y2]
[x2 y2 x3 y3]
...
```

## Output File Structure
The log file contains information about:
- `i=` iteration number
- `area=` total area seen
- each guard `g[n]=x y`, where `n` is the number of the guard, and `x y` are its coordinates in the current iteration
- `area[n]=` area seen by guard `n` in the current iteration
