This cpp file generates a ppm file (an image format), which shows a fractal similar to the Mandelbrot set, where the potentiation is applied, instead of a polynomial. The mathematical theory uses results from the study of the tetration function. (https://en.wikipedia.org/wiki/Tetration).

-----------------------------------------------------------------------

In the cpp file, several variables can be modified:

itermax  (maximal number of iterations)
escTol (tolerance for when the iteration escapes)
approxTol (tolerance for when two complex numbers are numerically the same)
name (name of the output file)
XLength (number of rows in the output file)
YLength (number of columns in the output file)
rmi (minimal real value)
rma (maximal real value)
imi (minimal imaginary value)
ima (maximal imaginary value)
YSymmetry (If we know that the output will be symmetric)
f (the underlying function; possible examples are commented out)

-----------------------------------------------------------------------

To open the ppm file we suggest GIMP. You might need the proper extension.

-----------------------------------------------------------------------

Things to improve:
 - Concurrency/Parallelization is not optimal.
 - The symmetry appliction might demand to much RAM for large images.
 - Involving the GPU might lead to faster run-time.