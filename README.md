# Weighted distance transform
A Python implementation for computing a weighted distance transform of an image.

A [distance transform](https://en.wikipedia.org/wiki/Distance_transform) is a map of an image that assign to each pixel its distance to the nearest boundary.
A weighted distance transform extends this by allowing for ... weighted distances. Useful for image analysis or path finding

This script implements a fast marching algorithm in Python with some acceleration in Fortran (if a Fortran compiler is present)

### Prerequisites

 * Python 3 (but no real obstacles for Python 2)
 * Numpy/Scipy and the python imaging library (PIL)
 * Matplotlib for plotting images (not required)
 * Fortran compiler for the acceleration (recommended, not required)

```
Give examples
```
Input: A PNG/JPG image with with boundaries/exits from which to compute the distance (red), obstacles (full black), normal space (white), less 'accessible space' (grey-ish)
Output: A distance transform. Image plotted in Matplotlib
### Installing

Installing uses Distutils with NumPy's extensions

Clone the repository

```bash
git clone https://github.com/0mar/weighted-distance-transform.git
```

If a Fortran compiler is present, install the modules:

```bash
cd weighted-distance-transform
python3 setup.py install
```

Run the script with some example figures (like `images/ex1.png`)

```python
wdt = WDT('images/ex1.png')
wdt.get_weighted_distance_transform()
wdt.plot()
```
End with an example of getting some data out of the system or using it for a little demo

## Authors

* **Omar Richardson**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## References

Algorithms by
* [Tsitsiklis](http://www.mit.edu/~jnt/dijkstra.html)
* [Treuille et al.](http://grail.cs.washington.edu/projects/crowd-flows/78-treuille.pdf)
