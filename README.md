# Weighted distance transform
Python and Fortran implementation for computing a weighted distance transform of an image.

A [distance transform](https://en.wikipedia.org/wiki/Distance_transform) is a map of an image that assign to each pixel its distance to the nearest boundary.
A weighted distance transform extends this by allowing for _weighted distances_, replacing the uniform Euclidian distance measure with a non-uniform marginal cost function. Useful for image analysis or computing potential functions in path finding.

This script implements a fast marching algorithm in Fortran, compiled to a Python library upon install (if Numpy and a fortran compiler are present). Alternatively, a Python implementation of the same algorithm is provided as fallback.

### Prerequisites

 * Python 3 (but no real obstacles for Python 2)
 * Numpy/Scipy and the Python Imaging library (PIL)
 * Matplotlib for plotting images (not required)
 * Fortran compiler, makes script significantly faster (recommended, not required)

### How does it work?

![Example image](/images/cover_example.png?raw=true "Example of WDT in action")
Input (left): A PNG/JPG image with with boundaries/exits from which to compute the distance (green), obstacles (full black), normal space (white), less 'accessible space' (grey-ish).

Output (right): A distance transform. Image plotted in Matplotlib.
### Installing

Installing uses `setuptools` from distutils with NumPy extensions.

Clone the repository

```bash
git clone https://github.com/0mar/weighted-distance-transform.git
```

If a Fortran compiler is present, install the modules:

```bash
cd weighted-distance-transform
python3 setup.py install
```

### Usage:

The original use case of this script is creating motion planners for 2D environments (in [crowd dynamics](https://symbols.hotell.kau.se/2016/11/30/mercurial/)). We identify goals with _green_, obstacles in _black_, and accessible space in _white_. The corresponding color mapping is included.
You can use any color (except green) to indicate an area that is less easily accessible (and preferably should be avoided). The darker the colour (mapped to a greyscale), the more difficult it is to move through.

Since this mapping is pretty arbitrary, you can provide your own in `wdt.py`.

Once you created an image, you can provide it to the script as shown below. 
`example.py` runs the code below.

```python
import wdt
cost_field = wdt.map_image_to_costs('images/ex2.png')
distance_transform = wdt.get_weighted_distance_transform(cost_field)
wdt.plot(weighted_distane_transform)
```

## Authors

* **Omar Richardson**

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## References

 - Algorithms inspired by
    * [Tsitsiklis](http://www.mit.edu/~jnt/dijkstra.html)
    * [Treuille et al.](http://grail.cs.washington.edu/projects/crowd-flows/78-treuille.pdf)
 - Heap implementation based on Daniel Pena's [mheap](https://github.com/trifling/mheap)
 - Uses `f2py` to create Python libraries
