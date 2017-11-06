# Weighted distance transform
A Python implementation for computing a weighted distance transform of an image.

A [distance transform](https://en.wikipedia.org/wiki/Distance_transform) is a map of an image that assign to each pixel its distance to the nearest boundary.
A weighted distance transform extends this by allowing for ... weighted distances. Useful for image analysis or computing potential functions in path finding.

This script implements a fast marching algorithm in Python with some acceleration in Fortran (if a Fortran compiler is present).

### Prerequisites

 * Python 3 (but no real obstacles for Python 2)
 * Numpy/Scipy and the python imaging library (PIL)
 * Matplotlib for plotting images (not required)
 * Fortran compiler, makes script significantly faster (recommended, not required)

### How does it work?

![Image with text](/images/example.png?raw=true "Example image")
Input (left): A PNG/JPG image with with boundaries/exits from which to compute the distance (red), obstacles (full black), normal space (white), less 'accessible space' (grey-ish).

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

The original use case of this script is creating motion planners for 2D maps (in [crowd dynamics](https://symbols.hotell.kau.se/2016/11/30/mercurial/)). We identify goals with _red_, obstacles in _black_, and accessible space in _white_. 
You can use any color (except full red) to indicate an area that is less easily accessible (and preferably should be avoided). The darker the colour (mapped to a greyscale), the more difficult it is to move through. 

Once you created an image, you can provide it to the script as shown below. 
Otherwise, run the script with some example figures (like `images/ex1.png`)

```python
wdt = WDT('images/ex1.png')
wdt.plot()
```

## Authors

* **Omar Richardson**

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## References

Algorithms inspired by
* [Tsitsiklis](http://www.mit.edu/~jnt/dijkstra.html)
* [Treuille et al.](http://grail.cs.washington.edu/projects/crowd-flows/78-treuille.pdf)

