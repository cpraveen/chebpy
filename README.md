# chebpy: Spectral and Chebyshev methods in Python

See this as a jupyter-book here

https://cpraveen.github.io/chebpy

These examples are based on the book 

> L. N. Trefethen, Spectral Methods in MATLAB, SIAM.

I have rewritten these codes in Python.

For nbviewer links to these files, see here

http://cpraveen.github.io/teaching/chebpy.html

Run the code in binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cpraveen/chebpy/HEAD)

Run the code in colab: [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/cpraveen/chebpy)

Matlab versions of the code are available here

https://github.com/cpraveen/spectral_matlab

## Building jupyter-book

Install 

```shell
conda install jupyter-book ghp-import
```

Build html and publish

```shell
git checkout git@github.com:cpraveen/chebpy
jb build chebpy
cd chebpy
ghp-import -n -p -f _build/html
```
