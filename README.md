# SARA-adaptive-ppd
Matlab code associated with the method developed in  
> A. Dabbech, A. Onose, A. Abdulaziz, R. A. Perley, O. M. Smirnov, Y. Wiaux, [Cygnus A super-resolved via convex optimization from VLA data](https://doi.org/10.1093/mnras/sty372), *Mon. Not. Roy. Astron. Soc.*, 476(3):2853–2866, 2018.

The code also uses a set of MATLAB scripts, published
> A. Onose, A. Dabbech and Y. Wiaux, [An accelerated splitting algorithm for radio-interferometric imaging: when natural and uniform weighting meet](http://dx.doi.org/10.1093/mnras/stx755), *Mon. Not. Roy. Astron. Soc.*, 469(1):938-949, 2017.

The version of the codes released when the paper was published are available on the `v0` branch.

**Setup**: if needed, configure the paths in `minimalist_script.m` before running.

**Dependencies:** the present codes depend on the content of the `measurement-operator` github repository, loaded as a github `submodule`. This module contains codes associated with the following publications

> J. A. Fessler and B. P. Sutton, Nonuniform Fast Fourier Transforms Using Min-Max Interpolation, *IEEE Trans. Image Process.*, vol. 51, n. 2, pp. 560--574, Feb. 2003.
>
> A. Dabbech, L. Wolz, L. Pratley, J. D. McEwen and Y. Wiaux, [The w-effect in interferometric imaging: from a fast sparse measurement operator to superresolution](http://dx.doi.org/10.1093/mnras/stx1775), *Mon. Not. Roy. Astron. Soc.*, 471(4):4300-4313, 2017.

**Installation** To properly clone the project with the content of the submodule, you may need to do follow one of set of instructions:

- cloning the repository

```bash
git clone --recurse-submodules https://github.com/basp-group/SARA-adaptive-ppd.git
```

- updating from an existing `SARA-adaptive-ppd` repository

```bash
git pull
git submodule sync --recursive # update submodule address, in case the url has changed
git submodule update --init --recursive # update the content of the submodules
git submodule update --remote --merge # fetch and merge latest state of the submodule
```
