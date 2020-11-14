# SARA-adaptive-ppd

This code depends on ./irt (Fessler 2003) and ./sara-ppd (Onose et al. 2017).
Make sure to change their associated  paths in minimalist_script.m

**Dependencies:** the present codes includes a slightly modified version of the MATLAB NUFFT algorithm available at http://web.eecs.umich.edu/~fessler/irt/fessler.tgz, described in

> J. A. Fessler and B. P. Sutton - 
<strong>Nonuniform Fast Fourier Transforms Using Min-Max Interpolation</strong>, <em>IEEE Trans. Image Process.</em>, vol. 51, n. 2, pp. 560--574, Feb. 2003.

**Installation** To properly clone the project with the submodules, you may need to do follow one of set of instructions:

- updating from an existing `SARA-adaptive-ppd` repository:

```bash
git pull
git submodule sync --recursive # update submodule address, in case the url has changed
git submodule update --init --recursive # update the content of the submodules
```

- cloning the repository from scratch

```bash
git clone --recurse-submodules https://github.com/basp-group-private/SARA-adaptive-ppd.git
```
