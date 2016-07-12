# Streaming Kernel Principal Component Analysis

Implementation of the approach presented in:

Mina Ghashami, Daniel J. Perry, and Jeff M. Phillips. "Streaming Kernel Principal Component Analysis." Proceedings of the 19th International Conference on Artificial Intelligence and Statistics. 2016.


## Installation

### Prerequisites
These three unregistered julia packages implement various prerequisites such as Frequent Directions (FD), Random Fourier Features (RFF), and the Gaussian (RBF) kernel.
They can be installed by launching the julia prompt and running the following commands:
> Pkg.clone("git@github.com:daniel-perry/Sketch.jl.git")
> Pkg.clone("git@github.com:daniel-perry/RandomFourierFeatures.jl.git")
> Pkg.clone("git@github.com:daniel-perry/Kernel.jl.git")

### Install StreamingKPCA
The streaming KPCA package itself can be installed by running this command at the julia prompt:
> Pkg.clone("git@github.com:daniel-perry/StreamingKPCA.jl.git")


