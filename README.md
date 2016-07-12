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

## Running the example

After installing the above packages, you can run the example script to compare between SKPCA, RNCA, and Nystrom similar to those done in the paper.

To run the examples run,

> julia comparison_exampe.jl

this will generate a results .mat file.

To plot the results run,

> julia comparison_plot.jl results.mat gaussian

(here gaussian refers to the dataset name which is set to "gaussian" by default).
The second command will plot the results for comparison.
