# How To Pick a Peak
This repository contains all the essential code necessary to reproduce the results presented in Dahlbom and Braasch (2020),
including a complete implementation of the pitch model presented there. Precise reproduction of the figures requires trial data
that are not included in this repository (they are too large), but I am happy to provide them on request.

The code here can be taken as a general framework for matching "adaptive templates" to one-dimensional data, using mathematical techniques similar to those
presented in the work of William Sethares (see, for example, his [book](https://www.springer.com/gp/book/9781846286391) on rhythm).  

The pitch model itself can be run by calling the function `estimate_τ0` from the PitchModelPeaks.jl module.  This is a rather
monstrous function, which includes many options for making comparisons with other models and generating graphs, but
it can be called in a straightforward way simply by providing a signal and a sample rate.

The paper also contains an appendix providing the results of the simulation of a bushy
cell in the cochlear nucleus recieving a nearly-harmonic input signal.  This material is contained in a [separate repository](https://github.com/analogouscircuit/rinzelmodel).

## Pitch Model
The pitch model proposed in the paper carries over all the basic features of the standard autocorrelation model (Meddis and Hewitt, 1991), incorporates a recent extension that permits the model to predict noise edge pitch (Hartmann, Cariani and Colburn, 2019), and further refines the peak selection process so that the model can also predict the pitch shifts due to mistuned harmonics. (Image from Dahlbom and Braasch (2020).)

![Model in action](/images/Figure5.png)

## References
Dahlbom, D. A. and Braasch, J. (2020).  "How to pick a peak: Pitch and peak shifting in temporal models of pitch perception," 
Journal of the Acoustical Society of America 147, 2713-2727.

Hartmann, W. M., Cariani, P. A., and Colburn, S. (2019). "Noise edge pitch and models of pitch perception," Journal of the Acoustical Society of America 145, 1993-2008.

Meddis, R. and Hewitt, M. J. (1991). "Virtual pitch and phase sensitivity of a computer model of the auditory periphery. I: Pitch identification," Journal of the Acoustical Society of America 89, 2866-2882.
