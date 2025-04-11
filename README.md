# CGM Signal Denoising

This project performs deconvolution of Continuous Glucose Monitoring (CGM) data to estimate the underlying glucose signal and its first derivative. The goal is to recover the glucose rate of change by applying regularized deconvolution techniques on noisy CGM measurements.

## Project Description

The analysis is based on the following model:
```matlab
y(t) = (h * x)(t) + v(t)
```

Where:
- `y(t)` is the measured CGM signal
- `x(t)` is the true (unknown) glucose input signal
- `h(t)` is the system’s impulse response
- `v(t)` is the measurement noise

The main steps of the analysis are:
1. **Convolution Model Setup**  
   The CGM signal is modeled as a convolution of the input glucose signal with an exponential impulse response.
2. **Regularized Deconvolution**  
   The problem is solved using gradient-based optimization with Tikhonov regularization.
3. **Derivative Estimation**  
   The deconvolved signal is used to estimate both glucose levels and their rate of change over time.

## Files
- main_script.m: Main script performing the full analysis.


## Requirements
- MATLAB
- Basic signal processing toolbox
