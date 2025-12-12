# gm2-histogramming-sample

This repository provides a representative C++/ROOT histogramming package developed
for sub-run–level data preparation in the Fermilab Muon g−2 experiment.
The code reads reconstructed ROOT TTrees, applies event selection, and constructs
time-dependent decay spectra (“wiggle plots”) used for extracting the anomalous
precession frequency ω_a.

## Structure

- `HistogramBase.hh`  
  Abstract interface defining histogram booking and filling.

- `Byu2Histograms.hh / .cc`  
  Experiment-specific histogram implementations derived from the base interface.

- `runHistogramming.cc`  
  Main driver handling TTree I/O, event selection, and histogram filling.

- `Makefile`  
  Minimal build configuration for compiling the package with ROOT.

This package illustrates my approach to scientific software development:
modular C++ design, clear separation of interfaces and implementation, and
reproducible data-processing workflows suitable for large experimental datasets.
