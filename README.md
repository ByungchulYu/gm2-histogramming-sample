# gm2-histogramming-sample
This repository provides a representative C++/ROOT histogramming package developed for sub-run–level data preparation in the Fermilab Muon g−2 experiment.
The code reads reconstructed ROOT TTrees, applies event selection, and constructs time-dependent decay spectra (“wiggle plots”) used for extracting the anomalous precession frequency ω_a.
## Structure
- HistogramBase.hh: abstract interface defining histogram booking and filling
- HistogramBase.hh: abstract interface defining histogram booking and filling
- runHistogramming.cc: main driver handling TTree I/O and event processing
- Makefile: minimal build configuration
This package illustrates my approach to scientific software development:
modular design, clear separation of interfaces and implementation, and reproducible data-processing workflows.
