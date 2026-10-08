Overview
This repository contains three R scripts related to the estimation and inference of multivariate fractional Brownian motion (mfBm). The codes are based on the methodology proposed in:

Markus Bibinger, Jun Yu, and Chen Zhang, Modeling and Forecasting Realized Volatility with Multivariate Fractional Brownian Motion.

The three files included are:

setup_code.R – Preliminary functions and definitions adapted from existing literature.

BYZestimator.R – Estimation procedures developed in our paper.

BYZinference.R – Inference procedures developed in our paper.

Setup
Before running either MYZestimator.R or MYZinference.R, please make sure to run setup_code.R to load the necessary functions and definitions.

The setup code is adapted from:

Amblard, P.-O., Coeurjolly, J.-F., Philippe, A., and Lavancier, F. (2012), and related work.

Notes
The estimation and inference codes are developed specifically for the models introduced in our study.

The inference code can be executed independently, once the setup code is loaded.

Please cite our work if you use these scripts in your research.