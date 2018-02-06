![logo](https://upload.wikimedia.org/wikipedia/commons/7/77/Trident_logo.svg) 

# HiGrad

This is the R package for HiGrad (Hierarchical Incremental Gradient Descent), 
an algorithm for statistical inference for online learning and stochastic approximation.
Also included in this repository are the R codes for reproducing the results and plots
of the simulation and real data study in Su & Zhu, 2018.

#### Description


Stochastic gradient descent (SGD) is an immensely popular approach for online learning in 
settings where data arrives in a stream or data sizes are very large. 
However, despite an ever-increasing volume of work on SGD, much less is known about the 
statistical inferential properties of SGD-based predictions. 
Taking a fully inferential viewpoint, this paper introduces a novel procedure termed 
HiGrad to conduct statistical inference for online learning, 
without incurring additional computational cost compared with the vanilla SGD. 
The HiGrad procedure begins by performing SGD iterations for a while and then split the single thread into a few, 
and this procedure hierarchically operates in this fashion along each thread. 
With predictions provided by multiple threads in place, 
a t-based confidence interval is constructed by decorrelating predictions 
using covariance structures given by the Ruppert–Polyak averaging scheme. 
Under certain regularity conditions, the HiGrad confidence interval 
is shown to attain asymptotically exact coverage probability. 

#### Installation

Follow the code below to install the package. 
You will need to install the [devtools](https://cran.r-project.org/package=devtools) package first in case you haven't. 

    devtools::install_github(repo = "captainyc/higrad", subdir="higrad")

#### Reference

Weijie Su and Yuancheng Zhu. (2018) *Statistical Inference for Online Learning and Stochastic Approximation via Hierarchical Incremental Gradient Descent*.
