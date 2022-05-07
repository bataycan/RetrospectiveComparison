# Retrospective Comparison

A Julia package that utilizes two prediction algorithms, Moving Average and Binomial, and allows the capabilities to make a comparison using there best finding parameter to do so. To briefly describe these two methods, both Moving Average and Binomial use a smoothing matrix to make predicts given an inputted vector. The generated smoothing matrix is where is these two algorithms differ, the Moving Average algorithm uses even weights for each step, where as Binomial Filtering uses the Pascal Triangle to determine the weights. 

This Package is completely independent of any existing Moving Average and Binomial Filtering packages previously developed in Julia, or any other programming language. Which is why two of the function that have been created are the algorithms for Moving Average and Binomial Filtering, both with output that return the predicted Y vector with the given inputted vector and tested m parameter. With this idea, the best_m algorithm was derived, instead of choose a single parameter the best_m algorithm use an interval of this parameter and test them individually. Additionally to the inputted interval, a step size parameter is also inputted to determine how to move throughout this interval. Although the actual values is typically unknown when working with real data, when making a comparison of these two algorithms since they are both data driven methods it is required to use the actual data to be compared against the predictions. 

```julia
using Retrospective Comparison
using Random
using Distributions

x=rand(Uniform(0.0,30.0), 30)
y = zeros(30)
for i in 1:length(x)
    y[i] = 25*cos(pi * x[i])
end

plotComparison(x, y, 2, 11, 3)
bM, test, ttable = binomial_best_m(x, y, 1, 6, 1)
mM, test, ttable = movingaverage_best_m(x, y, 1, 6, 1)

```

### TODO
* Parallelize the algorithm
* Enhance the algorithm to be more efficient
* Create restraints on the parameters limitlessness
* incorprate the cross validation function that has been included, into more of the other comparison function. 


**References**
Bezanson, J., Karpinski, S., Shah, V. B., & Edelman, A. (2012). Julia: A fast dynamic language for technical computing. ArXiv Preprint ArXiv:1209.5145.

Chatla, S. (2022). Smoothing. CPS 5320 Lecture Notes. El Paso ; University of Texas at El Paso. 
