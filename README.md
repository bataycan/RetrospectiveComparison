# Retrospective Comparison

A Julia package that utilizes two prediction algorithms, Moving Average and Binomial, and allows the capabilities to make a comparison using there best finding parameter to do so. To briefly describe these two methods, both Moving Average and Binomial use a smoothing matrix to make predicts given an inputted vector. The generated smoothing matrix is where is these two algorithms differ, the Moving Average algorithm uses even weights for each step, where as Binomial Filtering uses the Pascal Triangle to determine the weights. 

This Package is completely independent of any existing Moving Average and Binomial Filtering packages previously developed in Julia, or any other programming language. Which is why two of the function that have been created are the algorithms for Moving Average and Binomial Filtering, both with output that return the predicted Y vector with the given inputted vector and tested m parameter. With this idea, the best_m algorithm was derived, instead of choose a single parameter the best_m algorithm use an interval of this parameter and test them individually. Additionally to the inputted interval, a step size parameter is also inputted to determine how to move throughout this interval. Although it is typically unknown when working with real data, for this function 



### TODO
* automatic selection of m parater
* subsampling of design grid for higher efficiency 


**References**
