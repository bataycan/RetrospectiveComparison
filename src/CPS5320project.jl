module RetrospectiveComparison
    import Distributions
    import Statistics
    import StatsBase
    import LinearAlgebra
    import Plots
    import SmoothingSplines
    import Pkg
    using Pkg
    Pkg.add("DataFrames")

    using DataFrames
    using Statistics
    using SmoothingSplines
    using RDatasets
    using StatsBase
    using Plots
    using LinearAlgebra
    using Distributions



    function movingaverage(X::Vector,m::Int)
        n = length(X)
        H = zeros(n, n)
        for i in 1:n
            for j in i:i+m
                if j <= n
                    H[i,j] = 1/(2*m+1)
                end
            end
            for j in i-m:i-1
                if j > 0
                    H[i,j] = 1/(2*m+1)
                end
            end
            for j in i+n-m:n
                if j <= n
                    H[i,j] = 1/(2*m+1)
                end
            end
            for j in i-n:i-n+m
                if j > 0
                    H[i,j] = 1/(2*m+1)
                end
            end
        end
        return H * X
    end

    function pascalTri(m::Int)
        # Creating Pascals triangle and using m to determine the row
        # needed (i.e. the coefficients for the weighted H matrix)
        r=[1, 2, 1]
        pr = [1, 2, 1]
        i = 2
        while i < m
            append!(r, 1)
            for j in 2:length(r)-1
                r[j] = pr[j] + pr[j-1]
            end
            i = i+1
            pr = copy(r)
        end
        r = r/sum(r)
        return r
    end

    function binomial(X::Vector, m::Int)
        n = length(X)
        H = zeros(n, n)
        r = pascalTri(m)

        for i in 1:n
            for j in i:i+floor(Int, m/2)
                if j <= n
                    H[i,j] = r[j-i+floor(Int, (length(r)+1)/2)]
                    H[j,i] = H[i,j]
                end
            end
            for j in i+n-floor(Int, m/2):n
                if j <= n
                    H[i,j] = r[i-j+floor(Int, (length(r)+1)/2)+n]
                    H[j,i] = H[i,j]
                end
            end
        end
        return H * X
    end

    # Inputting the X and Y, this function will determine the best suited m using the moving average algorithm
    function movingaverage_best_m(X::Vector, Y::Vector, mInitial::Int, mFinal::Int, mStep::Int)
        
        n = length(X)
        H = zeros(n, n)
        steps = ((mFinal - mInitial) / mStep) + 1
        
        steps = Int.(steps)
        mVec = zeros(steps)
        Ms = zeros(steps)
        #mVec = Vector{Float64}()
        index = 1
        #Creation of the error table without inserting any values
        errTable = DataFrame()
        for m in mInitial:mStep:mFinal
            
            Ypred = movingaverage(X, m)

            errVec = zeros(length(Y))
            for ind in 1:length(errVec)
                errVec[ind]  = Y[ind] - Ypred[ind]
            end

            errTable[!,:"m = " * string(m)] = errVec

            # test the error between the predicted and actual
            # Using the Mean Standard Deviation
            error = 0
            error = msd(Y, Ypred)

            #add error to mVec
            Ms[index] = m
            mVec[index] = error
            index += 1
        end
        # Locate where the smallest error is and return the optimal m
        loc = argmin(mVec)
        optimalM = mInitial + loc*(mStep) - 1*(mStep)
        display(plot(Ms, mVec, labels = "Error Trends for Moving Average"))
        return optimalM
    end

    function binomial_best_m(X::Vector, Y::Vector, mInitial::Int, mFinal::Int, mStep::Int)
        
        n = length(X)
        H = zeros(n, n)
        steps = ((mFinal - mInitial) / mStep) + 1
        
        steps = Int.(steps)
        mVec = zeros(steps)
        Ms = zeros(steps)
        #mVec = Vector{Float64}()
        errTable = DataFrame()
        index = 1
        for m in mInitial:mStep:mFinal
            
            Ypred = binomial(X, m)

            errVec = zeros(length(Y))
            for ind in 1:length(errVec)
                errVec[ind]  = Y[ind] - Ypred[ind]
            end

            errTable[!,:"m = " * string(m)] = errVec

            # test the error between the predicted and actual
            # Using the Mean Standard Deviation
            error = 0
            error = msd(Y, Ypred)

            #add error to mVec and index the x axis for the plot
            Ms[index] = m
            mVec[index] = error
            index += 1
        end
        # Locate where the smallest error is and return the optimal m
        loc = argmin(mVec)
        optimalM = mInitial + loc*(mStep) - 1*(mStep)
        display(plot(Ms, mVec, labels = "Error Trends for Binomial Filtering"))
        return optimalM
    end


    function plotComparison(X::Vector, Y::Vector, mInitial::Int, mFinal::Int, mStep::Int)
        #First finds the best m for both bionomial and moving average 
        bM = binomial_best_m(X, Y, mInitial, mFinal, mStep)
        mM = movingaverage_best_m(X, Y, mInitial, mFinal, mStep)

        #Utilizes the best m and rerun program to perform a plot comparison
        MA = movingaverage(X, mM)
        B = binomial(X, bM)
        
        plot(MA, labels = "Moving Average: m = " * string(mM))
        plot!(B, labels = "Binomial: m = " * string(bM))
        display(plot!(X, Y, seriestype = :scatter, label = "Actual", title = "Binomial vs. Moving Average"))
    end
end