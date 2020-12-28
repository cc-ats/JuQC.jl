abstract type ObjectiveFunction end 

abstract type OptimizationAlgorithm end

struct DIISOptimizer
    tol::Number
end