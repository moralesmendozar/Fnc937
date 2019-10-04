## Basic RBC model with full depreciation
#
# Jesus Fernandez-Villaverde
# Haverford, July 29, 2013

# start julia with --depwarn=no
# Redefine log
#mylog(x::Float64)  = ccall((:log, :libm), Float64, (Float64,), x)
#include("/home/moralesmendozar/Dropbox/03_UPenn/classes/2019_fall/01_937-/PS01/ps01_01_Fast.jl")
mylog(x::Float64)  = ccall((:log, :libopenlibm), Float64, (Float64,), x)

function main()

    ##  1. Calibration

    α = 1/3     # Elasticity of output w.r.t. capital

    ρ = 0.8;   % ar(1) parameter of log(a) (productivity)
    r = 0.02;    % 1/(1+r) discount rate of firm
    β = 1/(1+r);   % Discount factor of the firm
    δ = 0.1; % depreciation
    θ1 = 0.3; % kapital elasticity (cobbdouglas)
    θ2 = 0.6; % labor elasticity (cobbdouglas)
    W =2;   % wage
    σ = 0.1;  % std dev of eps (ar(1)) of log(a) (productivity)
    abar = 1; % log(abar) is the mean of log(a), which is ar(1)
    b0 = 1;
    b1 = 0.5;

    # Productivity values
    vProductivity = [0.9792 0.9896 1.0000 1.0106 1.0212]

    # Transition matrix
    mTransition   = [0.9727 0.0273 0.0000 0.0000 0.0000;
                     0.0041 0.9806 0.0153 0.0000 0.0000;
                     0.0000 0.0082 0.9837 0.0082 0.0000;
                     0.0000 0.0000 0.0153 0.9806 0.0041;
                     0.0000 0.0000 0.0000 0.0273 0.9727]

    # 2. Steady State

    capitalSteadyState = (α*β)^(1/(1-α))
    outputSteadyState = capitalSteadyState^α
    consumptionSteadyState = outputSteadyState-capitalSteadyState

    println("Output = ",outputSteadyState," Capital = ",capitalSteadyState," Consumption = ",consumptionSteadyState)

    # We generate the grid of capital
    vGridCapital = collect(0.5*capitalSteadyState:0.000005:1.5*capitalSteadyState)

    nGridCapital = length(vGridCapital)
    nGridProductivity = length(vProductivity)

    # 3. Required matrices and vectors

    mOutput           = zeros(nGridCapital,nGridProductivity)
    mValueFunction    = zeros(nGridCapital,nGridProductivity)
    mValueFunctionNew = zeros(nGridCapital,nGridProductivity)
    mPolicyFunction   = zeros(nGridCapital,nGridProductivity)
    expectedValueFunction = zeros(nGridCapital,nGridProductivity)

    # 4. We pre-build output for each point in the grid

    mOutput = (vGridCapital.^α)*vProductivity;

    # 5. Main iteration

    maxDifference = 10.0
    tolerance = 0.00000001
    iteration = 0

    while(maxDifference > tolerance)
        expectedValueFunction = mValueFunction*mTransition';

        for nProductivity in 1:nGridProductivity

            # We start from previous choice (monotonicity of policy function)
            gridCapitalNextPeriod = 1

            for nCapital in 1:nGridCapital

                valueHighSoFar = -1000.0
                capitalChoice  = vGridCapital[1]

                for nCapitalNextPeriod in gridCapitalNextPeriod:nGridCapital

                    consumption = mOutput[nCapital,nProductivity]-vGridCapital[nCapitalNextPeriod]
                    @inbounds valueProvisional = (1-β)*mylog(consumption)+β*expectedValueFunction[nCapitalNextPeriod,nProductivity]

                    if (valueProvisional>valueHighSoFar)
                	       valueHighSoFar = valueProvisional
                	       capitalChoice = vGridCapital[nCapitalNextPeriod]
                	       gridCapitalNextPeriod = nCapitalNextPeriod
                    else
                	       break # We break when we have achieved the max
                    end

                end

                mValueFunctionNew[nCapital,nProductivity] = valueHighSoFar
                mPolicyFunction[nCapital,nProductivity] = capitalChoice

            end

        end

        maxDifference     = maximum(abs.(mValueFunctionNew-mValueFunction))
        mValueFunction, mValueFunctionNew = mValueFunctionNew, mValueFunction

        iteration = iteration+1
        if mod(iteration,20)==0 || iteration == 1
            println(" Iteration = ", iteration, " Sup Diff = ", maxDifference)
        end

    end

    println(" Iteration = ", iteration, " Sup Diff = ", maxDifference)
    println(" ")
    println(" My check = ", mPolicyFunction[1000,3])
    println(" My check = ", mValueFunction[1000,3])

end
