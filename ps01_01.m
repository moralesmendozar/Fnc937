%% Basic RBC model with full depreciation
%
% Jesus Fernandez-Villaverde
% Haverford, July 3, 2013

%% 0. Housekeeping

clear all
close all
clc

tic

%%  1. Calibration
clear
clc

rho = 0.8;   % ar(1) parameter of log(a) (productivity)
r = 0.02;    % 1/(1+r) discount rate of firm
bbeta = 1/(1+r);   % Discount factor of the firm
delta = 0.1; % depreciation
theta1 = 0.3; % kapital elasticity (cobbdouglas)
aalpha =theta1;     % Elasticity of output w.r.t. capital
theta2 = 0.6; % labor elasticity (cobbdouglas)
W = 2;   % wage
sigma = 0.1;  % std dev of eps (ar(1)) of log(a) (productivity)
abar = 1; % log(abar) is the mean of log(a), which is ar(1)
b0 = 1;
b1 = 0.5;

% Productivity values
na = 9;  %number of points for logaGrid = [loga0, loga1... loga9]
logabar = log(abar);
deltaLoga = 0.3;  %change between each log(a);

vProductivity = (logabar - deltaLoga*floor(na/2)):deltaLoga:(logabar + deltaLoga*floor(na/2));

% Transition matrix
P   = eye(na);
a_1 = vProductivity(1);
a_n = vProductivity(na);
for j = 1:na
    aj = vProductivity(j);
    upperBoundA = (a_1 - rho*aj - (1-rho)*abar +deltaLoga/2)/sigma;
    P(1,j) = normcdf(upperBoundA);
    lowerboundA = (a_n - rho*aj - (1-rho)*abar -deltaLoga/2)/sigma;
    P(na,j) = 1-normcdf(lowerboundA);
end
for i = 2:(na-1)
    for j = 1:(na)
        ai = vProductivity(i);
        aj = vProductivity(j);
        upperBoundA = (ai - rho*aj - (1-rho)*abar +deltaLoga/2)/sigma;
        lowerboundA = (ai - rho*aj - (1-rho)*abar -deltaLoga/2)/sigma;
        P(i,j) = normcdf(upperBoundA)-normcdf(lowerboundA);
    end
end
P'
mTransition   = P;


%% 2. Steady State

capitalSteadyState = (aalpha*bbeta)^(1/(1-aalpha));
outputSteadyState = capitalSteadyState^aalpha;
consumptionSteadyState = outputSteadyState-capitalSteadyState;

fprintf(' Output = %2.6f, Capital = %2.6f, Consumption = %2.6f\n', outputSteadyState, capitalSteadyState, consumptionSteadyState); 
fprintf('\n')

% We generate the grid of capital
capMiddle = (aalpha*abar/(r+delta))^(1/(1-aalpha));
kstep = 0.01; %0.00001
%vGridCapital = 0.5*capitalSteadyState:kstep:1.5*capitalSteadyState;
vGridCapital = 0.5*capMiddle:kstep:1.5*capMiddle;
%kapitalMax = (abar/(delta))^(1/(1-aalpha));
%vGridCapital = 0:kstep:kapitalMax;


nGridCapital = length(vGridCapital);
nGridProductivity = length(vProductivity);

%% Plot of the profit function
labor = @(a,k) (theta1*(k.^theta1)'*a/W).^(1/(1-theta2));
profit = @(a,k) ((k'.^theta1)*a).*(labor(a,k).^theta2) - W*labor(a,k);
investment = @(a,k,kprime) kprime - (1-delta)* k'*ones(size(a,1));
phi = @(a,k,kprime) b0 * k'*ones(size(a,1)) + b1*(( investment(a,k,kprime)./ ( k'*ones(size(a,1)) )- delta ).^2).*( k'*ones(size(a,1)) ) ;

plot(vGridCapital, profit(abar,vGridCapital));

%% 3. Required matrices and vectors

mOutput           = zeros(nGridCapital,nGridProductivity);
mValueFunction    = zeros(nGridCapital,nGridProductivity);
mValueFunctionNew = zeros(nGridCapital,nGridProductivity);
mPolicyFunction   = zeros(nGridCapital,nGridProductivity);
expectedValueFunction = zeros(nGridCapital,nGridProductivity);

%% 4. We pre-build output for each point in the grid

profitMatrix = profit(exp(vProductivity),vGridCapital);

%% 5. Main iteration

maxDifference = 10.0;
tolerance = 0.0000001;
iteration = 0;

while (maxDifference>tolerance)  
    
    expectedValueFunction = mValueFunction*mTransition';
    
    for nProductivity = 1:nGridProductivity
        
        % We start from previous choice (monotonicity of policy function)
        gridCapitalNextPeriod = 1;
        
        for nCapital = 1:nGridCapital
                        
            valueHighSoFar = -1000.0;
            capitalChoice  = vGridCapital(1);
            
            for nCapitalNextPeriod = gridCapitalNextPeriod:nGridCapital
                
                %consumption = mOutput(nCapital,nProductivity)-vGridCapital(nCapitalNextPeriod);
                prodctvt = exp(vProductivity(nProductivity));
                currentCapital = vGridCapital(nCapital);
                kprimeTomorow = vGridCapital(nCapitalNextPeriod);
                dividend = profitMatrix(nCapital,nProductivity) - investment(prodctvt,currentCapital,kprimeTomorow) - phi(prodctvt,currentCapital,kprimeTomorow);
                valueProvisional = (1-bbeta)*dividend+bbeta*expectedValueFunction(nCapitalNextPeriod,nProductivity);
            
                if (valueProvisional>valueHighSoFar)
                    valueHighSoFar = valueProvisional;
                    capitalChoice = vGridCapital(nCapitalNextPeriod);
                    gridCapitalNextPeriod = nCapitalNextPeriod;
                else
                    break; % We break when we have achieved the max
                end    
                  
            end
            
            mValueFunctionNew(nCapital,nProductivity) = valueHighSoFar;
            mPolicyFunction(nCapital,nProductivity) = capitalChoice;
            
        end
        
    end
    
    maxDifference = max(max(abs(mValueFunctionNew-mValueFunction)));
    mValueFunction = mValueFunctionNew;
    
    iteration = iteration+1;
    if (mod(iteration,10)==0 || iteration ==1)
        fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference); 
    end
           
end

fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference); 
fprintf('\n')
%fprintf(' My check = %2.6f\n', mPolicyFunction(1000,3));
fprintf('\n')

toc

%% 6. Plotting results

figure(1)

subplot(3,1,1)
plot(vGridCapital,mValueFunction)
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
xlabel('k')
title('Value Function')

subplot(3,1,2)
plot(vGridCapital,mPolicyFunction)
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
xlabel('k')
title('Policy Function')

vExactPolicyFunction = aalpha*bbeta.*(vGridCapital.^aalpha);

subplot(3,1,3)
plot((100.*(vExactPolicyFunction'-mPolicyFunction(:,3))./mPolicyFunction(:,3)))
title('Comparison of Exact and Approximated Policy Function')

%set(gcf,'PaperOrientation','landscape','PaperPosition',[-0.9 -0.5 12.75 9])
%print('-dpdf','Figure1.pdf')