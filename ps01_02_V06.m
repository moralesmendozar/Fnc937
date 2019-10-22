%% UPENN, @Wharton
% Finance 937. 
% Prof. Joao Gomes
% Student: Rodrigo A. Morales M. && Mr. Paw Bednarek
% Based on Moll HJB code
% Okt, 2019
% Problem Set 01. Exercise 2) continuous time model
%% 0. Housekeeping
clear all
close all
clc
tic
%%  1. Calibration
aalpha = 0.3;   % Elasticity of output w.r.t. capital
bb = 0.5;       % cost of investment
delta = 0.05;  % start with a very small one
r = 0.05;       % interest rate, return of capital
% Productivity values
vProductivity = [0.9 1.0 1.1];
% vProductivity = 1;
nGridProductivity = length(vProductivity);
na = nGridProductivity;  %number of points of Productivyt
if na >1
    mTransition = [1/2 1/2 0; 1/4 1/2 1/4; 0 1/2 1/2]';
else
    mTransition = 1;
end

%% 2. Steady State
%define the functions of labor, profit and investment...
abar = 1.0;
capitalSteadyState = (aalpha*1/(r+delta))^(1/(1-aalpha));
outputSteadyState = abar*capitalSteadyState^aalpha;
%fprintf(' Output = %2.6f, Capital = %2.6f\n', outputSteadyState, capitalSteadyState); 
fprintf('Code Began running...\n')
% Generate the grid of capital
kapitalMax = (abar/(delta))^(1/(1-aalpha));
kmin = 0.01;
kmax = 10;
vGridCapital = linspace(kmin,kmax,200);
% Get sizes of capital and prodctvt grid:
nGridCapital = length(vGridCapital);
nGridProductivity = length(vProductivity);
I = nGridCapital;
dk = (kmax-kmin)/(I-1);

%% 3. Required matrices and vectors
%mValueFunction    = linspace(1,vGridCapital(nGridCapital),nGridCapital)'*vProductivity;
mValueFunction    = (vGridCapital'.^aalpha)*vProductivity;
%mValueFunction(1,:) = 6*ones(1,nGridProductivity);
%mValueFunction(end,:) = (vGridCapital(end)^aalpha./r)*vProductivity;
%mValueFunctionNew = zeros(nGridCapital,nGridProductivity);
mValueFunctionNew = mValueFunction;
mPolicyFunction   = zeros(nGridCapital,nGridProductivity); %investment
% Matrices for the Implicit iteration in each productivity:
mdVf = zeros(nGridCapital,nGridProductivity);
mdVb = zeros(nGridCapital,nGridProductivity);
%% 4. We pre-build output for each point in the grid
%profitMatrix = profit(exp(vProductivity),vGridCapital);
%% 5. Main iteration
maxDifference = 10.0;
maxite = 100000;
tolerance = 0.00099;%0.0000001;
iteration = 0;
deltaV = 0.00001;
explicit = 1; % 1 -to run the explicit, 2- implicit

while (maxDifference>tolerance || iteration > maxite) 
%     iteration
%     maxDifference
    for nProductivity = 1:nGridProductivity
        dVf = mdVf(:,nProductivity);
        dVb = mdVb(:,nProductivity);
        prodctvt = vProductivity(nProductivity);
        V = mValueFunction(:,nProductivity);
        % forward difference
        dVf(1:I-1) = (V(2:I)-V(1:I-1))/dk;
        %dVf(I) = (aalpha*prodctvt.*kmax.^(aalpha-1) - delta)/r; %state constraint, for stability
        dVf(I) = 1;
        % backward difference
        dVb(2:I) = (V(2:I)-V(1:I-1))/dk;
        %dVb(1) = (aalpha*prodctvt.*kmin.^(aalpha-1) - delta)/r; %state constraint, for stability
        dVb(1) = 1; %state constraint, for stability
        %investment with forward difference
        investmentVectorF = ((dVf-1)/bb + delta).*vGridCapital';
        muf = investmentVectorF - delta*vGridCapital';
        %investment with backward difference
        investmentVectorb = ((dVb-1)/bb + delta).*vGridCapital';
        mub = investmentVectorb - delta*vGridCapital';
        %investment and derivative of value function at steady state
        investmentVector0 = delta*vGridCapital';
        dV0 = 1+bb*(investmentVector0./vGridCapital' -delta);%-bb/2*delta^2 +1;   %%%%%%%%%% NOT SURE %%%%%%%%%%%%%%%%%%
        % dV_upwind makes a choice of forward or backward differences based on
        % the sign of the drift    
        If = muf > 0; %below steady state
        Ib = mub < 0; %above steady state
        I0 = (1-If-Ib); %at steady state
        dV_Upwind = dVf.*If + dVb.*Ib + 0.*I0; %important to include third term
        invest = ((dV_Upwind-1)/bb + delta).*vGridCapital';
        %Get vectors and do implicit Euler...
        diffPoisson = mValueFunction*mTransition(:,nProductivity);
        dn = prodctvt*vGridCapital'.^aalpha - invest -...
            bb/2*(invest./vGridCapital' -delta).^2.*vGridCapital' + ...
            diffPoisson;
        %CONSTRUCT MATRIX
        X = -min(mub,0)/dk;
        Y = -max(muf,0)/dk + min(mub,0)/dk;
        Z = max(muf,0)/dk;
        %sparse matrix: faster
        AAA =spdiags(Y,0,I,I)+spdiags(X(2:I),-1,I,I)+spdiags([0;Z(1:I-1)],1,I,I);
        B = (r + 1/delta)*speye(I) - AAA + speye(I);
        if explicit ==1 %Explicit Euler
            Vnew = deltaV*(dn + AAA*V + V - r*V)+V;
        else
            vecb = dn + 1/deltaV*V;
            Vnew = B\vecb; %SOLVE SYSTEM OF EQUATIONS
        end
        mValueFunctionNew(:,nProductivity) = Vnew;
        mdVf(:,nProductivity) = dVf;
        mdVb(:,nProductivity) = dVb;
        mPolicyFunction(:,nProductivity) = invest;
%         pause
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

figure(2)
endPlotN = nGridCapital;
subplot(3,1,1)
plot(vGridCapital(1:endPlotN),mValueFunction(1:endPlotN,:))
xlim([vGridCapital(1) vGridCapital(endPlotN)])
xlabel('k')
title('Value Function, continuous time')

subplot(3,1,2)
plot(vGridCapital(1:endPlotN),mPolicyFunction(1:endPlotN,:))
xlim([vGridCapital(1) vGridCapital(endPlotN)])
ylim([min(mPolicyFunction(:,1)) (max(mPolicyFunction(:,nGridProductivity))+1)])
xlabel('k')
title('Policy Function = Investment (continuous time)')

% vExactPolicyFunction = aalpha*bbeta.*(vGridCapital.^aalpha);
% 
subplot(3,1,3)
plot(vGridCapital(1:endPlotN),mPolicyFunction(1:endPlotN,:)-delta*vGridCapital(1:endPlotN)')
hold on
plot(vGridCapital(1:endPlotN),zeros(endPlotN,1),'-k')
xlim([vGridCapital(1) vGridCapital(endPlotN)])
ylim([min(mPolicyFunction(:,1)-delta*vGridCapital(1:endPlotN)') (max(mPolicyFunction(:,nGridProductivity)-delta*vGridCapital(1:endPlotN)')+2)])
xlabel('i(k) - delta*k')
title(' dk/dt  = I(k) - delta*k')
hold off

%set(gcf,'PaperOrientation','landscape','PaperPosition',...
%[-0.9 -0.5 12.75 9])
versionp = 'v06_implicit';
if na == 1
    name =  ['Plots/PS01_02_a_',versionp];
    name2 =  ['Plots/PS01_02_a_',versionp,'.fig'];
else
    name =  ['Plots/PS01_02_b_',versionp];
    name2 =  ['Plots/PS01_02_a_',versionp,'.fig'];
end
print('-djpeg',name)
