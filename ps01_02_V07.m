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
kmax = 9;
vGridCapital = linspace(kmin,kmax,600);
% Get sizes of capital and prodctvt grid:
nGridCapital = length(vGridCapital);
nGridProductivity = length(vProductivity);
I = nGridCapital;
dk = (kmax-kmin)/(I-1);
kapMatrix = vGridCapital'*ones(1,nGridProductivity);

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
mc = zeros(nGridCapital,nGridProductivity);
%% 4. We pre-build output for each point in the grid
%profitMatrix = profit(exp(vProductivity),vGridCapital);
%% 5. Main iteration
maxDifference = 10.0;
maxite = 100000;
tolerance = 0.0000001;%0.0000001;
iteration = 0;
deltaV = 1000;
explicit = 2; % 1 -to run the explicit, 2- implicit
%Note: this matrix only works when there are three states in the economy
Achange = [zeros(I), mTransition(2,1)*speye(I) mTransition(3,1)*speye(I);...
    mTransition(1,2)*speye(I) zeros(I) mTransition(3,2)*speye(I);...
    mTransition(1,3)*speye(I) mTransition(2,3)*speye(I) zeros(I)];

while (maxDifference>tolerance || iteration > maxite) 
%     iteration
%     maxDifference
    %Not necessary to loop over Productivity:
    %for nProductivity = 1:nGridProductivity
        dVf = mdVf;
        dVb = mdVb;
        %prodctvt = vProductivity(nProductivity);
        V = mValueFunction;
        % forward difference
        dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/dk;
        %dVf(I) = (aalpha*prodctvt.*kmax.^(aalpha-1) - delta)/r; %state constraint, for stability
        dVf(I,:) = 1;
        % backward difference
        dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/dk;
        %dVb(1) = (aalpha*prodctvt.*kmin.^(aalpha-1) - delta)/r; %state constraint, for stability
        dVb(1,:) = 1; %state constraint, for stability
        %investment with forward difference
        investmentVectorF = ((dVf-1)/bb + delta).*kapMatrix;
        muf = investmentVectorF - delta*kapMatrix;
        %investment with backward difference
        investmentVectorb = ((dVb-1)/bb + delta).*kapMatrix;
        mub = investmentVectorb - delta*kapMatrix;
        %investment and derivative of value function at steady state
        investmentVector0 = delta*kapMatrix;
        dV0 = 1+bb*(investmentVector0./kapMatrix -delta);%-bb/2*delta^2 +1;   %%%%%%%%%% NOT SURE %%%%%%%%%%%%%%%%%%
        % dV_upwind makes a choice of forward or backward differences based on
        % the sign of the drift    
        If = muf > 0; %below steady state
        Ib = mub < 0; %above steady state
        I0 = (1-If-Ib); %at steady state
        dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0; %important to include third term
        invest = ((dV_Upwind-1)/bb + delta).*kapMatrix;
        %Get vectors and do implicit Euler...
%         diffPoisson = mValueFunction*mTransition(:,nProductivity);
        dn = vGridCapital'.^aalpha*vProductivity - invest -...
            bb/2*(invest./kapMatrix -delta).^2.*kapMatrix ;
        driftK = muf.*If + mub.*Ib + 0.*I0;
        %CONSTRUCT MATRIX
        X = -min(mub,0)/dk;
        Z = max(muf,0)/dk;
        % To build Y, we need to consider the diffPoisson effect
        Y = -max(muf,0)/dk + min(mub,0)/dk - ones(I,1)*sum(mTransition'-diag(mTransition'),2)';
        %sparse matrix: faster
        A1 = spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
        A2 = spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
        A3 = spdiags(Y(:,3),0,I,I)+spdiags(X(2:I,3),-1,I,I)+spdiags([0;Z(1:I-1,3)],1,I,I);
        AAA = [A1, sparse(I,I) sparse(I,I);...
            sparse(I,I) A2 sparse(I,I);...
            sparse(I,I) sparse(I,I) A3] + Achange;
        B = (r + 1/delta)*speye(3*I) - AAA;
        % We need to write down all in vector form, as we stacked all
        d_n = dn(:);
        V_n = V(:);
        vecb = d_n + V_n/deltaV;
        Vnew = B\vecb; %SOLVE SYSTEM OF EQUATIONS
        Vnew = reshape(Vnew, [I,3]); %rebuild the matrix
        
    mValueFunctionNew = Vnew;
    mdVf = dVf;
    mdVb = dVb;
    mPolicyFunction = invest;
%         pause
    %end
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
