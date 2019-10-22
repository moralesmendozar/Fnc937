%% UPENN, @Wharton
% Finance 937. 
% Prof. Joao Gomes
% Student: Rodrigo A. Morales M. && Mr. Paw Bednarek
% 
% Based on Jesus Fernandez-Villaverde RBC comparison code
% Okt, 2019
% Problem Set 01. Exercise 2) continuous time model

%% 0. Housekeeping
clear all
close all
clc
tic
%%  1. Calibration
clear
clc
aalpha = 0.3;   % Elasticity of output w.r.t. capital
bb = 0.5;       % cost of investment
delta = 0.05;  % start with a very small one
r = 0.05;       % interest rate, return of capital
options = optimset('Display', 'off');  %when solves for kss
% Productivity values
vProductivity = [0.9 1.0 1.1];
vProductivity = 1;
nGridProductivity = length(vProductivity);
na = nGridProductivity;  %number of points for logaGrid = [loga0, loga1... loga9]
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
fprintf(' Output = %2.6f, Capital = %2.6f\n', outputSteadyState, capitalSteadyState); 
fprintf('\n')

% with kss, generate the grid of capital
capMiddle = (aalpha*1/(r+delta))^(1/(1-aalpha));
kstep = 0.1; %0.00001
%vGridCapital = 0.5*capitalSteadyState:kstep:1.5*capitalSteadyState;
vGridCapital = 0.5*capMiddle:kstep:1.5*capMiddle;
kapitalMax = (abar/(delta))^(1/(1-aalpha));
vGridCapital = (0):kstep:(0.6*kapitalMax);
%vGridCapital = (0):kstep:(capitalSteadyState);
%vGridCapital = linspace(0.01*capitalSteadyState,capitalSteadyState*1.5,200);
%vGridCapital = linspace(0.01*capitalSteadyState,7,200);
% vGridCapital = linspace(0.01*capitalSteadyState,6,200);
kmin = 0.01;
kmax = 7;
vGridCapital = linspace(kmin,kmax,200);
% Get sizes of capital and prod grid:
nGridCapital = length(vGridCapital);
nGridProductivity = length(vProductivity);
I = nGridCapital;
dk = (kmax-kmin)/(I-1);

%% 3. Required matrices and vectors
% mValueFunction    = zeros(nGridCapital,nGridProductivity);
%mValueFunction    = linspace(1,vGridCapital(nGridCapital),nGridCapital)'*vProductivity;
% mValueFunction    = ones(nGridCapital,nGridProductivity);
mValueFunction    = (vGridCapital'.^aalpha./r)*vProductivity;
%mValueFunction(1,:) = 6*ones(1,nGridProductivity);
%mValueFunction(end,:) = (vGridCapital(end)^aalpha./r)*vProductivity;
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
tolerance = 0.0001;%0.0000001;
iteration = 0;
deltaV = 0.0001;
%explicit = 1; % 1 -to run the explicit, 2- implicit

while (maxDifference>tolerance) 
%     iteration
%     maxDifference
    for nProductivity = 1:nGridProductivity
        dVf = mdVf(:,nProductivity);
        dVb = mdVb(:,nProductivity);
        c = mPolicyFunction(:,nProductivity); %investment
        prodctvt = vProductivity(nProductivity);
        V = mValueFunctionNew(:,nProductivity);
        % forward difference
        dVf(1:I-1) = (V(2:I)-V(1:I-1))/dk;
        dVf(I) = (aalpha*prodctvt.*kmax.^(aalpha-1) - delta)/r; %state constraint, for stability
        dVf(I) = dVf(I-1);
        % backward difference
        dVb(2:I) = (V(2:I)-V(1:I-1))/dk;
        dVb(1) = (aalpha*prodctvt.*kmin.^(aalpha-1) - delta)/r; %state constraint, for stability
        %investment with forward difference
        cf = ((dVf-1)/bb + delta).*vGridCapital';
        investmentVectorF = cf;
        muf = investmentVectorF - delta*vGridCapital';
        %investment with backward difference
        cb = ((dVb-1)/bb + delta).*vGridCapital';
        investmentVectorb = cb;
        mub = investmentVectorb - delta*vGridCapital';
        %investment and derivative of value function at steady state
        c0 = delta*vGridCapital';
        %prodctvt.*vGridCapital.^aalpha - delta.*vGridCapital;
        investmentVector0 = c0;
        dV0 = -bb/2*delta^2 +1;   %%%%%%%%%% NOT SURE %%%%%%%%%%%%%%%%%%
        % dV_upwind makes a choice of forward or backward differences based on
        % the sign of the drift    
        If = muf > 0; %below steady state
        Ib = mub < 0; %above steady state
        I0 = (1-If-Ib); %at steady state
        dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0; %important to include third term
        invest = ((dV_Upwind-1)/bb + delta).*vGridCapital';
        %Get vectors and do implicit Euler...
        diffPoisson = mValueFunction*mTransition(:,nProductivity);
        dn = prodctvt*vGridCapital'.^aalpha - invest -...
                    bb/2*(invest./vGridCapital' -delta).^2.*vGridCapital' + ...
                    diffPoisson;
        
        vectorDf = prodctvt*vGridCapital'.^aalpha - cf -...
                    bb/2*(cf./vGridCapital' -delta).^2.*vGridCapital' + ...
                    diffPoisson;
        vectorDb = prodctvt*vGridCapital'.^aalpha - cb -...
                    bb/2*(cb./vGridCapital' -delta).^2.*vGridCapital' + ...
                    diffPoisson;
        vectorD0 = prodctvt*vGridCapital'.^aalpha-c0 + diffPoisson;
        
        invest = investmentVectorF.*If + investmentVectorb.*Ib +...
            investmentVector0.*I0;
        mPolicyFunction(:,nProductivity) = invest;
        dn = vectorDf.*If + vectorDb.*Ib +...
            vectorD0.*I0;
        %CONSTRUCT MATRIX
        X = -min(mub,0)/dk;
        Y = -max(muf,0)/dk + min(mub,0)/dk;
        Z = max(muf,0)/dk;
        %sparse matrix: faster
        AAA =spdiags(Y,0,I,I)+spdiags(X(2:I),-1,I,I)+spdiags([0;Z(1:I-1)],1,I,I);
        B = (r + 1/delta)*speye(I) - AAA + speye(I);
        vecb = dn + 1/deltaV*V;
        Vnew = B\vecb; %SOLVE SYSTEM OF EQUATIONS
        mValueFunctionNew(:,nProductivity) = Vnew;
        mdVf(:,nProductivity) = dVf;
        mdVb(:,nProductivity) = dVb;
        mPolicyFunction(:,nProductivity) = invest;
        pause
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

subplot(3,1,1)
plot(vGridCapital(1:nGridCapital/2),mValueFunction(1:nGridCapital/2,:))
xlim([vGridCapital(1) vGridCapital(nGridCapital/2)])
xlabel('k')
title('Value Function, continuous time')

subplot(3,1,2)
plot(vGridCapital(1:nGridCapital/2),mPolicyFunction(1:nGridCapital/2,:))
xlim([vGridCapital(1) vGridCapital(nGridCapital/2)])
xlabel('k')
title('Policy Function = Investment (continuous time)')

% vExactPolicyFunction = aalpha*bbeta.*(vGridCapital.^aalpha);
% 
subplot(3,1,3)
plot(vGridCapital(1:nGridCapital/2),mPolicyFunction(1:nGridCapital/2,:)-delta*vGridCapital(1:nGridCapital/2)')
xlim([vGridCapital(1) vGridCapital(nGridCapital/2)])
xlabel('i(k) - delta*k')
title(' dk/dt  = I(k) - delta*k')

%set(gcf,'PaperOrientation','landscape','PaperPosition',...
%[-0.9 -0.5 12.75 9])
%print('-dpdf','Figure1.pdf')

%% Plotting optimal investment and financing for lowest, avg and highest
% % mOptimalInvestment = zeros(nGridProductivity, nGridCapital);
% % mFinancing = zeros(nGridProductivity, nGridCapital);
% % for nProductivity = 1:nGridProductivity
% %     for nCapital = 1:nGridCapital
% %         prodctvt = exp(vProductivity(nProductivity));
% %         currentCapital = vGridCapital(nCapital);
% %         captomorrow = mPolicyFunction(nCapital,nProductivity);
% %         invstmnt = investment(prodctvt,currentCapital,captomorrow);
% %         mOptimalInvestment(nCapital,nProductivity) = invstmnt;
% %         profits = profitMatrix(nCapital,nProductivity);
% %         costOfInvstmnt = phi(prodctvt,currentCapital,captomorrow);
% %         dividend = profits - invstmnt - costOfInvstmnt;
% %         mFinancing(nCapital,nProductivity) = max(-dividend,0);
% %     end
% % end
% % %shock a
% % figure(3)
% % % subplot(2,1,1)
% % plot(vGridCapital(1:100),mOptimalInvestment(1:100,1))
% % xlabel('k')
% % title('d)01) Optimal Invesment for the mean shock only')

% hold on
% plot(vGridCapital,mOptimalInvestment(:,5))
% plot(vGridCapital,mOptimalInvestment(:,1))
% legend('high a', 'middle a', 'low a')
% hold off

% subplot(2,1,2)
% plot(vGridCapital(1:100),mFinancing(1:100,1))
% xlabel('k')
% title('d) Financing (-d(a,k) when positive) for the mean shock only')

% hold on
% plot(vGridCapital,mFinancing(:,5))
% plot(vGridCapital,mFinancing(:,1))
% legend('high a', 'middle a', 'low a')
% hold off