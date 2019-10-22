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
% vProductivity = 1;
nGridProductivity = length(vProductivity);
na = nGridProductivity;  %number of points for logaGrid = [loga0, loga1... loga9]
if na >1
    mTransition = [1/2 1/2 0; 1/4 1/2 1/4; 0 1/2 1/2]';
else
    mTransition = 1;
end

% % Transition matrix % Tauchen method..
% if nGridProductivity > 1
%     P   = eye(na);
%     a_1 = vProductivity(1);
%     a_n = vProductivity(na);
%     for j = 1:na
%         aj = vProductivity(j);
%         upperBoundA = (a_1 - rho*aj - (1-rho)*logabar +deltaLoga/2)/sigma;
%         P(1,j) = normcdf(upperBoundA);
%         lowerboundA = (a_n - rho*aj - (1-rho)*logabar -deltaLoga/2)/sigma;
%         P(na,j) = 1-normcdf(lowerboundA);
%     end
%     for i = 2:(na-1)
%         for j = 1:(na)
%             ai = vProductivity(i);
%             aj = vProductivity(j);
%             upperBoundA = (ai - rho*aj - (1-rho)*logabar +deltaLoga/2)/sigma;
%             lowerboundA = (ai - rho*aj - (1-rho)*logabar -deltaLoga/2)/sigma;
%             P(i,j) = normcdf(upperBoundA)-normcdf(lowerboundA);
%         end
%     end
%     P'
%     mTransition   = P;
% else
%     mTransition = 1;
% end


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
vGridCapital = linspace(0.01*capitalSteadyState,7,200);
% vGridCapital = linspace(0.01*capitalSteadyState,6,200);

nGridCapital = length(vGridCapital);
nGridProductivity = length(vProductivity);

%% 3. Required matrices and vectors

mOutput           = zeros(nGridCapital,nGridProductivity);
% mValueFunction    = zeros(nGridCapital,nGridProductivity);
mValueFunction    = linspace(1,vGridCapital(nGridCapital),nGridCapital)'*vProductivity;
% mValueFunction    = ones(nGridCapital,nGridProductivity);
mValueFunctionNew = mValueFunction;
mPolicyFunction   = zeros(nGridCapital,nGridProductivity); %investment

%% 4. We pre-build output for each point in the grid

%profitMatrix = profit(exp(vProductivity),vGridCapital);

%% 5. Main iteration

maxDifference = 10.0;
tolerance = 0.001;%0.0000001;
iteration = 0;
deltaV = 0.0001;

vectorD = zeros(nGridCapital,1);
%Sparse matrix for rv = d + B*v
BB = spalloc(nGridCapital-2,nGridCapital-2,4+3*(nGridCapital-4));
B = zeros(nGridCapital,nGridCapital);
explicit = 1; % 1 -to run the explicit, 2- implicit

while (maxDifference>tolerance) 
%     iteration
%     maxDifference
    %expectedValueFunction = mValueFunction*mTransition';
    for nProductivity = 1:nGridProductivity
        vn = mValueFunction(2:end-1,nProductivity);
        for nCapital = (2):(nGridCapital-1)
            prodctvt = vProductivity(nProductivity);
            currentCapital = vGridCapital(nCapital);
            %left approach = Backward Approach
            deltakBackward = currentCapital - vGridCapital(nCapital-1);
            derivV = (mValueFunction(nCapital,nProductivity) - ...
                mValueFunction(nCapital-1,nProductivity)) /...
                (deltakBackward);
            ViBackward = derivV;
            ViForward = 0;
            ViCentral = 0;
            obj = @(x) -prodctvt*currentCapital^aalpha + x +...
                bb/2*(x/currentCapital -delta)^2*currentCapital - ...
                derivV*(x-delta*currentCapital);
            [invest,maxD] = fmincon(obj,0,[],[],[],[],[],[],[],options);
            %one guess for fmincon : delta*currentCapital
            maxD = -maxD;
%             if derivV > 0
%                 display('got positive investment, yay!')
%                 maxD
%                 temporary = invest - delta*currentCapital
%                 pause;
%             end
            if( (invest -delta*currentCapital) <= 0) %use Backward
                kdot = (invest-delta*currentCapital)/ deltakBackward;
                B(nCapital-1,nProductivity) = -kdot;
                B(nCapital,nProductivity) = kdot;
                B(nCapital+1,nProductivity) = 0;
            else
            %elseif( (invest -delta*currentCapital) > 0) %could have used Forward
                ViBackward = 0;
                deltakForward = vGridCapital(nCapital+1)-currentCapital;
                derivV = (mValueFunction(nCapital+1,nProductivity) - ...
                    mValueFunction(nCapital,nProductivity)) /...
                    (deltakForward);
                ViForward = derivV;
                obj = @(x) -prodctvt*currentCapital^aalpha + x +...
                    bb/2*(x/currentCapital -delta)^2*currentCapital - ...
                    derivV*(x-delta*currentCapital);
                [invest,maxD] = fmincon(obj,0,[],[],[],[],[],[],[],options);
                maxD = -maxD;
                if ((invest -delta*currentCapital) >= 0) %correct to have used Forward
                    kdot = (invest-delta*currentCapital)/ deltakForward;
                    B(nCapital-1,nProductivity) = 0;
                    B(nCapital,nProductivity) = -kdot;
                    B(nCapital+1,nProductivity) = kdot;
                else%should have used Central
                    deltakCentral = vGridCapital(nCapital+1)-vGridCapital(nCapital-1);
                    derivV = (mValueFunction(nCapital+1,nProductivity) - ...
                        mValueFunction(nCapital-1,nProductivity)) /...
                        (deltakCentral);
                    ViCentral = derivV;
                    obj = @(x) -prodctvt*currentCapital^aalpha + x +...
                        bb/2*(x/currentCapital -delta)^2*currentCapital - ...
                        derivV*(x-delta*currentCapital);
                    [invest,maxD] = fmincon(obj,0,[],[],[],[],[],[],[],options);
                    maxD = -maxD;
                    kdot = (invest-delta*currentCapital)/ deltakCentral;
                    B(nCapital-1,nProductivity) = -kdot;
                    B(nCapital,nProductivity) = 0;
                    B(nCapital+1,nProductivity) = kdot;
                end
            end
            diffPoisson = mValueFunction(nCapital,:)*mTransition(:,nProductivity);
            vectorD(nCapital) = maxD-derivV*(invest-delta*currentCapital)+diffPoisson;
            mPolicyFunction(nCapital,nProductivity) = invest;
        end
        dn = vectorD(2:end-1);
        BB = sparse( B(2:nGridCapital-1,2:nGridCapital-1) - eye(nGridCapital-2) );
%         deltaV = 0.001;
        %Explicit Euler
        if explicit == 1
            vnext = deltaV*(dn + BB*vn - r*vn)+vn;
        else %implicit Euler
            vnext = ((r+1/deltaV)*eye(nGridCapital-2) - BB)\(dn + 1/deltaV*vn);
        end
            
        mValueFunctionNew(2:end-1,nProductivity) = vnext;
%         mValueFunctionNew(1,nProductivity) = vnext(1);
%         mValueFunctionNew(end,nProductivity) = vnext(end);
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