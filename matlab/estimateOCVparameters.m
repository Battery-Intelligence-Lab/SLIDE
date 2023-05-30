% This function is written to verify SLIDE's estimation algorithm.
% Param sp sn AMp AMn

clear variables; close all; clc;

OCVp = load('../data/OCVfit_cathode.csv'); % stochiometry vs cathode potential
OCVn = load('../data/OCVfit_anode.csv'); % stochiometry vs anode potential
OCV  = load('../data/OCVfit_cell.csv'); % Ah vs full cell OCV

OCV(:,1) = (16/OCV(end,1))*OCV(:,1);

Ah_cell = OCV(:,1);
As = Ah_cell*3600;
n = 1;
F = 96487;      %!< Faraday constant.
cmaxp = 51385;  %!< maximum li-concentration in the cathode [mol m-3]
cmaxn = 30555;  %!< maximum li-concentration in the anode [mol m-3]

cap = Ah_cell(end);

AMp_guess = cap*3600/(n*F*cmaxp);
AMn_guess = cap*3600/(n*F*cmaxn);


% ------ Constraints ------
% x0 = [sp, sn, 1/AMp, 1/AMn];
% Boundaries:
lb = [0;    0.2;  1/(5*AMp_guess); 1/(5*AMn_guess)];
ub = [0.8;    1;   1/(0.01*AMp_guess)   ;   1/(0.01*AMn_guess) ];

% Linear constraints:
% sp = x(1);  sn = x(2);
% AMp = 1/x(3); AMn = 1/x(4);
%
% sp(end)<1
% sp(end) = sp0 + As(end)/(n*F)/AMp/cmaxp
% x(1) + (As(end)/(n*F*cmaxp))*(1/AMp)  < 1
% x(1) + (As(end)/(n*F*cmaxp))*x(3)  < 1
%
%
% sn(end)>0
% sn(end) = sn0 - As/(n*F)/AMn/cmaxn
% x(2) - (As(end)/(n*F*cmaxn))*(1/AMn) > 0
% -x(2) + (As(end)/(n*F*cmaxn))*x(4) < 0

A = [1, 0, (As(end)/(n*F*cmaxp)), 0;
    0, -1, 0, (As(end)/(n*F*cmaxn))];

b = [1;-0.01];

% -------------------------

x0 = [0.38912; 0.603952; 1/3.59E-06; 1/5.77E-06];
x02 = [0.38912; 0.463952; 1/3.59E-06; 1/10.77E-06];
x03 = [0.25, 0.75, 10e-6, 10e-6];

f = @(x) norm(residualFun(x, OCVp, OCVn, OCV),2);

x3 = [0.3853; 0.5645; 3.5997e-06; 6.1982e-06]; % BEST solution

x_best_16Ah = [0.3853; 0.5647; 1/2.1171e-05; 1/3.6459e-05]; % Best solution if Cell was 16 Ah. AMp and AMn just scaled by 16/2.71 = 5.88
x_best_2 =  [0.385208057798616;   0.564328197343929; 1/2.116379738e-05; 1/3.648106966443356e-05];

x0 = x_best_2;
% Values to use:

% sp(1) = 0.385208057798616;
% sp(end) = 0.934146599395399;
% sp(50%) = 0.659677328597008;
%
% sn(1) = 0.564328197343929;
% sn(end) = 0.028773505826926;
% sn(50%) = 0.296550851585427;

gStats = @(x) guessStats(x, OCVp, OCVn, OCV);


%x0 = x3;

%%

%[x,resnorm,residual,exitflag,output] = lsqnonlin(f,x03([1,3]),lb([1,3]),ub([1,3]),optimoptions('lsqnonlin',Display='iter-detailed'));


% %x_fmincon = fmincon(f,x0,A,b,[],[],lb,ub,[],optimoptions('fmincon',Display='iter'));
gs = GlobalSearch();

problem = createOptimProblem('fmincon','x0',x0,'objective',f,'lb',lb,'ub',ub, 'Aineq',A, 'bineq',b, 'options', optimoptions('fmincon',Display='iter'));
x_fmincon = run(gs,problem);

%x_fmincon = lsqnonlin(f,x0,lb,ub,optimoptions('lsqnonlin',Display='iter-detailed'));




f(x_fmincon)
gStats(x_fmincon);
fprintf('\nStats for previous best:\n')
gStats(x_best_2);
%x_particle = particleswarm(f, length(x0), lb, ub, optimoptions('particleswarm','Display','iter','SwarmSize',20000));



function residuals = residualFun2(x, OCVp, OCVn, OCV)
% Only sp0 and 1/AMp.
n = 1;
F = 96487;      %!< Faraday constant.
cmaxp = 51385;  %!< maximum li-concentration in the cathode [mol m-3]
cmaxn = 30555;  %!< maximum li-concentration in the anode [mol m-3]

As = OCV(:,1)*3600;

sp0 = x(1);
AMp = 1/x(2);

% Then other values are found from boundary conditions.
% Assume initial values are same.

ocvpi0 = interp1(OCVp(:,1),OCVp(:,2), sp0);

OCV(:,2)





end


function residuals = residualFun(x, OCVp, OCVn, OCV)
n = 1;
F = 96487;      %!< Faraday constant.
cmaxp = 51385;  %!< maximum li-concentration in the cathode [mol m-3]
cmaxn = 30555;  %!< maximum li-concentration in the anode [mol m-3]

As = OCV(:,1)*3600;

sp0 = x(1);
sn0 = x(2);

AMp = 1/x(3);
AMn = 1/x(4);

sp = sp0 + As/(n*F)/AMp/cmaxp;
sn = sn0 - As/(n*F)/AMn/cmaxn;

if(sp(end)>1 || sn(end)<0)
    residuals = [OCV(:,2); diff(OCV(:,2))];
    residuals = OCV(:,2);
else
    ocvpi = interp1(OCVp(:,1),OCVp(:,2), sp);
    ocvni = interp1(OCVn(:,1),OCVn(:,2), sn);
    ocv = ocvpi - ocvni;

    %  residuals = [OCV(:,2)-ocv; diff(OCV(:,2)) - diff(ocv)];
    residuals = OCV(:,2)-ocv;

end

end


function [] = guessStats(x, OCVp, OCVn, OCV)

n = 1;
F = 96487;      %!< Faraday constant.
cmaxp = 51385;  %!< maximum li-concentration in the cathode [mol m-3]
cmaxn = 30555;  %!< maximum li-concentration in the anode [mol m-3]

As = OCV(:,1)*3600;

sp0 = x(1);
sn0 = x(2);

AMp = 1/x(3);
AMn = 1/x(4);

sp = sp0 + As/(n*F)/AMp/cmaxp;
sn = sn0 - As/(n*F)/AMn/cmaxn;

if(sp(end)>1 || sn(end)<0)
    residuals = OCV(:,2);
else
    ocvpi = interp1(OCVp(:,1),OCVp(:,2), sp);
    ocvni = interp1(OCVn(:,1),OCVn(:,2), sn);
    ocv = ocvpi - ocvni;

    residuals = OCV(:,2) - ocv;
end

figure;
plot(As/3600, OCV(:,2)); hold on;
plot(As/3600, ocv,'--');
xlabel('Charge throughput (Ah)');
ylabel('OCV (V)');
grid on;
legend('OCV','fitting')

RMSE = sqrt(mean(residuals.^2));

fprintf('sp range: [%4.4f, %4.4f]\n', sp(1), sp(end));
fprintf('sn range: [%4.4f, %4.4f]\n', sn(1), sn(end));

fprintf('RMSE is : %4.4f\n',RMSE);

end