% This function is written to verify SLIDE's estimation algorithm.
% Param sp sn AMp AMn

clear all; close all; clc;

OCVp = load('../data/OCVfit_cathode.csv'); % stochiometry vs cathode potential

OCVn = load('../data/OCVfit_anode.csv'); % stochiometry vs anode potential

OCV  = load('../data/OCVfit_cell.csv'); % Ah vs full cell OCV
OCV(:,1) = (16/OCV(end,1))*OCV(:,1);

Ah_cell = OCV(:,1);
n = 1;
F = 96487;      %!< Faraday constant.
cmaxp = 51385;  %!< maximum li-concentration in the cathode [mol m-3]
cmaxn = 30555;  %!< maximum li-concentration in the anode [mol m-3]

cap = Ah_cell(end);

AMp_guess = cap*3600/(n*F*cmaxp);
AMn_guess = cap*3600/(n*F*cmaxn);

lb = [0;    0.5;  0.01*AMp_guess; 0.01*AMn_guess];
ub = [0.5;    1;     5*AMp_guess;    5*AMn_guess];


x0 = [0.38912; 0.603952; 3.59E-06; 5.77E-06];
x02 = [0.38912; 0.463952; 3.59E-06; 10.77E-06];
x03 = [0.25, 0.75, 10e-6, 10e-6];

f = @(x) residualFun2(x, OCVp, OCVn, OCV);
f2 = @(x) norm(f(x))^2;

f3 = @(x) norm(residualFun(x, OCVp, OCVn, OCV))^2;

x3 = [0.3853; 0.5645; 3.5997e-06; 6.1982e-06]; % BEST solution

x_best_16Ah = [0.3853, 0.5647, 2.1171e-05, 3.6459e-05]; % Best solution if Cell was 16 Ah. AMp and AMn just scaled by 16/2.71 = 5.88

f3(x3)
%f2(x4);

%%

%[x,resnorm,residual,exitflag,output] = lsqnonlin(f,x03([1,3]),lb([1,3]),ub([1,3]),optimoptions('lsqnonlin',Display='iter-detailed'));


x_fmincon = fmincon(f2,[0.3853, 2.1171e-05],[],[],[],[],lb([1,3]),ub([1,3]),[],optimoptions('fmincon',Display='iter'));


f2(x_fmincon);
x_particle = particleswarm(f2, 2, lb([1,3]), ub([1,3]), optimoptions('particleswarm','Display','iter','SwarmSize',5000));


function residuals = residualFun(x, OCVp, OCVn, OCV)
n = 1;
F = 96487;      %!< Faraday constant.
cmaxp = 51385;  %!< maximum li-concentration in the cathode [mol m-3]
cmaxn = 30555;  %!< maximum li-concentration in the anode [mol m-3]

As = OCV(:,1)*3600;

sp0 = x(1);
sn0 = x(2);

AMp = x(3);
AMn = x(4);

sp = sp0 + As/(n*F)/AMp/cmaxp;
sn = sn0 - As/(n*F)/AMn/cmaxn;

if(sp(end)>1 || sn(end)<0)
    residuals = inf(size(sp));
else
   ocvpi = interp1(OCVp(:,1),OCVp(:,2), sp);
   ocvni = interp1(OCVn(:,1),OCVn(:,2), sn);
   ocv = ocvpi - ocvni;

   residuals = OCV(:,2) - ocv;
end



end


function residuals = residualFun2(x, OCVp, OCVn, OCV)
n = 1;
F = 96487;      %!< Faraday constant.
cmaxp = 51385;  %!< maximum li-concentration in the cathode [mol m-3]
cmaxn = 30555;  %!< maximum li-concentration in the anode [mol m-3]

As = OCV(:,1)*3600;

sp0 = x(1);
ocvpi0 = interp1(OCVp(:,1),OCVp(:,2), sp0);
sn0 = interp1(OCVn(:,2), OCVn(:,1), ocvpi0 - OCV(1,2)); 

if(isnan(sn0))
    residuals = inf(size(OCV(:,1)));
    return;
end

AMp = x(2);

sp_end = sp0 + As(end)/(n*F)/AMp/cmaxp;

if(sp_end >1)
   residuals = inf(size(OCV(:,1)));
    return;
end


ocvpi_end = interp1(OCVp(:,1),OCVp(:,2), sp_end);
sn_end = interp1(OCVn(:,2), OCVn(:,1), ocvpi_end - OCV(end,2)); 


if(isnan(sn_end))
    residuals = inf(size(OCV(:,1)));
    return;
end

AMn = (As(end)/(n*F)) /(cmaxn*(sn0 - sn_end));

sp = sp0 + As/(n*F)/AMp/cmaxp;
sn = sn0 - As/(n*F)/AMn/cmaxn;



ocvpi = interp1(OCVp(:,1),OCVp(:,2), sp,'linear',0);
ocvni = interp1(OCVn(:,1),OCVn(:,2), sn,'linear',0);
ocv = ocvpi - ocvni;

residuals = OCV(:,2) - ocv;




end


