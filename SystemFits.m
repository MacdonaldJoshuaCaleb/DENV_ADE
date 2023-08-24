function sol= SystemFits(option)
set(0,'defaultTextInterpreter','latex');
close all
% option - select serotype
% 1 -> DENV1
% 2 -> DENV2
% 3 -> DENV3


% for weighted fits, primary
function rr = paramfunPrim(p, t)
    r = p(1);
     a1 = a1_init;   % Viral growth enhancing rate induced by cross-reactiveantibody-virus binding (ADE)
    a2 = (3/2).*a1;   % Cross-reactive antibody-virus killing rate upon cooperativebinding
    d = p(2);   % Specific antibody-virus killing rate upon binding
    f1 = p(3);   % Cross-reactive antibody activation rate
    f2 = (5/4).*p(3);  % Specific antibody activation rate
    k1 = 0;   % Antibody interference competition coeficient
    k2 = 0;  % Antibody interference competition coeficient
   A1 = (1/10).*p(4);    % Saturation coeficients of Hill functions for cross-reactive antibody
    A2 = (1/2).*p(5);  % Saturation coeficients of Hill functions for cross-reactive antibody
    C1 = (1/5).*p(6);   % Saturation coeficients of Hill functions for cross-reactive antibody
    B = (1/20).*p(7);  % Saturation coefficients of Hill functions for specific antibody
    C2 = (1/10).*p(8);  % Saturation coefficients of Hill functions for specific antibody

    u0 = [v0,p(9),p(10)]; % initial conditions 
    f = @(t,u) [u(1)*(r + (a1*u(2))/(A1+u(2)) - (a2*u(2).^2)/(A2+u(2).^2) - (d*u(3))/(B+u(3)));
                  (f1*u(1)*u(2))/(C1+u(2)+k1*u(3));
                  (f2*u(1)*u(3))/(C2+k2*u(2)+u(3))];
    [~,rr] = ode45(f,t,u0);
    rr = rr(:,1)'.*swt1;
end

% for unweighted fits, primary  
function rr = paramfunPrimCheck(p, t)
    r = p(1);
     a1 = p(15);   % Viral growth enhancing rate induced by cross-reactiveantibody-virus binding (ADE)
    a2 = (3/2).*a1;   % Cross-reactive antibody-virus killing rate upon cooperativebinding
    d = p(2);   % Specific antibody-virus killing rate upon binding
    f1 = p(3);   % Cross-reactive antibody activation rate
    f2 = (5/4).*p(3);  % Specific antibody activation rate
    k1 = 0;   % Antibody interference competition coeficient
    k2 = 0;  % Antibody interference competition coeficient
   A1 = (1/10).*p(4);    % Saturation coeficients of Hill functions for cross-reactive antibody
    A2 = (1/2).*p(5);  % Saturation coeficients of Hill functions for cross-reactive antibody
    C1 = (1/5).*p(6);   % Saturation coeficients of Hill functions for cross-reactive antibody
    B = (1/20).*p(7);  % Saturation coefficients of Hill functions for specific antibody
    C2 = (1/10).*p(8);  % Saturation coefficients of Hill functions for specific antibody

    u0 = [p(14),p(9),p(10)]; % initial conditions 
    f = @(t,u) [u(1)*(r + (a1*u(2))/(A1+u(2)) - (a2*u(2).^2)/(A2+u(2).^2) - (d*u(3))/(B+u(3)));
                  (f1*u(1)*u(2))/(C1+u(2)+k1*u(3));
                  (f2*u(1)*u(3))/(C2+k2*u(2)+u(3))];
    [~,rr] = ode45(f,t,u0);
    %rr = [rr(:,1)'.*swt1,rr(end,2)'.*swt2,rr(end,3)'.*swt3];
end

% for weighted fits, secondary  
function rr = paramfunSec(p, t)   
    r = p(1);
     a1 = a1_init;   % Viral growth enhancing rate induced by cross-reactiveantibody-virus binding (ADE)
    a2 = (3/2).*a1;   % Cross-reactive antiboceil(.01*length(U))dy-virus killing rate upon cooperativebinding
    d = p(2);   % Specific antibody-virus killing rate upon binding
    f1 = p(3);   % Cross-reactive antibody activation rate
    f2 = (5/4).*p(3);  % Specific antibody activation rate
    k1 = 0;   % Antibody interference competition coeficient
    k2 = 0;  % Antibody interference competition coeficient
   A1 = (1/10).*p(4);    % Saturation coeficients of Hill functions for cross-reactive antibody
    A2 = (1/2).*p(5);  % Saturation coeficients of Hill functions for cross-reactive antibody
    C1 = (1/5).*p(6);   % Saturation coeficients of Hill functions for cross-reactive antibody
    B = (1/20).*p(7);  % Saturation coefficients of Hill functions for specific antibody
    C2 = (1/10).*p(8);  % Saturation coefficients of Hill functions for specific antibody


    u0 = [v0,p(11),p(10)]; % initial conditions 
    f = @(t,u) [u(1)*(r + (a1*u(2))/(A1+u(2)) - (a2*u(2).^2)/(A2+u(2).^2) - (d*u(3))/(B+u(3)));
                  (f1*u(1)*u(2))/(C1+u(2)+k1*u(3));
                  (f2*u(1)*u(3))/(C2+k2*u(2)+u(3))];
    [~,rr] = ode45(f,t,u0);
    rr = rr(:,1)'.*swt2;
end

% for unweighted fits secondary 
function rr = paramfunSecCheck(p, t)   
    r = p(1);
     a1 = p(15);   % Viral growth enhancing rate induced by cross-reactiveantibody-virus binding (ADE)
    a2 = (3/2).*a1;   % Cross-reactive antibody-virus killing rate upon cooperativebinding
    d = p(2);   % Specific antibody-virus killing rate upon binding
    f1 = p(3);   % Cross-reactive antibody activation rate
    f2 = (5/4).*p(3);  % Specific antibody activation rate
    k1 = 0;   % Antibody interference competition coeficient
    k2 = 0;  % Antibody interference competition coeficient
    A1 = (1/10).*p(4);    % Saturation coeficients of Hill functions for cross-reactive antibody
    A2 = (1/2).*p(5);  % Saturation coeficients of Hill functions for cross-reactive antibody
    C1 = (1/5).*p(6);   % Saturation coeficients of Hill functions for cross-reactive antibody
    B = (1/20).*p(7);  % Saturation coefficients of Hill functions for specific antibody
    C2 = (1/10).*p(8);  % Saturation coefficients of Hill functions for specific antibody

    u0 = [p(14),p(11),p(10)]; % initial conditions 
    f = @(t,u) [u(1)*(r + (a1*u(2))/(A1+u(2)) - (a2*u(2).^2)/(A2+u(2).^2) - (d*u(3))/(B+u(3)));
                  (f1*u(1)*u(2))/(C1+u(2)+k1*u(3));
                  (f2*u(1)*u(3))/(C2+k2*u(2)+u(3))];
    [~,rr] = ode45(f,t,u0);
end

function ret = objfun(p,t)
ret1 = paramfunPrim(p,t(idxP)); % primary strain fit
ret2 = paramfunSec(p,t(idxS)); % secondary strain fitt
ret = [ret1,ret2];
end

data = readtable('DengueDataRaw.csv');
% extract primary and secondary paitent data 
ind1p = find(string(table2array(data(:,3))) == 'DENV1' & string(table2array(data(:,4))) == 'primary');
ind2p = find(string(table2array(data(:,3))) == 'DENV2' & string(table2array(data(:,4))) == 'primary');
ind3p = find(string(table2array(data(:,3))) == 'DENV3' & string(table2array(data(:,4))) == 'primary');
ind1s = find(string(table2array(data(:,3))) == 'DENV1' & string(table2array(data(:,4))) == 'secondary');
ind2s = find(string(table2array(data(:,3))) == 'DENV2' & string(table2array(data(:,4))) == 'secondary');
ind3s = find(string(table2array(data(:,3))) == 'DENV3' & string(table2array(data(:,4))) == 'secondary');

data1p = data(ind1p,:);
data2p = data(ind2p,:);
data3p = data(ind3p,:);

data1s = data(ind1s,:);
data2s = data(ind2s,:);
data3s = data(ind3s,:);

%%

% extract primary and secondary viral load data and associated times
TT1p = table2array(data1p(:,6))./24;
TT2p = table2array(data2p(:,6))./24;
TT3p = table2array(data3p(:,6))./24;

TT1s = table2array(data1s(:,6))./24;
TT2s = table2array(data2s(:,6))./24;
TT3s = table2array(data3s(:,6))./24;

VV1p = log10(table2array(data1p(:,7)));
VV2p = log10(table2array(data2p(:,7)));
VV3p = log10(table2array(data3p(:,7)));

VV1s = log10(table2array(data1s(:,7)));
VV2s = log10(table2array(data2s(:,7)));
VV3s = log10(table2array(data3s(:,7)));

% extract individual paitents 
indD1Patprim = find(table2array(data1p(:,5)) == 1);
indD2Patprim = find(table2array(data2p(:,5)) == 1);
indD3Patprim = find(table2array(data3p(:,5)) == 1);


indD1Patsec = find(table2array(data1s(:,5)) == 1);
indD2Patsec = find(table2array(data2s(:,5)) == 1);
indD3Patsec = find(table2array(data3s(:,5)) == 1);

if option == 1
    indDPatprim = indD1Patprim;
    indDPatsec = indD1Patsec;
elseif option == 2
    indDPatprim = indD2Patprim;
    indDPatsec = indD2Patsec;
elseif option == 3
    indDPatprim = indD3Patprim;
    indDPatsec = indD3Patsec;
end

count = 1;
% loop over all combinations of hosts 
for rp = 1:1 %length(indDPatprim)
%for rp = 1:1
% select primary paitent an remove values w/o positive virus assay value
if rp < length(indDPatprim) && option == 1
    patTprim = TT1p(indDPatprim(rp):indDPatprim(rp+1)-1);
    patVprim = VV1p(indDPatprim(rp):indDPatprim(rp+1)-1);
    PtNosPrim = table2array(data1p(indDPatprim(rp),1));
elseif rp == length(indDPatprim) && option == 1
    patTprim = TT1p(indDPatprim(rp):end);
    patVprim = VV1p(indDPatprim(rp):end);
    PtNosPrim = table2array(data1p(indDPatprim(rp),1));
elseif rp < length(indDPatprim) && option == 2
    patTprim = TT2p(indDPatprim(rp):indDPatprim(rp+1)-1);
    patVprim = VV2p(indDPatprim(rp):indDPatprim(rp+1)-1);
    PtNosPrim = table2array(data2p(indDPatprim(rp),1));
elseif rp == length(indDPatprim) && option == 2
    patTprim = TT2p(indDPatprim(rp):end);
    patVprim = VV2p(indDPatprim(rp):end);
    PtNosPrim = table2array(data2p(indDPatprim(rp),1));
elseif rp < length(indDPatprim) && option == 3
    patTprim = TT2p(indDPatprim(rp):indDPatprim(rp+1)-1);
    patVprim = VV2p(indDPatprim(rp):indDPatprim(rp+1)-1);
    PtNosPrim = table2array(data3p(indDPatprim(rp),1));
elseif rp == length(indDPatprim) && option == 3
    patTprim = TT3p(indDPatprim(rp):end);
    patVprim = VV3p(indDPatprim(rp):end);
    PtNosPrim = table2array(data3p(indDPatprim(rp),1));
end

temp = find(patVprim > 0);
patVprim = patVprim(patVprim > 0);
patTprim = patTprim(temp);



for rs = 1:1 %length(indDPatsec)
% select seconary paitent and remove values w/o positive virus assay value
if rs < length(indDPatsec) && option == 1
    patTsec = TT1s(indDPatsec(rs):indDPatsec(rs+1)-1);
    patVsec = VV1s(indDPatsec(rs):indDPatsec(rs+1)-1);
    PtNosSec = table2array(data1s(indDPatsec(rs),1));
elseif rs == length(indDPatsec) && option == 1
    patTsec = TT1s(indDPatsec(rs):end);
    patVsec = VV1s(indDPatsec(rs):end);
    PtNosSec = table2array(data1s(indDPatsec(rs),1));
elseif rs < length(indDPatsec) && option == 2
    patTsec = TT2s(indDPatsec(rs):indDPatsec(rs+1)-1);
    patVsec = VV2s(indDPatsec(rs):indDPatsec(rs+1)-1);
    PtNosSec = table2array(data2s(indDPatsec(rs),1));
elseif rs == length(indDPatsec) && option == 2
    patTsec = TT2s(indDPatsec(rs):end);
    patVsec = VV2s(indDPatsec(rs):end);
    PtNosSec = table2array(data2s(indDPatsec(rs),1));
elseif rs < length(indDPatsec) && option == 3
    patTsec = TT3s(indDPatsec(rs):indDPatsec(rs+1)-1);
    patVsec = VV3s(indDPatsec(rs):indDPatsec(rs+1)-1);
    PtNosSec = table2array(data3s(indDPatsec(rs),1));
elseif rs == length(indDPatsec) && option == 3
    patTsec = TT3s(indDPatsec(rs):end);
    patVsec = VV3s(indDPatsec(rs):end);
    PtNosSec = table2array(data3s(indDPatsec(rs),1));
end


temp = find(patVsec > 0);
patVsec = patVsec(patVsec > 0);
patTsec = patTsec(temp);

count2 = 1;
while count2 <= 500
    % initalize fits by setting bite time, intial viral load, and IgG
    % activation rate 
    % run for 100 sets of initilization parameters 
t10 = normrnd(-5,1/3);
t20 = normrnd(-5,1/3);
v0 = lognrnd(log(.01),.5);
a1_init = normrnd(4,1/3);
T1s = [t10; patTsec; min(12,patTsec(end)+8)]';
T1p = [t20; patTprim; min(12,patTprim(end)+8)]';
Vs = [0; patVsec; 0]'; 
Vp = [0; patVprim; 0]';

t = union(T1s,T1p);
% get a single vector of all time points 
[sharedvals,idxS] = intersect(t,T1s,'stable');
[sharedvals,idxP] = intersect(t,T1p,'stable');

% fitting stuff

    
swt1=ones(1,length(Vp));
idx = find(Vp > 0);
swt1(idx) = 1;
%swt1(idx(end)+1) = .5;
[M1 I1] = max(Vp);

swt1(I1) = 1.5;
    swt1(I1-1) = 1.5;
    swt1(I1+1) = 1.5;


swt2 = ones(1,length(Vs));
idx = find(Vs > 0);
%swt2(idx(end)+1) = .5;
swt2(idx) = 1;
[M2 I2] = max(Vs);
swt2(I2) = 1.5;

    swt2(I2-1) = 1.5;
    swt2(I2+1) = 1.5;


S = sort([max(Vp),max(Vs)]);

% allowed parameter bounds 
lb = [1, 1,0.05,.75*S(1),.75*S(1),.75*S(1),.75*S(1),.75*S(1),.05,.05,.05];
ub = [5, 5, 2,1.25.*S(2),1.25*S(2),1.25*S(2),1.25*S(2),1.25*S(2),.5,.5,1.5];
% initalize on first iterate of pair 
% randomize intial values to test for effect of variation in minimization starting point  
a = lb;
b = ub;
p0 = (b-a).*rand(1,11) + a;
% fit the model
[pfit, resnorm(count2)] = lsqcurvefit(@objfun,p0,t,[Vp.*swt1,Vs.*swt2],lb,ub);

% store the estimates 
pfitInd(count2,:) = [pfit,t10,t20,v0,a1_init];
tt1 = pfitInd(count2,12):.01:20;
tt2 = pfitInd(count2,13):.01:20;

ccP = paramfunPrimCheck(pfitInd(count2,1:15),tt1);
ccS = paramfunSecCheck(pfitInd(count2,1:15),tt2);
MP(count2) = max(ccP(:,1));
MS(count2) = max(ccS(:,1));
MGC(count2) = max([max(ccP(:,2)),max(ccS(:,2))]);
MGS(count2) = max([max(ccP(:,3)),max(ccS(:,3))]);
ccPF = paramfunPrimCheck(pfitInd(count2,1:15),T1p);
ccSF = paramfunSecCheck(pfitInd(count2,1:15),T1s);
SS(count2,:) = sqrt(sum(abs([ccPF(:,1)', ccSF(:,1)']-[Vp,Vs]).^2));
count2 = count2+1;
end

% remove infeasibly high maximum viral loads 
feasIDX = find(MS < 20 & MP < 20 & MGC < 15 & MGS < 15);
if length(feasIDX) < 1
    feasIDX = find(MS < 25 & MP < 25 & MGC < 15 & MGS < 15);
end
if length(feasIDX) < 1
    feasIDX = find(MS < 30 & MP < 30 & MGC < 15 & MGS < 15);
end
SS = SS(feasIDX,:);
pfitInd = pfitInd(feasIDX,:);
resnorm = resnorm(feasIDX);
MP = MP(feasIDX);
MS = MS(feasIDX);
% select from among fits based on how well maximum viral laod is fit 




[m1 I] = min(SS);

if rp == 1 && rs == 1
     [mS, IS] = min(SS)
    
    fig = figure;
     fig.Position = [1 1 800 800];
    hold on
    for j = 1:length(feasIDX)
        ts = pfitInd(j,12):.01:20;
        sol = paramfunPrimCheck(pfitInd(j,:),ts);
        
        plot(ts,10.^sol(:,1),'linewidth',.5,'LineStyle',':','Color',[235 	175 	63]/256)
        plot(ts,10.^sol(:,2),'linewidth',.5,'LineStyle',':','Color',[153 	153 	153]/256)
        plot(ts,10.^sol(:,3),'linewidth',.5,'LineStyle',':','Color',[204 	204 	204]/256)
    end
     ts = pfitInd(I,12):.01:20;
        sol = paramfunPrimCheck(pfitInd(I,:),ts);
       
            plot(ts,10.^sol(:,1),'linewidth',4,'LineStyle','-.','Color','black')
            plot(ts,10.^sol(:,2),'linewidth',4,'LineStyle','-.','Color','black')
            plot(ts,10.^sol(:,3),'linewidth',4,'LineStyle','-.','Color','black')

    scatter(patTprim,10.^patVprim,80,'filled','black')
    xlim([min([pfitInd(:,12);pfitInd(:,13)]),20])
    ylim([1 10.^16])
    hold off
    set(gca, 'YScale', 'log')
    ylabel('Concentration')
    xlabel('Time since symptom onset (days)')
    title('Primary DENV 3 - ID 8')
    set(gca,'FontSize',24)


    fig = figure;
     fig.Position = [1 1 800 800];
    hold on
    for j = 1:length(feasIDX)
        ts = pfitInd(j,13):.01:20;
        sol = paramfunSecCheck(pfitInd(j,:),ts);
        
        plot(ts,10.^sol(:,1),'linewidth',.5,'LineStyle',':','Color',[233 	87 	87]./256)
        plot(ts,10.^sol(:,2),'linewidth',.5,'LineStyle',':','Color',[153 	153 	153]/256)
        plot(ts,10.^sol(:,3),'linewidth',.5,'LineStyle',':','Color',[204 	204 	204]/256)
    end
     ts = pfitInd(I,13):.01:20;
        sol = paramfunPrimCheck(pfitInd(I,:),ts);
       
            plot(ts,10.^sol(:,1),'linewidth',4,'LineStyle','-.','Color','black')
            plot(ts,10.^sol(:,2),'linewidth',4,'LineStyle','-.','Color','black')
            plot(ts,10.^sol(:,3),'linewidth',4,'LineStyle','-.','Color','black')

    scatter(patTprim,10.^patVprim,80,'filled','black')
    xlim([min([pfitInd(:,12);pfitInd(:,13)]),20])
    ylim([1 10.^16])
    hold off
    set(gca, 'YScale', 'log')
    xlabel('Time since symptom onset (days)')
    ylabel('Concentration')
    title('Secondary DENV 3 - ID 31')
    set(gca,'FontSize',24)

    figure
    hold on
    histogram(pfitInd(:,12))
    xline(pfitInd(I,12),'linewidth',2,'color','red')
    hold off
    set(gca,'FontSize',24)
    xlabel('$\tau_0^p$')
    ylabel('Count')

    figure
    hold on
    histogram(pfitInd(:,13))
    xline(pfitInd(I,13),'linewidth',2,'color','red')
    hold off
    set(gca,'FontSize',24)
    xlabel('$\tau_0^s$')
    ylabel('Count')

   

end
    [mS, IS] = min(SS);
% extract maximum viral loads 
tt1F = pfitInd(I,12):.01:20;
tt2F = pfitInd(I,13):.01:20;
ccPF = paramfunPrimCheck(pfitInd(I,1:15),tt1F);
ccSF = paramfunSecCheck(pfitInd(I,1:15),tt2F);
MPF = max(ccPF(:,1));
MSF = max(ccSF(:,1));

% for each pair of hosts store parameter estimates, max viral loads, and
% paitent IDs
pfits(count,:) = [pfitInd(I,:),MPF,MSF,PtNosPrim,PtNosSec];

check = ceil(.2*length(indDPatprim)*length(indDPatsec));
if mod(count,check) == 0
   fprintf('-----------------------------------\n')
   fprintf('-----------------------------------\n')
   fprintf('percent complete: %d\n',round(100*count./(length(indDPatprim)*length(indDPatsec)),2))
   fprintf('-----------------------------------\n')
   fprintf('-----------------------------------\n')
end

% iterate the total number of fits 
count = count+1;

% print indices 
[rp,rs]


end
end
% return the parameter estimates after all combinations are fit 
sol = pfits;
end