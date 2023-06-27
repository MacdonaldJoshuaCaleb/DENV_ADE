function sol= SystemFits(option)
close all
% its - random paring iteration 
% option - select serotype
% 1 -> DENV1
% 2 -> DENV2
% 3 -> DENV3
% option2 - select return output 
% 1 -> paramters 
% 2 -> time trajectories


function rr = paramfunPrim(p, t)
    r = p(1);
     a1 = 4;   % Viral growth enhancing rate induced by cross-reactiveantibody-virus binding (ADE)
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

    u0 = [.01,.1,.1]; % initial conditions 
    f = @(t,u) [u(1)*(r + (a1*u(2))/(A1+u(2)) - (a2*u(2).^2)/(A2+u(2).^2) - (d*u(3))/(B+u(3)));
                  (f1*u(1)*u(2))/(C1+u(2)+k1*u(3));
                  (f2*u(1)*u(3))/(C2+k2*u(2)+u(3))];
    [~,rr] = ode45(f,t,u0);
    rr = rr(:,1)'.*swt1;
end

function rr = paramfunPrimCheck(p, t)
    r = p(1);
     a1 = 4;   % Viral growth enhancing rate induced by cross-reactiveantibody-virus binding (ADE)
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

    u0 = [.01,.1,.1]; % initial conditions 
    f = @(t,u) [u(1)*(r + (a1*u(2))/(A1+u(2)) - (a2*u(2).^2)/(A2+u(2).^2) - (d*u(3))/(B+u(3)));
                  (f1*u(1)*u(2))/(C1+u(2)+k1*u(3));
                  (f2*u(1)*u(3))/(C2+k2*u(2)+u(3))];
    [~,rr] = ode45(f,t,u0);
    %rr = [rr(:,1)'.*swt1,rr(end,2)'.*swt2,rr(end,3)'.*swt3];
end

function rr = paramfunSec(p, t)   
    r = p(1);
     a1 = 4;   % Viral growth enhancing rate induced by cross-reactiveantibody-virus binding (ADE)
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


    u0 = [.01,p(9),p(10)]; % initial conditions 
    f = @(t,u) [u(1)*(r + (a1*u(2))/(A1+u(2)) - (a2*u(2).^2)/(A2+u(2).^2) - (d*u(3))/(B+u(3)));
                  (f1*u(1)*u(2))/(C1+u(2)+k1*u(3));
                  (f2*u(1)*u(3))/(C2+k2*u(2)+u(3))];
    [~,rr] = ode45(f,t,u0);
    rr = rr(:,1)'.*swt2;
end

function rr = paramfunSecCheck(p, t)   
    r = p(1);
     a1 = 4;   % Viral growth enhancing rate induced by cross-reactiveantibody-virus binding (ADE)
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

    u0 = [.01,p(9),p(10)]; % initial conditions 
    f = @(t,u) [u(1)*(r + (a1*u(2))/(A1+u(2)) - (a2*u(2).^2)/(A2+u(2).^2) - (d*u(3))/(B+u(3)));
                  (f1*u(1)*u(2))/(C1+u(2)+k1*u(3));
                  (f2*u(1)*u(3))/(C2+k2*u(2)+u(3))];
    [~,rr] = ode45(f,t,u0);
end

function ret = objfun(p,t)
%p1 = p(1:2); % primary strain 1

ret1 = paramfunPrim(p,t(idxP)); % primary strain fit
ret2 = paramfunSec(p,t(idxS)); % secondary strain fitt
ret = [ret1,ret2];
end

function ans = plotfunP(p,tt,opt)
        
    r = p(1);
     a1 = 4;   % Viral growth enhancing rate induced by cross-reactiveantibody-virus binding (ADE)
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


    u0 = [.01,.1,.1]; % initial conditions 
    f = @(t,u) [u(1)*(r + (a1*u(2))/(A1+u(2)) - (a2*u(2).^2)/(A2+u(2).^2) - (d*u(3))/(B+u(3)));
                  (f1*u(1)*u(2))/(C1+u(2)+k1*u(3));
                  (f2*u(1)*u(3))/(C2+k2*u(2)+u(3))];
    [~,rr] = ode45(f,tt,u0);
    DP = diff(rr(:,1));
    
    if opt == 1
    figure
    hold on
    plot(tt,rr(:,1),'r')
    plot(tt,rr(:,2),'b')
    plot(tt,rr(:,3),'g')
    scatter(t(idxP(2:end-1)),Vp(2:end-1),'r')
    hold off
    title('Primary')
    xlabel('Time since symptom onset')
    ylabel('Concentration')
    ylim([0 15])
    xlim([-5 20])
    hold off
    id = find(abs(DP) > .025);
    Vstart = tt(id(1)+1);
    Vend = tt(id(end)+1);
    Vdur = Vend-Vstart;
    Q = trapz(tt,rr(:,1));
    [MV IMV] = max(rr(:,1));
    ans = [Q, MV, tt(IMV),Vstart,Vend,Vdur];
    end
    if opt == 2
        id = find(abs(DP) > .025);
    Vstart = tt(id(1)+1);
    Vend = tt(id(end)+1);
    Vdur = Vend-Vstart;
        Q = trapz(tt,rr(:,1));
    [MV IMV] = max(rr(:,1));
    ans = [Q, MV, tt(IMV),Vstart,Vend,Vdur];
    end
    if opt == 3
        ans = rr;
    end
end

function ans = plotfunS(p,tt,opt)
    
    r = p(1);
     a1 = 4;   % Viral growth enhancing rate induced by cross-reactiveantibody-virus binding (ADE)
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


    u0 = [.01,p(9),p(10)]; % initial conditions 
    f = @(t,u) [u(1)*(r + (a1*u(2))/(A1+u(2)) - (a2*u(2).^2)/(A2+u(2).^2) - (d*u(3))/(B+u(3)));
                  (f1*u(1)*u(2))/(C1+u(2)+k1*u(3));
                  (f2*u(1)*u(3))/(C2+k2*u(2)+u(3))];
    [~,rr] = ode45(f,tt,u0);
    DP = diff(rr(:,1));
    id = find(abs(DP) > .025);
    Vstart = tt(id(1)+1);
    Vend = tt(id(end)+1);
    Vdur = Vend-Vstart;
    if opt == 1
    figure
    hold on
    plot(tt,rr(:,1),'r')
    plot(tt,rr(:,2),'b')
    plot(tt,rr(:,3),'g')
    scatter(t(idxS(2:end-1)),Vs(2:end-1),'r')
    hold off
    title('Secondary')
    xlabel('Time since symptom onset')
    ylabel('Concentration')
    ylim([0 15])
    xlim([-5 20])
    
    
        Q = trapz(tt,rr(:,1));
   [MV IMV] = max(rr(:,1));
    ans = [Q, MV, tt(IMV),Vstart,Vend,Vdur];
    end
     if opt == 2
        Q = trapz(tt,rr(:,1));
    [MV IMV] = max(rr(:,1));
    ans = [Q, MV, tt(IMV),Vstart,Vend,Vdur];
     end
        if opt == 3
        ans = rr;
    end
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

for rp = 1:length(indDPatprim)
%for rp = 1:1
% select primary paitent an remove values w/o positive virus assay value
if rp < length(indDPatprim) && option == 1
    patTprim = TT1p(indDPatprim(rp):indDPatprim(rp+1)-1);
    patVprim = VV1p(indDPatprim(rp):indDPatprim(rp+1)-1);
elseif rp == length(indDPatprim) && option == 1
    patTprim = TT1p(indDPatprim(rp):end);
    patVprim = VV1p(indDPatprim(rp):end);
elseif rp < length(indDPatprim) && option == 2
    patTprim = TT2p(indDPatprim(rp):indDPatprim(rp+1)-1);
    patVprim = VV2p(indDPatprim(rp):indDPatprim(rp+1)-1);
elseif rp == length(indDPatprim) && option == 2
    patTprim = TT2p(indDPatprim(rp):end);
    patVprim = VV2p(indDPatprim(rp):end);
elseif rp < length(indDPatprim) && option == 3
    patTprim = TT2p(indDPatprim(rp):indDPatprim(rp+1)-1);
    patVprim = VV2p(indDPatprim(rp):indDPatprim(rp+1)-1);
elseif rp == length(indDPatprim) && option == 3
    patTprim = TT3p(indDPatprim(rp):end);
    patVprim = VV3p(indDPatprim(rp):end);
end

temp = find(patVprim > 0);
patVprim = patVprim(patVprim > 0);
patTprim = patTprim(temp);



for rs = 1:length(indDPatsec)
% select seconary paitent and remove values w/o positive virus assay value
if rs < length(indDPatsec) && option == 1
    patTsec = TT1s(indDPatsec(rs):indDPatsec(rs+1)-1);
    patVsec = VV1s(indDPatsec(rs):indDPatsec(rs+1)-1);
elseif rs == length(indDPatsec) && option == 1
    patTsec = TT1s(indDPatsec(rs):end);
    patVsec = VV1s(indDPatsec(rs):end);
elseif rs < length(indDPatsec) && option == 2
    patTsec = TT2s(indDPatsec(rs):indDPatsec(rs+1)-1);
    patVsec = VV2s(indDPatsec(rs):indDPatsec(rs+1)-1);
elseif rs == length(indDPatsec) && option == 2
    patTsec = TT2s(indDPatsec(rs):end);
    patVsec = VV2s(indDPatsec(rs):end);
elseif rs < length(indDPatsec) && option == 3
    patTsec = TT3s(indDPatsec(rs):indDPatsec(rs+1)-1);
    patVsec = VV3s(indDPatsec(rs):indDPatsec(rs+1)-1);
elseif rs == length(indDPatsec) && option == 3
    patTsec = TT3s(indDPatsec(rs):end);
    patVsec = VV3s(indDPatsec(rs):end);
end


temp = find(patVsec > 0);
patVsec = patVsec(patVsec > 0);
patTsec = patTsec(temp);

count2 = 1;
while count2 <= 100
t10 = normrnd(-5,1/3);
t20 = normrnd(-5,1/3);
% a = -6;
% b = -4;
% % b(3) = .8;
% % a(3) = .2;
% t10 = (b-a).*rand(1,1) + a;
% t20 = (b-a).*rand(1,1) + a;
T1s = [t10; patTsec; 15]';
T1p = [t20; patTprim; 15]';
Vs = [0; patVsec; 0]'; 
Vp = [0; patVprim; 0]';

t = union(T1s,T1p);

[sharedvals,idxS] = intersect(t,T1s,'stable');
[sharedvals,idxP] = intersect(t,T1p,'stable');

% fitting stuff

    
swt1=ones(1,length(Vp));
idx = find(Vp > 0);
swt1(idx) = 1;
%swt1(idx(end)+1) = .5;
[M1 I1] = max(Vp);
if option == 2
swt1(I1) = 4;
    swt1(I1-1) = 2;
    swt1(I1+1) = 2;
end

swt2 = ones(1,length(Vs));
idx = find(Vs > 0);
%swt2(idx(end)+1) = .5;
swt2(idx) = 1;
[M2 I2] = max(Vs);
if option == 2
swt2(I2) = 3;

    swt2(I2-1) = 1.5;
    swt2(I2+1) = 1.5;
end

S = sort([max(Vp),max(Vs)]);

lb = [1, 1,0.05,S(1),S(1),S(1),S(1),S(1),.05,.05];
ub = [3, 5, 2,S(2),S(2),S(2),S(2),S(2),1,.25];
if rs == 1
    p0 = [1.5,3.5,.1,mean(S),mean(S),mean(S),mean(S),mean(S),.6,.1];
elseif  rs > 1
    p0 = pfits(count-1,1:10);
end
a = .75*p0;
b = 1.25*p0;
% b(3) = .8;
% a(3) = .2;
p0 = (b-a).*rand(1,10) + a;
%p0(3) = .6;
[pfit, resnorm(count2)] = lsqcurvefit(@objfun,p0,t,[Vp.*swt1,Vs.*swt2],lb,ub);
pfitInd(count2,:) = [pfit,t10,t20]; 
tt1 = pfitInd(count2,end-1):.01:20;
tt2 = pfitInd(count2,end):.01:20;
ccP = paramfunPrimCheck(pfitInd(count2,1:10),tt1);
ccS = paramfunSecCheck(pfitInd(count2,1:10),tt2);
MP(count2) = max(ccP(:,1));
MS(count2) = max(ccS(:,1));
count2 = count2+1;
end
if rp == 1 && rs == 1
    save('example_sols.mat','pfitInd')
    data_prim = [patTprim,patVprim];
    data_sec = [patTsec,patVsec];
    save('prmary_data.mat','data_prim')
    save('sec_data.mat','data_sec')
    
end
feasIDX = find(MS < 20 & MP < 20);
if length(feasIDX) < 1
    feasIDX = find(MS < 25 & MP < 25);
end
if length(feasIDX) < 1
    feasIDX = find(MS < 30 & MP < 30);
end

pfitInd = pfitInd(feasIDX,:);
resnorm = resnorm(feasIDX);
MP = MP(feasIDX);
MS = MS(feasIDX);
[m1 I] = min(sqrt((MP-max(Vp)).^2+(MS-max(Vs)).^2));
tt1F = pfitInd(I,end-1):.01:20;
tt2F = pfitInd(I,end):.01:20;
ccPF = paramfunPrimCheck(pfitInd(I,1:10),tt1F);
ccSF = paramfunSecCheck(pfitInd(I,1:10),tt2F);
MPF = max(ccPF(:,1));
MSF = max(ccSF(:,1));

pfits(count,:) = [pfitInd(I,:),MPF,MSF];

check = ceil(.2*length(indDPatprim)*length(indDPatsec));
if mod(count,check) == 0
   fprintf('-----------------------------------\n')
   fprintf('-----------------------------------\n')
   fprintf('percent complete: %d\n',round(100*count./(length(indDPatprim)*length(indDPatsec)),2))
   fprintf('-----------------------------------\n')
   fprintf('-----------------------------------\n')
end

count = count+1;


[rp,rs]
if rp == 1 && rs == 1
    figure
    scatter(patTprim,patVprim)
    xlim([min([min(pfitInd(:,end-1)),min(pfitInd(:,end))]),15])
    figure
    scatter(patTsec,patVsec)
    xlim([min([min(pfitInd(:,end-1)),min(pfitInd(:,end))]),15])
end
end
end
sol = pfits;
end