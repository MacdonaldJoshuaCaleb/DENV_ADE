function ret = makePlots(fits,option)
%fits = SystemFits(3)
clear devs
clear y0s
clear M
clear bin
% close all
clear bins
bins = .1:.1:max(fits(:,9));
for j = 1:length(bins)
if j == 1
idx = find(fits(:,9) <= bins(j));
elseif j==length(bins)
idx = find(fits(:,9) > bins(j));
end
if j > 1 && j < length(bins)
idx = find(fits(:,9) > bins(j-1) & fits(:,9) <= bins(j));
end
M(j) = mean(fits(idx,end));
y0s(j) = mean(fits(idx,9));
devs(j) = std(fits(idx,end));
end

figure
hold on
scatter(fits(:,9),fits(:,end))
for j = 1:length(M)
plot([y0s(j), y0s(j)],[M(j)-2*devs(j),M(j)+2*devs(j)],'r','linewidth',2)
end
plot(y0s,M,'k-^','MarkerFaceColor','k','linewidth',2)
hold off

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

function rr = dudt(p, u)   
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

     
    f = @(u) [u(1).*(r + (a1.*u(2))./(A1+u(2)) - (a2.*u(2).^2)./(A2+u(2).^2) - (d.*u(3))./(B+u(3)));
                  (f1.*u(1).*u(2))./(C1+u(2)+k1.*u(3));
                  (f2.*u(1).*u(3))./(C2+k2.*u(2)+u(3))];
    tempv = f(u);
    rr = tempv(1);
end

ss = size(fits);
tt = linspace(-5,20,200);

for j = 1:ss(1)
    ttsP(j,:) = [linspace(fits(j,11),max([max(fits(:,12)),max(fits(:,11))])+.01,4),linspace(max([max(fits(:,12)),max(fits(:,11))])+.02,20,196)];
    ttsS(j,:) = [linspace(fits(j,12),max([max(fits(:,12)),max(fits(:,11))])+.01,4),linspace(max([max(fits(:,12)),max(fits(:,11))])+.02,20,196)];
    UP = paramfunPrimCheck(fits(j,1:10),ttsP(j,:));
    US = paramfunSecCheck(fits(j,1:10),ttsS(j,:));
    Vp(:,j) = UP(:,1);
    Vs(:,j) = US(:,1);
    GcP(:,j) = UP(:,2);
    GcS(:,j) = US(:,2);
    GsP(:,j) = UP(:,3);
    GsS(:,j) = US(:,3);
 
end


if option ~= 4
Vp = Vp';
VpM = mean(Vp(:,5:end));
VpDev= std(Vp(:,5:end));
Vplb = min(Vp(:,5:end));
Vpub = max(Vp(:,5:end));


GcP = GcP';
GcPM = mean(GcP(:,5:end));
GcPDev= std(GcP(:,5:end));
GcPlb = min(GcP(:,5:end));
GcPub = max(GcP(:,5:end));


GsP = GsP';
GsPM = mean(GsP(:,5:end));
GsPDev= std(GsP(:,5:end));
GsPlb = GsPM - 3*GsPDev;
GsPub = GsPM + 3*GsPDev;
GsPlb = min(GsP(:,5:end));
GsPub = max(GsP(:,5:end));

Vs = Vs';
VsM = mean(Vs(:,5:end));
VsDev= std(Vs(:,5:end));
Vslb = min(Vs(:,5:end));
Vsub = max(Vs(:,5:end));


GcS = GcS';
GcSM = mean(GcS(:,5:end));
GcSDev= std(GcS(:,5:end));
GcSlb = min(GcS(:,5:end));
GcSub = max(GcS(:,5:end));


GsS = GsS';
GsSM = mean(GsS(:,5:end));
GsSDev= std(GsS(:,5:end));
GsSlb = GsSM - 3*GsSDev;
GsSub = GsSM + 3*GsSDev;
GsSlb = min(GsS(:,5:end));
GsSub = max(GsS(:,5:end));
ret = []
elseif option == 4
Vp = Vp';
VpM = mean(Vp(:,5:end));
VpM1 = mean(Vp(1:2232,5:end));
VpM2 = mean(Vp(2233:2232+270,5:end));
VpM3 = mean(Vp(2232+270+1:end,5:end));
VpDev= std(Vp(:,5:end));
Vplb = min(Vp(:,5:end));
Vpub = max(Vp(:,5:end));


GcP = GcP';
GcPM = mean(GcP(:,5:end));
GcPM1 = mean(GcP(1:2232,5:end));
GcPM2 = mean(GcP(2233:2232+270,5:end));
GcPM3 = mean(GcP(2232+270+1:end,5:end));
GcPDev= std(GcP(:,5:end));
GcPlb = min(GcP(:,5:end));
GcPub = max(GcP(:,5:end));


GsP = GsP';
GsPM = mean(GsP(:,5:end));
GsPM1 = mean(GsP(1:2232,5:end));
GsPM2 = mean(GsP(2233:2232+270,5:end));
GsPM3 = mean(GsP(2232+270+1:end,5:end));
GsPDev= std(GsP(:,5:end));
GsPlb = GsPM - 3*GsPDev;
GsPub = GsPM + 3*GsPDev;
GsPlb = min(GsP(:,5:end));
GsPub = max(GsP(:,5:end));

Vs = Vs';
VsM = mean(Vs(:,5:end));
VsDev= std(Vs(:,5:end));
Vslb = min(Vs(:,5:end));
Vsub = max(Vs(:,5:end));


GcS = GcS';
GcSM = mean(GcS(:,5:end));
GcSDev= std(GcS(:,5:end));
GcSlb = min(GcS(:,5:end));
GcSub = max(GcS(:,5:end));


GsS = GsS';
GsSM = mean(GsS(:,5:end));
GsSDev= std(GsS(:,5:end));
GsSlb = GsSM - 3*GsSDev;
GsSub = GsSM + 3*GsSDev;
GsSlb = min(GsS(:,5:end));
GsSub = max(GsS(:,5:end));

for j = 1:ss(1)
cum_viral_prim(j) = trapz(ttsP(j,:),Vp(j,:));
cum_viral_sec(j) = trapz(ttsS(j,:),Vs(j,:));

end
%ret = [cum_viral_prim';cum_viral_sec'];
for k = 1:ss(1)
for j = 1:length(ttsP(1,:))
    dVpdt(k,j) = dudt(fits(j,1:10),[Vp(k,j),GcP(k,j), GsP(k,j)]);
    dVsdt(k,j) = dudt(fits(j,1:10),[Vs(k,j),GcS(k,j), GsS(k,j)]);
end
idx = find(abs(dVpdt(k,:)) >= 0.05);
retP(k,:) = [ttsP(k,idx(1))-ttsP(k,1), ttsP(k,idx(end))-ttsP(k,1)];
idx = find(abs(dVsdt(k,:)) >= 0.05);
retS(k,:) = [ttsS(k,idx(1))-ttsS(k,1), ttsS(k,idx(end))-ttsS(k,1)];
[MP tMP] = max(Vp(k,:));
[MS tMS] = max(Vs(k,:));
retMP(k,:) = ttsP(k,tMP)-ttsP(k,1);
retMS(k,:) = ttsS(k,tMS)-ttsS(k,1);
end 
% retM =[retMP;retMS];
% ret = [retP;retS];
% ret = [ret,retM];

ret = GcP(:,end);
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
if option == 4
    indDPatprim = find(table2array([data1p(:,5);data2p(:,5);data3p(:,5)]) == 1);
    indDPatsec = find(table2array([data1s(:,5);data2s(:,5);data3s(:,5)]) == 1);
    TTp = [TT1p;TT2p;TT3p];
    VVp = [VV1p;VV2p;VV3p];
    TTs = [TT1s;TT2s;TT3s];
    VVs =[VV1s;VV2s;VV3s];
elseif option == 1
    indDPatprim = indD1Patprim;
    indDPatsec = indD1Patsec;
    TTp = TT1p;
    TTs = TT1s;
    VVp = VV1p;
    VVs = VV1s;
elseif option == 2
    indDPatprim = indD2Patprim;
    indDPatsec = indD2Patsec;
    TTp = TT2p;
    TTs = TT2s;
    
    VVp = VV2p;
    VVs = VV2s;
elseif option == 3
    indDPatprim = indD3Patprim;
    indDPatsec = indD3Patsec;
    TTp = TT3p;
    TTs = TT3s;
    
    VVp = VV3p;
    VVs = VV3s;
end



mm = min([min(ttsP(1,:)),min(ttsS(1,:))]);

x = [ttsP(1,5:end),fliplr(ttsP(1,5:end))];

set(0,'defaultTextInterpreter','latex');

   fig = figure;
   if opt == 1
   red = [17 	48 	66]/256;
   elseif opt == 2
         red = [102 	192 	156]/256;
   elseif opt == 3
          red= [235 	175 	63]/256;
   elseif opt == 4
       red = [0 0 0];
   end
   blue = [0 0 0];
   left_color = red;
   right_color = blue;
   set(fig,'defaultAxesColorOrder',[left_color; right_color]);
  box on
  grid on
yyaxis right
hold on
inBetween = [GcPlb,fliplr(GcPub)];
h = fill(x, inBetween,'k','LineStyle','none','Marker','none');
set(h,'facealpha',.5)
if option ~= 4
plot(ttsP(1,5:end),GcPM,'k-','linewidth',2)
elseif option == 4
    plot(ttsP(1,5:end),GcPM1,'r-','linewidth',2)
    plot(ttsP(1,5:end),GcPM2,'b-','linewidth',2)
    plot(ttsP(1,5:end),GcPM3,'g-','linewidth',2)
end

inBetween = [GsPlb,fliplr(GsPub)];
h = fill(x, inBetween,'k','LineStyle','none','Marker','none');
set(h,'facealpha',.2)

if option ~= 4
plot(ttsP(1,5:end),GsPM,'k-','linewidth',2)
elseif option == 4
    plot(ttsS(1,5:end),GsPM1,'r-','linewidth',2)
    plot(ttsS(1,5:end),GsPM2,'b-','linewidth',2)
    plot(ttsS(1,5:end),GsPM3,'g-','linewidth',2)
end
for j = 1:ss(1)
plot(ttsP(j,1:4),GsP(j,1:4),'Marker','none','Color',[right_color,.2],'linewidth',.5)
plot(ttsP(j,1:4),GcP(j,1:4),'Marker','none','Color',[right_color,.5],'linewidth',.5)
end
xline(ttsS(1,5),'k--')
hold off
%ylim([0 5])
ylabel('Antibody concentration')
ylim([0 8])
xlim([mm, 20])
yyaxis left
hold on


inBetween = [Vplb,fliplr(Vpub)];
h = fill(x, inBetween,left_color,'LineStyle','none','Marker','none');
set(h,'facealpha',.4)
for r1p = 1:length(indDPatprim)
   
if r1p < length(indDPatprim)
    pat1Tprim = TTp(indDPatprim(r1p):indDPatprim(r1p+1)-1);
    pat1Vprim = VVp(indDPatprim(r1p):indDPatprim(r1p+1)-1);
elseif r1p == length(indDPatprim)
    pat1Tprim = TTp(indDPatprim(r1p):end);
    pat1Vprim = VVp(indDPatprim(r1p):end);
end
temp = find(pat1Vprim > 0);
pat1Vprim = pat1Vprim(pat1Vprim > 0);
pat1Tprim = pat1Tprim(temp);
plot(pat1Tprim,pat1Vprim,'-','Color',left_color,'Marker','d')
end
if option ~= 4
plot(ttsP(1,5:end),VpM,'Color','k','linestyle','-','Marker','none','linewidth',2)
elseif option == 4
plot(ttsP(1,5:end),VpM1,'r','linestyle','-','Marker','none','linewidth',2)
plot(ttsP(1,5:end),VpM2,'b','linestyle','-','Marker','none','linewidth',2)
plot(ttsP(1,5:end),VpM3,'g','linestyle','-','Marker','none','linewidth',2)
end
for j = 1:ss(1)
plot(ttsP(j,1:4),Vp(j,1:4),'Marker','none','Color',left_color,'linewidth',.5)
end
hold off
ylabel('Viral load')
ylim([0 20])
%xlim([min([min(fits(:,11)),min(fits(:,12))]) 15])
xlabel('Time since symptom onset (days)')
if option == 1
title('DENV1 Primary Infection')
elseif option == 2
    title('DENV2 Primary Infection')
elseif option == 3
    title('DENV3 Primary Infection')
elseif option == 4
    title('Primary Infection')
end
xlim([mm, 20])
 set(gca,'FontSize',18)



 
x = [ttsS(1,5:end),fliplr(ttsS(1,5:end))];

set(0,'defaultTextInterpreter','latex');

   fig = figure;
   if option == 1
   red = [1 	132 	140]./256;
   elseif option == 2
       red = [71 	137 	189]./256;
   elseif option == 3
       red = [233 	87 	87]./256;
   end
   blue = [0 0 0];
   left_color = red;
   right_color = blue;
   set(fig,'defaultAxesColorOrder',[left_color; right_color]);
  box on
  grid on
yyaxis right
hold on
inBetween = [GcSlb,fliplr(GcSub)];
h = fill(x, inBetween,'k','LineStyle','none','Marker','none');
set(h,'facealpha',.5)
plot(ttsS(1,5:end),GcSM,'k-','linewidth',2)

inBetween = [GsSlb,fliplr(GsSub)];
h = fill(x, inBetween,'k','LineStyle','none','Marker','none');
set(h,'facealpha',.2)

plot(ttsS(1,5:end),GsSM,'k-','linewidth',2)
for j = 1:ss(1)
plot(ttsS(j,1:4),GsS(j,1:4),'Marker','none','Color',[right_color,.2],'linewidth',.5,'LineStyle','-')
plot(ttsS(j,1:4),GcS(j,1:4),'Marker','none','Color',[right_color,.5],'linewidth',.5,'LineStyle','-')
end
xline(ttsS(1,5),'k--')
hold off
%ylim([0 5])
ylabel('Antibody concentration')
ylim([0 8])
xlim([mm, 20])
yyaxis left
hold on
for r1p = 1:length(indDPatsec)
   
if r1p < length(indDPatsec)
    pat1Tsec = TTs(indDPatsec(r1p):indDPatsec(r1p+1)-1);
    pat1Vsec = VVs(indDPatsec(r1p):indDPatsec(r1p+1)-1);
elseif r1p == length(indDPatsec)
    pat1Tsec = TTs(indDPatsec(r1p):end);
    pat1Vsec = VVs(indDPatsec(r1p):end);
end
temp = find(pat1Vsec > 0);
pat1Vsec = pat1Vsec(pat1Vsec > 0);
pat1Tsec = pat1Tsec(temp);
plot(pat1Tsec,pat1Vsec,'-','Color',left_color,'Marker','d')
end

inBetween = [Vslb,fliplr(Vsub)];
h = fill(x, inBetween,left_color,'LineStyle','none','Marker','none');
set(h,'facealpha',.4)
for r1p = 1:length(indDPatsec)
   
if r1p < length(indDPatsec)
    pat1Tsec = TTs(indDPatsec(r1p):indDPatsec(r1p+1)-1);
    pat1Vsec = VVs(indDPatsec(r1p):indDPatsec(r1p+1)-1);
elseif r1p == length(indDPatsec)
    pat1Tsec = TTs(indDPatsec(r1p):end);
    pat1Vsec = VVs(indDPatsec(r1p):end);
end
temp = find(pat1Vsec > 0);
pat1Vsec = pat1Vsec(pat1Vsec > 0);
pat1Tsec = pat1Tsec(temp);
plot(pat1Tsec,pat1Vsec,'-','Color',left_color,'Marker','d')
end
plot(ttsS(1,5:end),VsM,'Color','black','linestyle','-','Marker','none','linewidth',2)
for j = 1:ss(1)
plot(ttsS(j,1:4),Vs(j,1:4),'Marker','none','Color',left_color,'linewidth',.5,'LineStyle','-')
end
hold off
ylabel('Viral load')
ylim([0 20])
xlim([mm, 20])
%xlim([min([min(fits(:,11)),min(fits(:,12))]) 15])
xlabel('Time since symptom onset (days)')
if option == 1
title('DENV1 Secondary Infection')
elseif option == 2
    title('DENV2 Secondary Infection')
elseif option == 3
    title('DENV3 Secondary Infection')
end
 set(gca,'FontSize',18)
 

 
 %ret= ttsP;
 
%  if option == 2
%      
%      for m = 1:length(indDPatprim)
%          figure
%      hold on
%      for j = length(indD2Patsec)+1:2*length(indD2Patsec)
%          plot(ttsP(j,:),Vp(j,:),'r')
%      end
%      for r1p = m:m
%    
%     if r1p < length(indDPatprim)
%         pat1Tprim = TTp(indDPatprim(r1p):indDPatprim(r1p+1)-1);
%         pat1Vprim = VVp(indDPatprim(r1p):indDPatprim(r1p+1)-1);
%     elseif r1p == length(indDPatprim)
%         pat1Tprim = TTp(indDPatprim(r1p):end);
%         pat1Vprim = VVp(indDPatprim(r1p):end);
%     end
%      end
%     temp = find(pat1Vprim > 0);
%     pat1Vprim = pat1Vprim(pat1Vprim > 0);
%     pat1Tprim = pat1Tprim(temp);
%     plot(pat1Tprim,pat1Vprim,'-','Color','black','Marker','d')
%      hold off
%      end
%  end
 
 

 
 
 
%      figure
%      hold on
%      plot(ttsS(1,5:end),median(fits(:,1))+(4.*GcSM)./(median(fits(:,4))+GcSM) - ((3/2).*4.*GcSM.^2)./(median(fits(:,5))+GcSM.^2),'linewidth',2)
%      
%      plot(ttsP(1,5:end),median(fits(:,1))+(4.*GcPM)./(median(fits(:,4))+GcPM) - ((3/2).*4.*GcPM.^2)./(median(fits(:,5))+GcPM.^2),'linewidth',2)
%      hold off
%      xlabel('Time since symptom onset (days)')
%      ylabel('Net viral replication rate')
%      legend('Secondary', 'Primary')
%       set(gca,'FontSize',18)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% primary %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    red = [17 	48 	66]/256;
%    elseif option == 2
%          red = [102 	192 	156]/256;
%    elseif option == 3
%           red= [235 	175 	63]/256;
%    elseif option == 4
%        red = [0 0 0];
%    end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SECONDARY %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   if option == 1
%    red = [1 	132 	140]./256;
%    elseif option == 2
%        red = [71 	137 	189]./256;
%    elseif option == 3
%        red = [233 	87 	87]./256;
%    end



figure
box on
grid on
hold on
for r1p = 1:length(indDPatsec)
   
if r1p < length(indDPatsec)
    pat1Tsec = TTs(indDPatsec(r1p):indDPatsec(r1p+1)-1);
    pat1Vsec = VVs(indDPatsec(r1p):indDPatsec(r1p+1)-1);
elseif r1p == length(indDPatsec)
    pat1Tsec = TTs(indDPatsec(r1p):end);
    pat1Vsec = VVs(indDPatsec(r1p):end);
end
temp = find(pat1Vsec > 0);
pat1Vsec = pat1Vsec(pat1Vsec > 0);
pat1Tsec = pat1Tsec(temp);
if option == 1
plot(pat1Tsec,pat1Vsec,'-','Color',[1 	132 	140]./256,'Marker','d')
elseif option == 2
plot(pat1Tsec,pat1Vsec,'-','Color',[71 	137 	189]./256,'Marker','d')
elseif option == 3
plot(pat1Tsec,pat1Vsec,'-','Color',[233 	87 	87]./256,'Marker','d')
end
end
for r1p = 1:length(indDPatprim)
   
if r1p < length(indDPatprim)
    pat1Tprim = TTp(indDPatprim(r1p):indDPatprim(r1p+1)-1);
    pat1Vprim = VVp(indDPatprim(r1p):indDPatprim(r1p+1)-1);
elseif r1p == length(indDPatprim)
    pat1Tprim = TTp(indDPatprim(r1p):end);
    pat1Vprim = VVp(indDPatprim(r1p):end);
end
temp = find(pat1Vprim > 0);
pat1Vprim = pat1Vprim(pat1Vprim > 0);
pat1Tprim = pat1Tprim(temp);
if option == 1
plot(pat1Tprim,pat1Vprim,'-','Color',[17 	48 	66]/256,'Marker','d')
elseif option == 2
plot(pat1Tprim,pat1Vprim,'-','Color', 'k','Marker','d')
elseif option == 3
    plot(pat1Tprim,pat1Vprim,'-','Color','k','Marker','d')
end
end
hold off
ylim([0 20])
xlim([-6 15])
xlabel('Time since symptom onset (days)')
ylabel('$\log_{10}$(cDNA equiv./mL')
if option == 1
    title('DENV1 data')
elseif option == 2
    title('DENV2 data')
elseif option == 3
    title('DENV3 data')
end
set(gca,'FontSize',18)
end