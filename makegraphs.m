function [] = makegraphs(opt)
if opt == 1
    DENV1=readtable("DENV1.csv");
    DENV1 = table2array(DENV1);
    sample = DENV1;
elseif opt == 2
      DENV2=readtable("DENV2.csv");
    DENV2 = table2array(DENV2);
    sample = DENV2;
elseif opt == 3
       DENV3=readtable("DENV3.csv");
    DENV3 = table2array(DENV3);
    sample = DENV3;
elseif opt == 4
    DENV1=readtable("DENV1.csv");
    DENV2=readtable("DENV2.csv");
    DENV3=readtable("DENV3.csv");
    rand_samp1 = sort(randsample(10000,4000));
    rand_samp2 = sort(randsample(10000,4000));
    rand_samp3 = sort(randsample(10000,4000));
    DENV1 = table2array(DENV1(rand_samp1,:));
    DENV2 = table2array(DENV2(rand_samp2,:));
    DENV3 = table2array(DENV3(rand_samp3,:));
    sample = [DENV1;DENV2;DENV3];
end

set(0,'defaultTextInterpreter','latex');
function rr = paramfunPrim(p, t)
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
end

function rr = paramfunSec(p, t)   
    r = p(1);
     a1 = p(15);   % Viral growth enhancing rate induced by cross-reactiveantibody-virus binding (ADE)
    a2 = (3/2).*a1;   % Cross-reactive antiboceil(.01*length(U))dy-virus killing rate upon cooperativebinding
    d = p(2);   % Specific antibody-virus killing rate upon binding
    f1 = p(3);   % Cross-reactive antibody activation rate
    f2 = (5/4).*p(3);  % Specific antibody activation rate
    k1 = 0;   % Antibody interference competition coeficient
    k2 = 0;  % Antibody interference competition coeficient
   A1 = (1/10).*p(4);    % Saturation coeficients of Hill functions for cross-reactive antibody
    A2 = (1/2).*p(5);  % Saturation coeficients of Hill functions for cross-reactive antibody
    C1 = (1/5).*p(6);   % Saturation coeficients of Hill functions for cross-reactive antibody
    B = (1/20).*p(7);  % Saturation coefficients of Hill functions for specific antibody)
    C2 = (1/10).*p(8);  % Saturation coefficients of Hill functions for specific antibody


    u0 = [p(14),p(11),p(10)]; % initial conditions 
    f = @(t,u) [u(1)*(r + (a1*u(2))/(A1+u(2)) - (a2*u(2).^2)/(A2+u(2).^2) - (d*u(3))/(B+u(3)));
                  (f1*u(1)*u(2))/(C1+u(2)+k1*u(3));
                  (f2*u(1)*u(3))/(C2+k2*u(2)+u(3))];
    [~,rr] = ode45(f,t,u0);
end





prim_starts = unique(sort(sample(:,12)));
sec_starts = unique(sort(sample(:,13)));
starts = union(prim_starts,sec_starts);
t_all = [max(starts)+.05:.05:20,20];
t_all = unique(t_all);
virus_prim = zeros([size(sample,1),length(t_all)]);
cross_prim = zeros([size(sample,1),length(t_all)]);
spec_prim = zeros([size(sample,1),length(t_all)]);
virus_sec = zeros([size(sample,1),length(t_all)]);
cross_sec = zeros([size(sample,1),length(t_all)]);
spec_sec = zeros([size(sample,1),length(t_all)]);
rand_samp = sort(randsample(10000,100));
count = 1;
for j = 1:size(sample,1)
    k_prim = find(sample(j,12) == starts);
    k_sec = find(sample(j,13) == starts);
    if k_prim < length(starts)
    t_prim = [starts(k_prim):.05,starts(end),t_all];
    elseif k_prim == length(starts)
        t_prim = [starts(end),t_all];
    end
    if k_sec < length(starts)
    t_sec = [starts(k_sec):.05:starts(end),t_all];
    elseif k_sec == length(starts)
        t_sec = [starts(end),t_all];
    end
    t_prim = unique(t_prim);
    t_sec = unique(t_sec);
    sol_prim = paramfunPrim(sample(j,:),t_prim);
    sol_sec = paramfunSec(sample(j,:),t_sec);
    if count <= 100
    if j == rand_samp(count) 
        t_prim_r(count,:) = linspace(t_prim(1),t_all(1),20);
        t_sec_r(count,:) = linspace(t_sec(1),t_all(1),20);
        sol_prim_r=  paramfunPrim(sample(j,:),t_prim_r(count,:));
        sol_sec_r= paramfunSec(sample(j,:),t_sec_r(count,:));
        virus_prim_r(count,:) = 10.^abs(sol_prim_r(:,1)');
        virus_sec_r(count,:) = 10.^abs(sol_sec_r(:,1)');
        cross_prim_r(count,:) = 10.^abs(sol_prim_r(:,2)');
        cross_sec_r(count,:) = 10.^abs(sol_sec_r(:,2)');
        spec_prim_r(count,:) = 10.^abs(sol_prim_r(:,3)');
        spec_sec_r(count,:) = 10.^abs(sol_sec_r(:,3)');
        count = count+1;
    end
%      r = p(1);
%      a1 = p(15);   % Viral growth enhancing rate induced by cross-reactiveantibody-virus binding (ADE)
%     a2 = (3/2).*a1;   % Cross-reactive antibody-virus killing rate upon cooperativebinding
%     d = p(2);   % Specific antibody-virus killing rate upon binding
%     f1 = p(3);   % Cross-reactive antibody activation rate
%     f2 = (5/4).*p(3);  % Specific antibody activation rate
%     k1 = 0;   % Antibody interference competition coeficient
%     k2 = 0;  % Antibody interference competition coeficient
%    A1 = (1/10).*p(4);    % Saturation coeficients of Hill functions for cross-reactive antibody
%     A2 = (1/2).*p(5);  % Saturation coeficients of Hill functions for cross-reactive antibody
%     C1 = (1/5).*p(6);   % Saturation coeficients of Hill functions for cross-reactive antibody
%     B = (1/20).*p(7);  % Saturation coefficients of Hill functions for specific antibody
%     C2 = (1/10).*p(8);  % Saturation coefficients of Hill functions for specific antibody
% 
%     u0 = [p(14),p(9),p(10)]; % initial conditions 
%     f = @(t,u) [u(1)*(r + (a1*u(2))/(A1+u(2)) - (a2*u(2).^2)/(A2+u(2).^2) - (d*u(3))/(B+u(3)));
%                   (f1*u(1)*u(2))/(C1+u(2)+k1*u(3));
%                   (f2*u(1)*u(3))/(C2+k2*u(2)+u(3))];

    end
    virus_prim(j,:) = abs(sol_prim(end-length(t_all)+1:end,1)');
    cross_prim(j,:) = abs(sol_prim(end-length(t_all)+1:end,2)');
    spec_prim(j,:) = abs(sol_prim(end-length(t_all)+1:end,3)');
    
    virus_sec(j,:) = abs(sol_sec(end-length(t_all)+1:end,1)');
    cross_sec(j,:) = abs(sol_sec(end-length(t_all)+1:end,2)');
    spec_sec(j,:) =  abs(sol_sec(end-length(t_all)+1:end,3)');

    
net_rep_prim(j,:) = sample(j,1) + (sample(j,15).*cross_prim(j,:))./(sample(j,4).*.1+cross_prim(j,:)) - ((3/2).*sample(j,15)*cross_prim(j,:).^2)./((1/2).*sample(j,5)+cross_prim(j,:).^2);
net_rep_sec(j,:) = sample(j,1) + (sample(j,15).*cross_sec(j,:))./(sample(j,4).*.1+cross_sec(j,:)) - ((3/2).*sample(j,15)*cross_sec(j,:).^2)./((1/2).*sample(j,5)+cross_sec(j,:).^2);
cross_act_prim(j,:) = (sample(j,3).*cross_prim(j,:))./(sample(j,6).*(1/5)+cross_prim(j,:));
cross_act_sec(j,:) = (sample(j,3).*cross_sec(j,:))./(sample(j,6).*(1/5)+cross_sec(j,:));
spec_act_prim(j,:) = ((5/4).*sample(j,3).*spec_prim(j,:))./(sample(j,8).*(1/10)+spec_prim(j,:));
spec_act_sec(j,:) = ((5/4).*sample(j,3).*spec_sec(j,:))./(sample(j,8).*(1/10)+spec_sec(j,:));

% net_cross_prim(j,:) =
% net_cross_sec(j,:) = 
        virus_prim(j,:) = 10.^abs(sol_prim(end-length(t_all)+1:end,1)');
    cross_prim(j,:) = 10.^abs(sol_prim(end-length(t_all)+1:end,2)');
    spec_prim(j,:) = 10.^abs(sol_prim(end-length(t_all)+1:end,3)');
    
    virus_sec(j,:) = 10.^abs(sol_sec(end-length(t_all)+1:end,1)');
    cross_sec(j,:) = 10.^abs(sol_sec(end-length(t_all)+1:end,2)');
    spec_sec(j,:) = 10.^abs(sol_sec(end-length(t_all)+1:end,3)');
        [MV I] = max(sol_sec(:,1));
    max_sec(j) = MV;
    max_sec_time(j) = t_sec(I);
    j
end
% figure
% hold on
% plot(t_all,quantile(virus_prim,.025),'color','red')
% plot(t_all,quantile(virus_prim,1-.025),'color','red')
% plot(t_all,quantile(virus_prim,.5),'color','red')
% for j = 1:20
% plot(t_prim_r(j,:),virus_prim_r(j,:),'r:')
% end
% hold off
% end
if opt < 4
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
if opt == 4
    indDPatprim = find(table2array([data1p(:,5);data2p(:,5);data3p(:,5)]) == 1);
    indDPatsec = find(table2array([data1s(:,5);data2s(:,5);data3s(:,5)]) == 1);
    TTp = [TT1p;TT2p;TT3p];
    VVp = [VV1p;VV2p;VV3p];
    TTs = [TT1s;TT2s;TT3s];
    VVs =[VV1s;VV2s;VV3s];
elseif opt == 1
    indDPatprim = indD1Patprim;
    indDPatsec = indD1Patsec;
    TTp = TT1p;
    TTs = TT1s;
    VVp = VV1p;
    VVs = VV1s;
elseif opt == 2
    indDPatprim = indD2Patprim;
    indDPatsec = indD2Patsec;
    TTp = TT2p;
    TTs = TT2s;
    
    VVp = VV2p;
    VVs = VV2s;
elseif opt == 3
    indDPatprim = indD3Patprim;
    indDPatsec = indD3Patsec;
    TTp = TT3p;
    TTs = TT3s;
    
    VVp = VV3p;
    VVs = VV3s;
end


x = [t_all,fliplr(t_all)];



   fig = figure;
   fig.Position = [1 1 800 800];
  
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
hold on
yyaxis left
set(gca, 'YScale', 'log')
% for r1p = 1:length(indDPatsec)
   
% if r1p < length(indDPatsec)
%     pat1Tsec = TTs(indDPatsec(r1p):indDPatsec(r1p+1)-1);
%     pat1Vsec = VVs(indDPatsec(r1p):indDPatsec(r1p+1)-1);
% elseif r1p == length(indDPatsec)
%     pat1Tsec = TTs(indDPatsec(r1p):end);
%     pat1Vsec = VVs(indDPatsec(r1p):end);
% end
% temp = find(pat1Vsec > 0);
% pat1Vsec = pat1Vsec(pat1Vsec > 0);
% pat1Tsec = pat1Tsec(temp);
% if opt == 1
% plot(pat1Tsec,pat1Vsec,'-','Color',[1 	132 	140]./256,'Marker','d')
% elseif opt == 2
% plot(pat1Tsec,pat1Vsec,'-','Color',[71 	137 	189]./256,'Marker','d')
% elseif opt == 3
% plot(pat1Tsec,pat1Vsec,'-','Color',[233 	87 	87]./256,'Marker','d')
% end
% end

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
  if opt == 1
   red = [17 	48 	66]/256;
   elseif opt== 2
         red = [102 	192 	156]/256;
   elseif opt == 3
          red= [235 	175 	63]/256;
  end

h=plot(pat1Tprim,10.^pat1Vprim,'-','Color',red,'Marker','d');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
% plot(t_all,virus_prim(remove+1,:),'Color',red,'LineStyle','-','Marker','None')
% plot(t_all,virus_prim(end-remove,:),'Color',red,'LineStyle','-','Marker','None')
inBetween = [quantile(virus_prim,.025),fliplr(quantile(virus_prim,1-.025))];
h = fill(x, inBetween,left_color,'LineStyle','none','Marker','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',.4)
plot(t_all,median(virus_prim),'Color','black','linewidth',2,'LineStyle','-','Marker','None')
plot(t_all,quantile(virus_prim,.667),'Color','black','linewidth',2,'LineStyle','--','Marker','None')
% hold off
ylim([1 10.^25])
xlim([min([min(t_prim_r(:,1)),min(t_sec_r(:,1))]),t_all(end)])
h=xline(t_all(1),'k--','linewidth',2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
[mprim Iprim] = min(abs(t_all));
Threshold_prim = virus_prim(:,Iprim);
[msec Isec] = min(abs(t_all));
Threshold_sec = virus_sec(:,Isec);
yline(quantile([Threshold_sec;Threshold_prim],.667),'Color',red,'LineStyle','-.','linewidth',2)
for j = 1:20
h = plot(t_prim_r(j,:),virus_prim_r(j,:),'Color',red,'LineStyle',':','Marker','None','linewidth',2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
ylabel('Virus [cDNA eq./mL]')
yyaxis right
set(gca, 'YScale', 'log')
% plot(t_all,cross_prim(remove+1,:),'Color',blue,'LineStyle','-','Marker','None')
% plot(t_all,cross_prim(end-remove,:),'Color',blue,'LineStyle','-','Marker','None')
inBetween = [quantile(cross_prim,.025),fliplr(quantile(cross_prim,1-.025))];
h = fill(x, inBetween,right_color,'LineStyle','none','Marker','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',.4)
h=plot(t_all,median(cross_prim),'Color','black','linewidth',2,'LineStyle','-','Marker','None');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% plot(t_all,spec_prim(remove+1,:),'Color',blue,'LineStyle','-','Marker','None')
% plot(t_all,spec_prim(end-remove,:),'Color',blue,'LineStyle','-','Marker','None')
for j = 1:20
h = plot(t_prim_r(j,:),cross_prim_r(j,:),'Color','black','LineStyle',':','Marker','None','linewidth',2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
inBetween = [quantile(spec_prim,.025),fliplr(quantile(spec_prim,1-.025))];
h = fill(x, inBetween,right_color,'LineStyle','none','Marker','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',.2)
h=plot(t_all,median(spec_prim),'Color','black','linewidth',2,'LineStyle','-','Marker','None');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
for j = 1:20
h = plot(t_prim_r(j,:),spec_prim_r(j,:),'Color','black','LineStyle',':','Marker','None','linewidth',2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
xlabel('Time since symptom onset (days)')
ylabel('IgG [$\mu g /mL$]')
legend('Median', '66.7%','symptom threshold')
ylim([1 10.^10])
if opt == 1
    title('DENV1 Primary')
elseif opt == 2
    title('DENV2 Primary')
elseif opt == 3
    title('DENV3 Primary')
end
set(gca,'FontSize',24)

if opt == 1
    str2 = 'DENV1PrimTime';
elseif opt == 2
    str2 = 'DENV2PrimTime';
elseif opt == 3
    str2 = 'DENV3PrimTime';
end
  baseFileName = sprintf(str2);
        fname = '~/Documents/MATLAB/DENVPlots';
     saveas(gca, fullfile(fname, baseFileName), 'png');
%      close all

   fig = figure;
   fig.Position = [1 1 800 800];
 
   if opt == 1
   red = [1 	132 	140]./256;
   elseif opt == 2
       red = [71 	137 	189]./256;
   elseif opt == 3
       red = [233 	87 	87]./256;
   end
   blue = [0 0 0];
   left_color = red;
   right_color = blue;
      set(fig,'defaultAxesColorOrder',[left_color; right_color]);
box on
grid on
hold on

yyaxis left
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
   if opt == 1
   red = [1 	132 	140]./256;
   elseif opt == 2
       red = [71 	137 	189]./256;
   elseif opt == 3
       red = [233 	87 	87]./256;
   end

h=plot(pat1Tsec,10.^pat1Vsec,'-','Color',red,'Marker','d');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
% for r1p = 1:length(indDPatprim)
%    
% if r1p < length(indDPatprim)
%     pat1Tprim = TTp(indDPatprim(r1p):indDPatprim(r1p+1)-1);
%     pat1Vprim = VVp(indDPatprim(r1p):indDPatprim(r1p+1)-1);
% elseif r1p == length(indDPatprim)
%     pat1Tprim = TTp(indDPatprim(r1p):end);
%     pat1Vprim = VVp(indDPatprim(r1p):end);
% end
% temp = find(pat1Vprim > 0);
% pat1Vprim = pat1Vprim(pat1Vprim > 0);
% pat1Tprim = pat1Tprim(temp);
% 
% end
% plot(t_all,virus_prim(remove+1,:),'Color',red,'LineStyle','-','Marker','None')
% plot(t_all,virus_prim(end-remove,:),'Color',red,'LineStyle','-','Marker','None')
inBetween = [quantile(virus_sec,.025),fliplr(quantile(virus_sec,1-.025))];
h = fill(x, inBetween,left_color,'LineStyle','none','Marker','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',.4)
plot(t_all,median(virus_sec),'Color','black','linewidth',2,'LineStyle','-','Marker','None')
plot(t_all,quantile(virus_sec,.667),'Color','black','linewidth',2,'LineStyle','--','Marker','None')
[msec Isec] = min(abs(t_all));
Threshold_sec = virus_sec(:,Isec);
yline(quantile([Threshold_sec;Threshold_prim],.667),'Color',red,'LineStyle','-.','linewidth',2)
for j = 1:20
h = plot(t_sec_r(j,:),virus_sec_r(j,:),'Color',red,'LineStyle',':','Marker','None','linewidth',2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
% hold off
ylim([1 10.^25])
xlim([min([min(t_prim_r(:,1)),min(t_sec_r(:,1))]),t_all(end)])
h=xline(t_all(1),'k--','linewidth',2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
ylabel('Virus [cDNA eq./mL]')
legend('Median', '66.7%','symptom threshold')
set(gca, 'YScale', 'log')
yyaxis right
% plot(t_all,cross_prim(remove+1,:),'Color',blue,'LineStyle','-','Marker','None')
% plot(t_all,cross_prim(end-remove,:),'Color',blue,'LineStyle','-','Marker','None')
inBetween = [quantile(cross_sec,.025),fliplr(quantile(cross_sec,1-.025))];
h = fill(x, inBetween,right_color,'LineStyle','none','Marker','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',.4)
h=plot(t_all,median(cross_sec),'Color','black','linewidth',2,'LineStyle','-','Marker','None');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% plot(t_all,spec_prim(remove+1,:),'Color',blue,'LineStyle','-','Marker','None')
% plot(t_all,spec_prim(end-remove,:),'Color',blue,'LineStyle','-','Marker','None')
inBetween = [quantile(spec_sec,.025),fliplr(quantile(spec_sec,1-.025))];
h = fill(x, inBetween,right_color,'LineStyle','none','Marker','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',.2)
for j = 1:20
h = plot(t_sec_r(j,:),cross_sec_r(j,:),'Color','black','LineStyle',':','Marker','None','linewidth',2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
h=plot(t_all,median(spec_sec),'Color','black','linewidth',2,'LineStyle','-','Marker','None');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
for j = 1:20
h = plot(t_sec_r(j,:),spec_sec_r(j,:),'Color','black','LineStyle',':','Marker','None','linewidth',2);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
xlabel('Time since symptom onset (days)')
ylabel('IgG [$\mu g /mL$]')
set(gca, 'YScale', 'log')
ylim([1 10.^10])
if opt == 1
    title('DENV1 Secondary')
elseif opt == 2
    title('DENV2 Secondary')
elseif opt == 3
    title('DENV3 Secondary')
end
set(gca,'FontSize',24)

   if opt == 1
   red = [1 	132 	140]./256;
   elseif opt == 2
       red = [71 	137 	189]./256;
   elseif opt == 3
       red = [233 	87 	87]./256;
   end

if opt == 1
    str2 = 'DENV1SecTime';
elseif opt == 2
    str2 = 'DENV2SecTime';
elseif opt == 3
    str2 = 'DENV3SecTime';
end
  baseFileName = sprintf(str2);
        fname = '~/Documents/MATLAB/DENVPlots';
     saveas(gca, fullfile(fname, baseFileName), 'png');
%      close all
elseif opt == 4
    x = [t_all,fliplr(t_all)];
color_prim = [76/256, 114/256, 176/256];
color_sec = [221/256, 132/256, 82/256];
 fig = figure;
   fig.Position = [1 1 800 800];
hold on
inBetween = [quantile(net_rep_prim,.025),fliplr(quantile(net_rep_prim,1-.025))];
h = fill(x, inBetween,color_prim,'LineStyle','none','Marker','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',.4)
plot(t_all,median(net_rep_prim),'Color',color_prim,'linewidth',2)
inBetween = [quantile(net_rep_sec,.025),fliplr(quantile(net_rep_sec,1-.025))];
h = fill(x, inBetween,color_sec,'LineStyle','none','Marker','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',.4)
plot(t_all,median(net_rep_sec),'Color',color_sec,'linewidth',2)
hold off
xlabel('Time since symptom onset (days)')
ylabel('Net viral replication rate [$r_n$]')
set(gca,'FontSize',24)

 fig = figure;
   fig.Position = [1 1 800 800];
hold on
inBetween = [quantile(cross_act_prim,.025),fliplr(quantile(cross_act_prim,1-.025))];
h = fill(x, inBetween,color_prim,'LineStyle','none','Marker','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',.4)
plot(t_all,median(cross_act_prim),'Color',color_prim,'linewidth',2)
inBetween = [quantile(cross_act_sec,.025),fliplr(quantile(cross_act_sec,1-.025))];
h = fill(x, inBetween,color_sec,'LineStyle','none','Marker','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',.4)
plot(t_all,median(cross_act_sec),'Color',color_sec,'linewidth',2)
hold off
xlabel('Time since symptom onset (days)')
ylabel('Cross-reactive IgG activ. rate [$y_A$]')
set(gca,'FontSize',24)

 fig = figure;
   fig.Position = [1 1 800 800];
hold on
inBetween = [quantile(spec_act_prim,.025),fliplr(quantile(spec_act_prim,1-.025))];
h = fill(x, inBetween,color_prim,'LineStyle','none','Marker','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',.4)
plot(t_all,median(spec_act_prim),'Color',color_prim,'linewidth',2)
inBetween = [quantile(spec_act_sec,.025),fliplr(quantile(spec_act_sec,1-.025))];
h = fill(x, inBetween,color_sec,'LineStyle','none','Marker','none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(h,'facealpha',.4)
plot(t_all,median(spec_act_sec),'Color',color_sec,'linewidth',2)
hold off
xlabel('Time since symptom onset (days)')
ylabel('Spec IgG activ. rate [$z_A$]')
set(gca,'FontSize',24)
end
end
