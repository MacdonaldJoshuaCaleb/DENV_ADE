function ret = multivariateplot(opt,n_samps)
% opt 1,2,3 for DENV1,DENV2,DENV3
% opt1 select return, 1 -> parameters, 2-> trajectories 
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



if opt == 1
load('params1.mat')
params1(:,12) = -params1(:,12);
params1(:,13) = -params1(:,13);
mu = mean(log(params1(:,1:15)));
Sigma = cov(log(params1(:,1:15)));
elseif opt == 2
    load('params2.mat')
params2(:,12) = -params2(:,12);
params2(:,13) = -params2(:,13);
mu = mean(log(params2(:,1:15)));
Sigma = cov(log(params2(:,1:15)));
elseif opt == 3
    
    load('params3.mat')
params3(:,12) = -params3(:,12);
params3(:,13) = -params3(:,13);
mu = mean(log(params3(:,1:15)));
Sigma = cov(log(params3(:,1:15)));
end

for j = 1:n_samps
temp_sample = exp(mvnrnd(mu,Sigma,1));
temp_sample(:,12) = -temp_sample(:,12);
temp_sample(:,13) = -temp_sample(:,13);
t1 = temp_sample(:,12):.05:20;
t2 = temp_sample(:,13):.05:20;
t_sol_p = paramfunPrim(temp_sample, t1);
t_sol_s = paramfunSec(temp_sample, t2);
Mp = max(t_sol_p(:,1));
Ms = max(t_sol_s(:,1));
while Mp > 30 || Ms > 30 || Mp <= 1 || Ms <= 1 || t_sol_p(end,2) < 0.25 || t_sol_s(end,2) < 0.25 || t_sol_s(end,2) > 4 || t_sol_p(end,2) > 4
    temp_sample = exp(mvnrnd(mu,Sigma,1));
    temp_sample(:,12) = -temp_sample(:,12); 
    temp_sample(:,13) = -temp_sample(:,13);
    t1 = temp_sample(:,12):.05:20;
    t2 = temp_sample(:,13):.05:20;
    t_sol_p = paramfunPrim(temp_sample, t1);
    t_sol_s = paramfunSec(temp_sample, t2);
    Mp = max(t_sol_p(:,1));
    Ms = max(t_sol_s(:,1));
end
MaxPrims(j) = Mp;
MSecs(j) = Ms;
Gc0s(j) = t_sol_s(1,2);
sample(j,:) = temp_sample;
end





prim_starts = unique(sort(sample(:,12)));
sec_starts = unique(sort(sample(:,13)));
starts = union(prim_starts,sec_starts);
t_all = [starts',starts(end)+.05:.05:20];
virus_prim = zeros([size(sample,1),length(t_all)]);
cross_prim = zeros([size(sample,1),length(t_all)]);
spec_prim = zeros([size(sample,1),length(t_all)]);
virus_sec = zeros([size(sample,1),length(t_all)]);
cross_sec = zeros([size(sample,1),length(t_all)]);
spec_sec = zeros([size(sample,1),length(t_all)]);
repli_prim = zeros([size(sample,1),length(t_all)]);
for j = 1:length(sample(:,1))
    k_prim = find(sample(j,12) == starts);
    k_sec = find(sample(j,13) == starts);
    t_prim = [starts(k_prim:end)',starts(end)+.05:.1:20];
    t_sec = [starts(k_sec:end)',starts(end)+.05:.1:20];
    
    sol_prim = paramfunPrim(sample(j,:),t_prim);
    sol_sec = paramfunSec(sample(j,:),t_sec);
    
    virus_prim(j,end-length(t_prim)+1:end) = sol_prim(:,1)';
    cross_prim(j,end-length(t_prim)+1:end) = sol_prim(:,2)';
    spec_prim(j,end-length(t_prim)+1:end) = sol_prim(:,3)';
    
    virus_sec(j,end-length(t_sec)+1:end) = sol_sec(:,1)';
    cross_sec(j,end-length(t_sec)+1:end) = sol_sec(:,2)';
    spec_sec(j,end-length(t_sec)+1:end) = sol_sec(:,3)';
    [MV I] = max(sol_sec(:,1));
    max_sec(j) = MV;
    max_sec_time(j) = t_sec(I);
    Q = trapz(t_all,virus_prim(j,:));
    cum_prim(j,:) = Q;
    Q = trapz(t_all,virus_sec(j,:));
    cum_sec(j,:) = Q;
    [MV I] = max(sol_prim(:,1));
    max_prim(j) = MV;
    max_prim_time(j) = t_prim(I);


    id_s = find(abs(sol_prim(:,1)) >= 1);
    if length(id_s) > 0
        id_e = find(abs(sol_prim(:,1)) >= 1);
    Vstart_prim(j) = t_prim(id_s(1))-t_prim(1);
    Vend_prim(j) = t_prim(id_e(end));
    Vdur_prim(j) = Vend_prim(j)-Vstart_prim(j)-t_prim(1);
    
    elseif length(id_s) == 0
        Vstart_prim(j) = NaN;
        Vend_prim(j) = NaN;
        Vdur_prim(j) = 0;
    end

    id_s = find(abs(sol_sec(:,1)) >= 1);
    if length(id_s) > 0
        id_e = find(abs(sol_sec(:,1)) >= 1);
    Vstart_sec(j) = t_sec(id_s(1))-t_sec(1);
    Vend_sec(j) = t_sec(id_e(end));
    Vdur_sec(j) = Vend_sec(j)-Vstart_sec(j)-t_sec(1);
    elseif length(id_s) == 0
        
        Vstart_sec(j) = NaN;
        Vend_sec(j) = NaN;
        Vdur_sec(j) = 0;
    end
    satu_prim(j) = sol_prim(end,2);
    repli_prim(j)

    j
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
if opt == 1
ret1 = [sample,max_prim',max_prim_time',max_sec',max_sec_time',cum_prim,cum_sec,Vstart_prim',Vstart_sec',Vend_prim',Vend_sec',Vdur_prim',Vdur_sec',satu_prim'];
% save('retDENV1.mat','ret1')
% save('t_all1.mat','t_all')
% save('virus_prim1.mat','virus_prim')
% save('virus_sec1.mat',"virus_sec")
% save('spec_prim1.mat','spec_prim')
% save('spec_sec1.mat','spec_sec')
% save('cross_prim1.mat','cross_prim')
% save('cross_sec1.mat','cross_sec')
ret = ret1;
elseif opt == 2
    ret2 = [sample,max_prim',max_prim_time',max_sec',max_sec_time',cum_prim,cum_sec,Vstart_prim',Vstart_sec',Vend_prim',Vend_sec',Vdur_prim',Vdur_sec',satu_prim'];
% save('retDENV2.mat','ret2')
% save('t_all2.mat','t_all')
% save('virus_prim2.mat','virus_prim')
% save('virus_sec2.mat',"virus_sec")
% save('spec_prim2.mat','spec_prim')
% save('spec_sec2.mat','spec_sec')
% save('cross_prim2.mat','cross_prim')
% save('cross_sec2.mat','cross_sec')
ret = ret2;
elseif opt == 3
    ret3 = [sample,max_prim',max_prim_time',max_sec',max_sec_time',cum_prim,cum_sec,Vstart_prim',Vstart_sec',Vend_prim',Vend_sec',Vdur_prim',Vdur_sec',satu_prim'];
% save('retDENV3.mat','ret3')
% save('t_all3.mat','t_all')
% save('virus_prim3.mat','virus_prim')
% save('virus_sec3.mat',"virus_sec")
% save('spec_prim3.mat','spec_prim')
% save('spec_sec3.mat','spec_sec')
% save('cross_prim3.mat','cross_prim')
% save('cross_sec3.mat','cross_sec')
ret = ret3;
end


% x = [t_all,fliplr(t_all)];
%    fig = figure;
%    if opt == 1
%    red = [17 	48 	66]/256;
%    elseif opt == 2
%          red = [102 	192 	156]/256;
%    elseif opt == 3
%           red= [235 	175 	63]/256;
%    elseif opt == 4
%        red = [0 0 0];
%    end
%    blue = [0 0 0];
%    left_color = red;
%    right_color = blue;
%       set(fig,'defaultAxesColorOrder',[left_color; right_color]);
% box on
% grid on
% hold on
% yyaxis left
% % for r1p = 1:length(indDPatsec)
%    
% % if r1p < length(indDPatsec)
% %     pat1Tsec = TTs(indDPatsec(r1p):indDPatsec(r1p+1)-1);
% %     pat1Vsec = VVs(indDPatsec(r1p):indDPatsec(r1p+1)-1);
% % elseif r1p == length(indDPatsec)
% %     pat1Tsec = TTs(indDPatsec(r1p):end);
% %     pat1Vsec = VVs(indDPatsec(r1p):end);
% % end
% % temp = find(pat1Vsec > 0);
% % pat1Vsec = pat1Vsec(pat1Vsec > 0);
% % pat1Tsec = pat1Tsec(temp);
% % if opt == 1
% % plot(pat1Tsec,pat1Vsec,'-','Color',[1 	132 	140]./256,'Marker','d')
% % elseif opt == 2
% % plot(pat1Tsec,pat1Vsec,'-','Color',[71 	137 	189]./256,'Marker','d')
% % elseif opt == 3
% % plot(pat1Tsec,pat1Vsec,'-','Color',[233 	87 	87]./256,'Marker','d')
% % end
% % end
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
%   if opt == 1
%    red = [17 	48 	66]/256;
%    elseif opt== 2
%          red = [102 	192 	156]/256;
%    elseif opt == 3
%           red= [235 	175 	63]/256;
%   end
% 
% h=plot(pat1Tprim,pat1Vprim,'-','Color',red,'Marker','d');
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% end
% % plot(t_all,virus_prim(remove+1,:),'Color',red,'LineStyle','-','Marker','None')
% % plot(t_all,virus_prim(end-remove,:),'Color',red,'LineStyle','-','Marker','None')
% inBetween = [virus_prim(remove+1,:),fliplr(virus_prim(end-remove,:))];
% h = fill(x, inBetween,left_color,'LineStyle','none','Marker','none');
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% set(h,'facealpha',.4)
% plot(t_all,median(virus_prim),'Color','black','linewidth',2,'LineStyle','-','Marker','None')
% plot(t_all,quantile(virus_prim,.667),'Color','black','linewidth',2,'LineStyle','--','Marker','None')
% % hold off
% ylim([0 22])
% xlim([t_all(1),t_all(end)])
% 
% [mprim Iprim] = min(abs(t_all));
% Threshold_prim = virus_prim(:,Iprim);
% [msec Isec] = min(abs(t_all));
% Threshold_sec = virus_sec(:,Isec);
% yline(quantile([Threshold_sec;Threshold_prim],.667),'Color',red,'LineStyle','-.','linewidth',2)
% 
% ylabel('Virus [$\log_{10}$(cDNA eq./mL)]')
% yyaxis right
% % plot(t_all,cross_prim(remove+1,:),'Color',blue,'LineStyle','-','Marker','None')
% % plot(t_all,cross_prim(end-remove,:),'Color',blue,'LineStyle','-','Marker','None')
% inBetween = [cross_prim(remove+1,:),fliplr(cross_prim(end-remove,:))];
% h = fill(x, inBetween,right_color,'LineStyle','none','Marker','none');
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% set(h,'facealpha',.4)
% h=plot(t_all,median(cross_prim),'Color','black','linewidth',2,'LineStyle','-','Marker','None');
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% % plot(t_all,spec_prim(remove+1,:),'Color',blue,'LineStyle','-','Marker','None')
% % plot(t_all,spec_prim(end-remove,:),'Color',blue,'LineStyle','-','Marker','None')
% inBetween = [spec_prim(remove+1,:),fliplr(spec_prim(end-remove,:))];
% h = fill(x, inBetween,right_color,'LineStyle','none','Marker','none');
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% set(h,'facealpha',.2)
% h=plot(t_all,median(spec_prim),'Color','black','linewidth',2,'LineStyle','-','Marker','None');
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% xlabel('Time since symptom onset (days)')
% ylabel('IgG [$\log_{10}(\mu g /mL)$]')
% legend('Median', '66.7%','symptom threshold')
% ylim([0 10])
% if opt == 1
%     title('DENV1 Primary')
% elseif opt == 2
%     title('DENV2 Primary')
% elseif opt == 3
%     title('DENV3 Primary')
% end
% set(gca,'FontSize',18)
% 
% if opt == 1
%     str2 = 'DENV1PrimTime';
% elseif opt == 2
%     str2 = 'DENV2PrimTime';
% elseif opt == 3
%     str2 = 'DENV3PrimTime';
% end
%   baseFileName = sprintf(str2);
%         fname = '~/Documents/MATLAB/DENVPlots';
%      saveas(gca, fullfile(fname, baseFileName), 'png');
%      close all
% 
% 
% 
%    fig = figure;
%    if opt == 1
%    red = [1 	132 	140]./256;
%    elseif opt == 2
%        red = [71 	137 	189]./256;
%    elseif opt == 3
%        red = [233 	87 	87]./256;
%    end
%    blue = [0 0 0];
%    left_color = red;
%    right_color = blue;
%       set(fig,'defaultAxesColorOrder',[left_color; right_color]);
% box on
% grid on
% hold on
% yyaxis left
% for r1p = 1:length(indDPatsec)
%    
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
%    if opt == 1
%    red = [1 	132 	140]./256;
%    elseif opt == 2
%        red = [71 	137 	189]./256;
%    elseif opt == 3
%        red = [233 	87 	87]./256;
%    end
% 
% h=plot(pat1Tsec,pat1Vsec,'-','Color',red,'Marker','d');
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% end
% % for r1p = 1:length(indDPatprim)
% %    
% % if r1p < length(indDPatprim)
% %     pat1Tprim = TTp(indDPatprim(r1p):indDPatprim(r1p+1)-1);
% %     pat1Vprim = VVp(indDPatprim(r1p):indDPatprim(r1p+1)-1);
% % elseif r1p == length(indDPatprim)
% %     pat1Tprim = TTp(indDPatprim(r1p):end);
% %     pat1Vprim = VVp(indDPatprim(r1p):end);
% % end
% % temp = find(pat1Vprim > 0);
% % pat1Vprim = pat1Vprim(pat1Vprim > 0);
% % pat1Tprim = pat1Tprim(temp);
% % 
% % end
% % plot(t_all,virus_prim(remove+1,:),'Color',red,'LineStyle','-','Marker','None')
% % plot(t_all,virus_prim(end-remove,:),'Color',red,'LineStyle','-','Marker','None')
% inBetween = [virus_sec(remove+1,:),fliplr(virus_sec(end-remove,:))];
% h = fill(x, inBetween,left_color,'LineStyle','none','Marker','none');
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% set(h,'facealpha',.4)
% plot(t_all,median(virus_sec),'Color','black','linewidth',2,'LineStyle','-','Marker','None')
% plot(t_all,quantile(virus_sec,.667),'Color','black','linewidth',2,'LineStyle','--','Marker','None')
% [msec Isec] = min(abs(t_all));
% Threshold_sec = virus_sec(:,Isec);
% yline(quantile([Threshold_sec;Threshold_prim],.667),'Color',red,'LineStyle','-.','linewidth',2)
% 
% % hold off
% ylim([0 22])
% xlim([t_all(1),t_all(end)])
% ylabel('Virus [$\log_{10}$(cDNA eq./mL)]')
% legend('Median', '66.7%','symptom threshold')
% yyaxis right
% % plot(t_all,cross_prim(remove+1,:),'Color',blue,'LineStyle','-','Marker','None')
% % plot(t_all,cross_prim(end-remove,:),'Color',blue,'LineStyle','-','Marker','None')
% inBetween = [cross_sec(remove+1,:),fliplr(cross_sec(end-remove,:))];
% h = fill(x, inBetween,right_color,'LineStyle','none','Marker','none');
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% set(h,'facealpha',.4)
% h=plot(t_all,median(cross_sec),'Color','black','linewidth',2,'LineStyle','-','Marker','None');
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% % plot(t_all,spec_prim(remove+1,:),'Color',blue,'LineStyle','-','Marker','None')
% % plot(t_all,spec_prim(end-remove,:),'Color',blue,'LineStyle','-','Marker','None')
% inBetween = [spec_sec(remove+1,:),fliplr(spec_sec(end-remove,:))];
% h = fill(x, inBetween,right_color,'LineStyle','none','Marker','none');
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% set(h,'facealpha',.2)
% h=plot(t_all,median(spec_sec),'Color','black','linewidth',2,'LineStyle','-','Marker','None');
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% xlabel('Time since symptom onset (days)')
% ylabel('IgG [$\log_{10}(\mu g /mL)$]')
% ylim([0 10])
% if opt == 1
%     title('DENV1 Secondary')
% elseif opt == 2
%     title('DENV2 Secondary')
% elseif opt == 3
%     title('DENV3 Secondary')
% end
% set(gca,'FontSize',18)
% 
%    if opt == 1
%    red = [1 	132 	140]./256;
%    elseif opt == 2
%        red = [71 	137 	189]./256;
%    elseif opt == 3
%        red = [233 	87 	87]./256;
%    end
% 
% if opt == 1
%     str2 = 'DENV1SecTime';
% elseif opt == 2
%     str2 = 'DENV2SecTime';
% elseif opt == 3
%     str2 = 'DENV3SecTime';
% end
%   baseFileName = sprintf(str2);
%         fname = '~/Documents/MATLAB/DENVPlots';
%      saveas(gca, fullfile(fname, baseFileName), 'png');
%      close all
% 
% % ret = [sample,max_prim',max_prim_time',max_sec',max_sec_time',cum_prim,cum_sec,Vstart_prim',Vstart_sec',Vend_prim',Vend_sec',Vdur_prim',Vdur_sec',satu_prim'];
end