function ret_ = getDecay
data = readtable('~/Documents/MATLAB/Dengue/Data_DENV3_anti_NS1_Mar_2016.csv');
patitents = unique(table2array(data(:,2)));

close all
options=odeset('RelTol',1e-12,'AbsTol',1e-12);

    function ret = objfxn(p,t)
        u0 = p(2);
f = @(t,u) -p(1).*u.^2;
[~,ret] = ode45(f,t,u0,options);
ret = ret';
if length(igg) == 2
    ret = [ret(1),ret(end)];
end
    end

    function ret = plotfxn(p,t)
        u0 = p(2);
f = @(t,u) -p(1).*u.^2;
[~,ret] = ode45(f,t,u0,options);
ret = ret';
% if length(igg) == 2
%     ret = [ret(1),ret(end)];
% end
    end

for j = 1:length(patitents)
    idx = find(table2array(data(:,2)) == patitents(j));
    times = table2array(data(idx,7));
    igg = table2array(data(idx,8));
    
    [M I] = max(igg);
%     lb(1) = .95*M;
%     ub(1) = 1.05*M;
%     lb(2) = 0;
%     ub(2) = 40;

    
%          figure
%          hold on
        times = times(I:end)
        igg = igg(I:end);
        igg_ = igg(I:end);
        idx2 = find(igg > 0);
        igg = igg(idx2);
        igg_plot = igg(idx2);
        if  length(igg) >= 2
%         igg = log10(igg);
%         M = log10(M);
        times = times(idx2);
        tt = times;
%         times = log10(times);
%         times = times/365;
        start = times(1);
        times = times - start;
        lb = [0,.9*M];
        ub = [inf,1.1*M];
        p0 = [1,M];
%      igg, objfxn(p0,times);
        tplot = times;
    sqwt = ones(size(igg'));
    sqwt(end) = 3;
%     objfxn(p0,times)
    times
    [pfit(j,:), resnorm] = lsqcurvefit(@objfxn,p0,times,igg',lb,ub);
%     array2table(pfit)

figure
hold on
    plot(tt(1):.01:tt(end)+1,plotfxn(pfit(j,:),times(1):.01:times(end)+1))
    scatter(tt,igg')
    hold off
    ylim([0, ceil(M)])
%     plot(tt,objfxn(pfit(j,:),times))
%     scatter(tt,igg')
%     hold off
%     plot(times_+start_,plotfun(pfit(j,:),times_))
%     figure
%     plot(tt,
    
    
% figure
% hold on
% plot(tt(1):1:tt(end),plotfun(pfit(j,:),tt(1)-tt(1):1:tt(end)-tt(1)))
% scatter(tt,igg_plot,'filled')
% hold off
%     
%     hold off
    end
end
hold off
% xlabel('Time post symptom onset (years)')
% ylabel('IGG level ($\log_{10}(\mu g/mL)$')
% set(gca,'FontSize',18)
decayrate = pfit(:,2);
decayrate = decayrate(decayrate > 0);
%decayrate = table2array(decayrate)
x = linspace(min(decayrate),max(decayrate),500);
[pHat,pCI] = gamfit(decayrate)
% ret_ = pHat;
figure
yyaxis right
plot(x,gampdf(x,pHat(1),pHat(2)),'linewidth',2)
ylabel('pdf')
yyaxis left
histogram(decayrate,20)
xlabel('IGG decay rate')
ylabel('Count')
xlim([0 max(x)])
set(gca,'FontSize',18)

figure
hold on
for j = 1:length(patitents)
    idx = find(table2array(data(:,2)) == patitents(j));
    times = table2array(data(idx,7));
    igg = table2array(data(idx,8));
    
    [M I] = max(igg);
    if I < length(igg) 
        times = times(I:end);
        igg = igg(I:end);
        idx2 = find(igg > 0);
        igg = igg(idx2);
        times = times(idx2);
        times = times./365;
        p0 = [M,.1];
    %[pfit(j,:), resnorm] = lsqcurvefit(@objfxn,p0,times,igg,lb,ub);
    plot(times,igg,'k-^')
    tt(j) = min(times);
    iggs(j) = igg(1);
    end
end
hold off
ret = pfit;
idx = find(ret(:,2) > 0);
ret_ = ret(idx,:);
%ret = pHat;
end