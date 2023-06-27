function ret = exampleplot
close all
set(0,'defaultTextInterpreter','latex');
T = readtable('PfitInd.csv');
prim_data = readtable('primary_example.csv');
prim_data = table2array(prim_data);
sec_data = readtable('secondary_example.csv');
sec_data = table2array(sec_data);
pfitInd = table2array(T);
optimum=[1.14854386016169,3.54590720891863,0.0808580962063910,9.41611228272528,8.45747187382417,8.45463745468353,9.07900359201643,9.43037030143069,0.330035141883313,0.188935644936385,-5.26991253692192,-5.41227278427428,8.50155457357576,9.48102309635976]
function rr = paramfunPrim(p,tt)
        
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
end


function rr = paramfunSec(p, t)   
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

figure 
hold on

for j = 1:100
    color = [1,0,0,.1];
    p = pfitInd(j,:);
    tp = p(end-1):.01:15;
    sol = paramfunPrim(p,tp);
    plot(tp,sol(:,1),'Color',color)
end
for j = 1:100
    color = [0,0,1,.1];
    p = pfitInd(j,:);
    tp = p(end-1):.01:15;
    sol = paramfunPrim(p,tp);
    plot(tp,sol(:,2),'Color',color,'Linestyle','-')
end
for j = 1:100
    color = [0,1,0,.1];
    p = pfitInd(j,:);
    tp = p(end-1):.01:15;
    sol = paramfunPrim(p,tp);
    plot(tp,sol(:,3),'Color',color)
end
xline(-5.2699,'r','linewidth',2)

tp = -5.2699:.01:15
sol = paramfunPrim(optimum,tp);
    plot(tp,sol(:,1),'r-.','linewidth',3)
    plot(tp,sol(:,2),'g-.','linewidth',3)
    plot(tp,sol(:,3),'b-.','linewidth',3)
   scatter(prim_data(:,1),prim_data(:,2),120,'filled','red')
hold off
title('Primary Denv 3 - Paitent 8')
xlabel('Time since symptom onset (days)')
ylabel('Concentration')
set(gca,'Fontsize',20)
xlim([-6, 15])
ylim([0, 10])

figure 
hold on

for j = 1:100
    color = [1,0,0,.1];
    p = pfitInd(j,:);
    tp = p(end):.01:15;
    sol = paramfunSec(p,tp);
    plot(tp,sol(:,1),'Color',color)
end
for j = 1:100
    color = [0,0,1,.1];
    p = pfitInd(j,:);
    tp = p(end):.01:15;
    sol = paramfunSec(p,tp);
    plot(tp,sol(:,2),'Color',color,'Linestyle','-')
end
for j = 1:100
    color = [0,1,0,.1];
    p = pfitInd(j,:);
    tp = p(end):.01:15;
    sol = paramfunSec(p,tp);
    plot(tp,sol(:,3),'Color',color)
end
xline(-5.2699,'r','linewidth',2)

tp = -5.4123:.01:15
sol = paramfunSec(optimum,tp);
    plot(tp,sol(:,1),'r-.','linewidth',3)
    plot(tp,sol(:,2),'g-.','linewidth',3)
    plot(tp,sol(:,3),'b-.','linewidth',3)
  scatter(sec_data(:,1),sec_data(:,2),120,'filled','red')
hold off
title('Secondary Denv 3 - Paitent 31')
xlabel('Time since symptom onset (days)')
ylabel('Concentration')
set(gca,'Fontsize',20)
xlim([-6, 15])
ylim([0 10])


figure
hold on
hist(pfitInd(:,end-1))
xline(-5.2699,'r','linewidth',2)
title('Primary Denv 3 - Paitent 8')
ylabel('Count')
xlabel('Mosquito bite time')
set(gca,'Fontsize',20)
hold off

figure
hold on
hist(pfitInd(:,end))
xline(-5.4123,'r','linewidth',2)
title('Secondary Denv 3 - Paitent 31')
ylabel('Count')
xlabel('Mosquito bite time')
set(gca,'Fontsize',20)
hold off
end