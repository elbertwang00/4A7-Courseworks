%
%  Script file for transonic airfoil design exercise
%
%
global Xupper CPupper Xlower CPlower em cl cm cdv cd1 cd2 HTE CPuold CPlold Xuold ...
    Xlold CPuold2 CPlold2 Xuold2 Xlold2 Msh;

airload Autorun.BRF

disp(sprintf(['========== CL: %0.5f, CM: %0.5f, CDV: %0.5f\n' ...
              '========== CD1: %0.5f, CD2: %0.5f\n' ...
              '========== %0.2f counts, H_TE: %0.3f, Msh: %0.3f'], cl, cm, cdv, cd1, cd2, cd2*10000, HTE, Msh));

if cd2*10000 < 12 && HTE < 2.2
    disp('=== PASS ===');
else
    disp('=== FAIL ===');
end


% Begin by calculating Cp*
cpstar = 2/1.4/em/em*(((2+0.4*em*em)/2.4)^3.5-1);

% Calculate the M=1.5-1.2 line in cp
machArr = linspace(1.5,1.2,length(Xlower));
cpBuffLine = 2*((1+0.2*machArr.^2).^(-1.4/0.4)/(1+0.2*em^2)^(-1.4/0.4) - 1) / 1.4 / em^2;
cpBuffLineLow = 2*((1+0.2*(machArr-0.05).^2).^(-1.4/0.4)/(1+0.2*em^2)^(-1.4/0.4) - 1) / 1.4 / em^2;
cpBuffLineHigh = 2*((1+0.2*(machArr+0.05).^2).^(-1.4/0.4)/(1+0.2*em^2)^(-1.4/0.4) - 1) / 1.4 / em^2;

% close all;

figure('Position',[30 30 1000 800]);
hold on;

% plot(Xuold,CPuold,'Color','#7E2F8E','DisplayName','Upper-Old')
% plot(Xlold,CPlold,'--','Color','#7E2F8E','DisplayName','Lower-Old');

% Aerofoil cp plots
plot(Xupper,CPupper,'Color','#036b18','DisplayName','Upper Surface c_p','LineWidth',1.3);
plot(Xlower,CPlower,'-.','Color','#036b18','DisplayName','Lower Surface c_p','LineWidth',1.3);

% Cp* line
plot([0,1],[cpstar,cpstar],'--','Color','black','DisplayName','c_p*');

% Buffet conditions
plot(Xupper, cpBuffLine,'--','Color','#960606','DisplayName','Mach Buffet Condition','LineWidth',1.3);
plot(Xupper, cpBuffLineLow,'-.','Color','#cc6002','DisplayName','Mach Buffet Condition - Lower Boundary','LineWidth',1.3);
plot(Xupper, cpBuffLineHigh,'-.','Color','#8c02cc','DisplayName','Mach Buffet Condition - Upper Boundary','LineWidth',1.3);

hold off

xlabel('x/c', 'FontSize', 20);
ylabel('-c_p', 'FontSize', 20);

% axis([0,1,-2.75,2]);
set (gca,'Ydir','reverse');
axis(axis);
grid on;
legend('FontSize', 14, 'Position', [0.5 0.25 0 0]);

Xuold=Xupper;
Xlold=Xlower;
CPuold=CPupper;
CPlold=CPlower;