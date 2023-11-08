resolution = 100;

MachRange = linspace(0.65,0.78,resolution);
AoARange = linspace(0.5,5,resolution);

% aerofoilBehaviour(0.75, AoARange, 1, 1);
% aerofoilBehaviour(MachRange, 3, 1, 1);

% [failAoA, clCPTE, clHTE, clMsh, clCPBound] = aerofoilBehaviour(0.62, AoARange, 1, 0);

AoALimits = zeros(length(MachRange),1);
clLimitsCPTE = zeros(length(MachRange),1);
clLimitsHTE = zeros(length(MachRange),1);
clLimitsMsh = zeros(length(MachRange),1);
clLimitsCPBound = zeros(length(MachRange),1);


machWB = waitbar(0);

for ii=1:1:length(MachRange)
    Mach = MachRange(ii);

    waitbar(ii/length(MachRange), machWB, sprintf('Iteration: %d/%d, M = %0.2f', ii, length(MachRange), Mach));

    [AoALimits(ii), clCPTE, clHTE, clMsh, clCPBound] = aerofoilBehaviour(Mach, AoARange, 1, 0);

    clLimitsCPTE(ii) = clCPTE(2);
    clLimitsHTE(ii) = clHTE(2);
    clLimitsMsh(ii) = clMsh(2);
    clLimitsCPBound(ii) = clCPBound(2);
end

close(machWB);
close all

MLimitsXsh = [0.7445 0.7393 0.7354 0.7314 0.7262 0.7222 0.7183 0.7143 0.7104 0.7275 0.7209 0.7117];
clLimitsXsh = [0.2081 0.2497 0.2914 0.3323 0.3714 0.4111 0.4504 0.4893 0.5279 0.6053 0.6400 0.6669];

figure('Position', [30 30 1000 800]);
plot(MachRange(clLimitsCPBound~=-1), clLimitsCPBound(clLimitsCPBound~=-1), 'DisplayName', 'Mach boundary limit', 'LineWidth', 1.3);
hold on
plot(MachRange(clLimitsMsh~=-1), clLimitsMsh(clLimitsMsh~=-1), 'DisplayName', 'M_{shock} limit', 'LineWidth', 1.3);
plot(MachRange(clLimitsCPTE~=-1), clLimitsCPTE(clLimitsCPTE~=-1), 'DisplayName', 'c_{p(TE)} limit', 'LineWidth', 1.3);
plot(MachRange(clLimitsHTE~=-1), clLimitsHTE(clLimitsHTE~=-1), 'DisplayName', 'H_{TE} limit', 'LineWidth', 1.3);
plot(MLimitsXsh, clLimitsXsh, '-.', 'DisplayName', 'X_{shock} limit', 'LineWidth', 1.3);

scatter(0.73, 0.636, 40, 'o', 'DisplayName', 'Design Operation Point', 'MarkerFaceColor', 'red');

hold off
grid on
xlabel('Mach', 'FontSize', 16);
ylabel('Limit of c_l', 'FontSize', 16);
legend('FontSize', 16, 'Position', [0.23 0.225 0.12 0.05]);