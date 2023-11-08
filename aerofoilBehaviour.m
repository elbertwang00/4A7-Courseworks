function [failAoA, clCPTE, clHTE, clMsh, clCPBound] = aerofoilBehaviour(MachRange, AoARange, showPlots, showSecondaryPlots)
   
    global Xupper CPupper Xlower CPlower em cl cm cdv cd1 cd2 HTE Xsh Msh
    
    systemCmd = '/usr/local/lib64/airfoil-hb/vgk < ./Autorun.SIR';
    
    MachRes = length(MachRange);
    AoARes = length(AoARange);
    
    CLsurf = zeros(MachRes, AoARes);
    CDsurf = zeros(MachRes, AoARes);
    CPTEsurf = zeros(MachRes, AoARes);
    HTEsurf = zeros(MachRes, AoARes);
    Xshsurf = zeros(MachRes, AoARes);
    Mshsurf = zeros(MachRes, AoARes);
    cpBoundMinsurf = zeros(MachRes, AoARes);
    
    scriptLines = readlines('.script');
    VINLines = readlines('Autorun.VIN');
    
    loading = waitbar(0, 'aerofoilBehvaiour Waitbar');
    runTotal = MachRes * AoARes;
    runCnt = 0;

    for ii=1:1:MachRes
        Mach = MachRange(ii);
        
        for jj = 1:1:AoARes
            waitbar(runCnt/runTotal, loading, sprintf('Runs completed: %d/%d', runCnt, runTotal));
            incidence = AoARange(jj);
            
            MachStr = sprintf('%0.4f', Mach);
            incStr = sprintf('%0.4f', incidence);
    
            scriptLines(9) = MachStr;
            scriptLines(10) = incStr;
    
            VINLines{3}(15:20) = MachStr;
            VINLines{7}(15:20) = MachStr;
            VINLines{12}(15:20) = MachStr;
            VINLines{3}(25:30) = incStr;
            VINLines{7}(25:30) = incStr;
            VINLines{12}(25:30) = incStr;
            
            fid = fopen('.script', 'w');
            for kk = 1:1:length(scriptLines)
                fprintf(fid, '%s\n', scriptLines(kk));
            end
            fclose(fid);
    
            fid = fopen('Autorun.VIN', 'w');
            for kk = 1:1:length(VINLines)
                fprintf(fid, '%s\n', VINLines(kk));
            end
            fclose(fid);
    
            runCmd(systemCmd);
            
            % Read results of run
            airload Autorun.BRF
    
            machArr = linspace(1.5,1.2,81)';
            cpBuffLineLow = 2*((1+0.2*(machArr-0.05).^2).^(-1.4/0.4)/(1+0.2*em^2)^(-1.4/0.4) - 1) / 1.4 / em^2;
    
            CLsurf(ii,jj) = cl;
            CDsurf(ii,jj) = cd2*10000;
            CPTEsurf(ii,jj) = CPupper(end);
            HTEsurf(ii,jj) = HTE;
            Xshsurf(ii,jj) = Xsh;
            Mshsurf(ii,jj) = Msh;
            cpBoundMinsurf(ii,jj) = -min(CPupper - cpBuffLineLow);

            runCnt = runCnt + 1;
        end
    end

    close(loading);
    
    %% Set up variables for plot-making

    [MachMesh, AoAMesh] = meshgrid(MachRange, AoARange);
    
    if MachRes==1 && AoARes==1
        disp(sprintf('CL: %0.5f, Drag count:%0.2f', CLsurf, CDsurf));
    
    elseif AoARes==1
        constVar = 'AoA';
        independent = 'Mach';
        constVal = AoARange;
    elseif MachRes == 1
        constVar = 'Mach';
        independent = 'AoA';
        constVal = MachRange;
    else
        surf(MachMesh, AoAMesh, CLsurf);
        title('CL')
        axis vis3d
        surf(MachMesh, AoAMesh, CDsurf);
        title('CD')
        axis vis3d
    end
    
    %% Buffet Condition Plots
    figures = findall(groot, 'Type', 'figure');
    % for kk=1:1:length(figures)
    %     if contains(figures(kk).Name, {'Mach', 'AoA'})
    %         close(figures(kk));
    %     end
    % end
    
    independentArr = eval([independent 'Range']);
    
    if showPlots
        figure('Position', [0 0 1980 1080], 'Name', [constVar ' = ' num2str(constVal)]);
        title(sprintf('Variables vs. %s, %s = %0.2f', independent, constVar, MachRange), 'FontSize', 20);
        
        % Distance to Mach Buffet boundary
        ax1 = subplot(2,2,1);
        fill([independentArr flip(independentArr,2)], [5*ones(1,length(independentArr)) zeros(1,length(independentArr))], ...
            'red','FaceAlpha',0.2,'EdgeAlpha',0,'DisplayName', 'Local Mach Condition');
        hold on
        plot(independentArr(cpBoundMinsurf~=-1), cpBoundMinsurf(cpBoundMinsurf~=-1),'Color',[0, 0.4470, 0.7410],'LineWidth',1.3,'HandleVisibility','off');
        yline(0, '--', 'Color', 'red', 'HandleVisibility', 'off');
        hold off
        grid on
        ylim([-1.5, 1]);
        title('(a)');
        xlabel(independent, 'FontSize', 16);
        ylabel('-min(c_{p(Upper)} - c_{p(Local M condition)})', 'FontSize', 16);
        legend;
        
        % Mshock
        ax2 = subplot(2,2,2);
        fill([independentArr flip(independentArr,2)], [5*ones(1,length(independentArr)) 1.2*ones(1,length(independentArr))], ...
            'red','FaceAlpha',0.2,'EdgeAlpha',0,'DisplayName', 'Shock Strength Condition');
        hold on
        plot(independentArr(Mshsurf~=-1), Mshsurf(Mshsurf~=-1),'Color',[0.8500, 0.3250, 0.0980],'LineWidth',1.3,'HandleVisibility','off');
        yline(1.2, '--', 'Color', 'red', 'HandleVisibility', 'off');
        hold off
        grid on
        ylim([0.9 1.4]);
        title('(b)');
        xlabel(independent, 'FontSize', 16);
        ylabel('M_{shock}', 'FontSize', 16);
        legend;
        
        % Cp at trailing edge plot
        ax3 = subplot(2,2,3);
        fill([independentArr flip(independentArr,2)], [5*ones(1,length(independentArr)) (-CPTEsurf(1)+0.04)*ones(1,length(independentArr))], ...
            'red','FaceAlpha',0.2,'EdgeAlpha',0,'DisplayName', 'Trailing Edge c_p Condition');
        hold on
        plot(independentArr(CPTEsurf~=-1), -CPTEsurf(CPTEsurf~=-1),'Color',[0.9290, 0.6940, 0.1250],'LineWidth',1.3,'HandleVisibility','off');
        yline(-CPTEsurf(1)+0.04, '--', 'Color', 'red', 'HandleVisibility', 'off');
        hold off
        grid on
        ylim([-0.25, -0.15]);
        title('(c)');
        xlabel(independent, 'FontSize', 16);
        ylabel('-c_{p(TE)}', 'FontSize', 16);
        legend;
        
        % HTE
        ax4 = subplot(2,2,4);
        fill([independentArr flip(independentArr,2)], [5*ones(1,length(independentArr)) 2.2*ones(1,length(independentArr))], ...
            'red','FaceAlpha',0.2,'EdgeAlpha',0,'DisplayName', 'Trailing Edge H Condition');
        hold on
        plot(independentArr(HTEsurf~=-1), HTEsurf(HTEsurf~=-1),'Color',[0.4940, 0.1840, 0.5560],'LineWidth',1.3,'HandleVisibility','off');
        yline(2.2, '--', 'Color', 'red', 'HandleVisibility', 'off');
        hold off
        grid on
        ylim([1.6, 2.4]);
        title('(d)');
        xlabel(independent, 'FontSize', 16);
        ylabel('H_{TE}', 'FontSize', 16);
        legend;
        
        linkaxes([ax1,ax2,ax3,ax4],'x');
        xlim([min(independentArr), max(independentArr)]);
    end

    %% Calculate outputs for Buffet boundary
    
    failAoA = -1;
    clCPTE = [-1, -1];
    clHTE = [-1, -1];
    clMsh = [-1, -1];
    clCPBound = [-1, -1];

    for ii=1:1:length(independentArr)
        if clCPTE(1) == -1 && CPTEsurf(ii) ~= -1 && -CPTEsurf(ii) > -CPTEsurf(1) + 0.04
            if failAoA == -1    
                failAoA = independentArr(ii);
            end
            clCPTE = [independentArr(ii), CLsurf(ii)];
        end
        if clHTE(1) == -1 && HTEsurf(ii) > 2.2
            if failAoA == -1    
                failAoA = independentArr(ii);
            end
            clHTE = [independentArr(ii), CLsurf(ii)];
        end
        if clMsh(1) == -1 && Mshsurf(ii) > 1.2
            if failAoA == -1    
                failAoA = independentArr(ii);
            end
            clMsh = [independentArr(ii), CLsurf(ii)];
        end
        if clCPBound(1) == -1 && cpBoundMinsurf(ii) > 0
            if failAoA == -1    
                failAoA = independentArr(ii);
            end
            clCPBound = [independentArr(ii), CLsurf(ii)];
        end
    end
    
    %% Secondary Plots
    if showSecondaryPlots
        % Xshock
        figure('Position', [30 30 600 400], 'Name', [constVar ' = ' num2str(constVal)]);
        plot(independentArr(Xshsurf~=-1), Xshsurf(Xshsurf~=-1),'black','LineWidth',1.3);
        grid on
        title('Shock position');
        xlabel(independent, 'FontSize', 16);
        ylabel('X_{shock}', 'FontSize', 16);
    
        % CL
        figure('Position', [30 30 600 400], 'Name', [constVar ' = ' num2str(constVal)]);
        plot(independentArr, CLsurf,'black','LineWidth',1.3);
        grid on
        title('Shock position');
        xlabel(independent, 'FontSize', 16);
        ylabel('c_l', 'FontSize', 16);
    
        % % Drag Count
        % figure('Position', [30 30 600 400]);
        % plot(independentArr(CDsurf~=-1), CDsurf(CDsurf~=-1),'black','LineWidth',1.3);
        % grid on
        % title('Drag Count');
        % xlabel(independent, 'FontSize', 16);
        % ylabel('Drag Count (c_d*1000)', 'FontSize', 16);
        % 
        % % Drag Polar
        % figure('Position', [30 30 600 400]);
        % plot(CLsurf(CDsurf~=-1), CDsurf(CDsurf~=-1),'black','LineWidth',1.3);
        % grid on
        % title('Drag Polar');
        % xlabel('c_l', 'FontSize', 16);
        % ylabel('Drag Count (c_d*1000)', 'FontSize', 16);
    end

    
    function runCmd(cmd)
        [status, cmdout] = system(cmd);
        
        if status ~= 0
            error('Running VGK resulted in an error');
        end
    end
end