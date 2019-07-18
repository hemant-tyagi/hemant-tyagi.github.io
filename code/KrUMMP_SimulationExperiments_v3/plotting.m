clear
clc
close all


% Defaults 
width = 4.5;     % Width in inches
height = 4.5;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 15;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize

% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz

% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);

root_folder = 'results/';
filename = ['Results_MC400_C4_GminSep0.05_K2_5_uminmax3_10_L4_noise0_mumax0.01'];
load([root_folder,filename,'.mat']); % load the input

%------------------------
% max_loc_error plots
%------------------------
exp_name = '/max_loc_error/';
mkdir([root_folder, filename], exp_name); % create subfolder to dump plots

    for K=Kmin:Kmax
        for l=1:L
            h = figure;
            set(gca, 'FontSize', fsz, 'LineWidth', alw);
            plot(squeeze(min_sep_t(:,K,l)),squeeze(max_loc_wrap_err(:,K,l)),'.')
            %axis auto
            xlabel('$\triangle_l$','interpreter','latex')
            ylabel('maximum $d_w$','interpreter','latex')
            perc_good = (sum(sign(find(max_loc_wrap_err(:,K,l) <= min_sep_t(:,K,l))))/MC)*100;
            title(['(', num2str(perc_good), ')'])
            %pause
            saveas(h,strcat(root_folder,filename,exp_name,'K_',num2str(K),'_l_',num2str(l),'.png'),'png')
        end
    end
    close all
    
%------------------
% mean loc error
%------------------
exp_name = '/mean_loc_error/';
mkdir([root_folder, filename], exp_name); % create subfolder to dump plots

    for K=Kmin:Kmax
        for l=1:L
            h = figure;
            set(gca, 'FontSize', fsz, 'LineWidth', alw);
            plot(squeeze(min_sep_t(:,K,l)), squeeze(mean_loc_wrap_err(:,K,l)), '.')
            axis auto
            xlabel('$\triangle_l$','interpreter','latex')
            ylabel('mean $d_w$','interpreter','latex') 
            perc_good = (sum(sign(find(mean_loc_wrap_err(:,K,l) <= min_sep_t(:,K,l))))/MC)*100;
            title(['(', num2str(perc_good), ')'])
            %pause
            saveas(h,strcat(root_folder,filename,exp_name,'K_',num2str(K),'_l_',num2str(l),'.png'),'png')
        end
    end
    close all