clc 
clear all 
close all
% stream = RandStream('mt19937ar','Seed',28);
% RandStream.setDefaultStream(stream);
 
%--------------------------------
% Simulation related parameters
%--------------------------------
MC = 400; % number of monte carlo runs
Kmin = 2; Kmax = 5; % Run from K = Kmin:Kmax
global_min_sep = 0.05; % Global minimum separation between t_i's across all levels.

%-----------------------------
% Set problem parameters here
%-----------------------------
L = 4; % Number of levels
noise_level = 0.0002; % Standard deviation of gaussian noise (added to real and imaginary parts of fourier samples)
mu_largest = 0.01; % Largest variance of a kernel
umin = 3; umax = 10; % Min and Max amplitudes of spikes 

% Set kernel variances (good choice for L=3: mu(1) = .001; mu(2) = .005; mu(3) = .009)
% mu(1) < mu(2) < ----- < mu(L)
mu(L) = mu_largest;
if(L > 1)
    for l=L-1:-1:1
        mu(l) = mu(l+1)/2;  
    end
end


%--------------------------------
% Set algorithm parameters here
%--------------------------------

%------------------------------------------------------------------------------------------------------
% Sampling window width along a direction. For level l, we sample at s0(l) + i where -m  <= i < m - 1
%------------------------------------------------------------------------------------------------------
m = ceil(1/global_min_sep) + 5; 

%--------------------------------------------------------------------------
% Sampling offsets (for L = 3, a good choice is: s0(1) = 80; s0(2) = 40;
% s0(3) = 0)
%--------------------------------------------------------------------------
s0(L) = 0; % Last window
C = 0.6; % constant used below
eps_base = 1e-2;
if(L > 1)
    for l=L-1:-1:1
        eps = eps_base^2;
        s0(l) = m + ceil((C/(sqrt(2*pi^2*(mu(l+1)^2 - mu(l)^2))))*(sqrt(log (mu(L)/(mu(l)*eps))))); 
    end
end


%---------------------------------------------
% Monte Carlo runs for different values of K
%---------------------------------------------
for mc = 1:MC    
        for K = Kmin:Kmax
            
            %--------------------------------------------------------------
            % Generate source parameters (spike amplitudes and locations) 
            % randomly
            %--------------------------------------------------------------
            sgnseq = sign(randn(L,K)); sgnseq(sgnseq == 0) = 1; % random sign sequence for amplitudes
            u = umin + (umax-umin)*rand(L,K); u = u.*sgnseq; % real-valued random spike amplitudes
            t = zeros(L,K); % variable for storing spike locations for all levels
            for l = 1:L
                flag = 0;
                while(flag == 0)
                    temp_seq = rand(1,K); % sample K locations at random in (0,1)
                    temp1 = sort(temp_seq);
                    temp2 = [temp1-1,temp1,temp1+1];
                    minsep_temp_seq = min(abs(diff(temp2)));
                    if(minsep_temp_seq >= global_min_sep) % if min separation >= global min separation then accept 
                        t(l,:) = temp_seq; 
                        flag = 1;
                    end
                end    
            end
            
            %----------------------------------------------------
            % Variables for storing estimated source parameters 
            %----------------------------------------------------
            test = zeros(L,K);
            uest = zeros(L,K);            

            %--------------------------------
            % The kernel unmixing algorithm
            %--------------------------------
            for l = 1:L
                
                pastFourTot = zeros(2*m,1); % stores the sum of fourier samples corresponding to levels estimated in the past
                svals = s0(l)-m:s0(l)+m-1; svals = svals';
                gbar = (sqrt(2*pi)*mu(l))*exp((-1)*(svals.^2)*2*(pi^2)*(mu(l)^2));
                
                %--------------------------------------------------
                % Obtain noiseless Fourier samples in
                % s0(l) -  m : s0(l) + m - 1 and then add complex 
                % valued noise
                %--------------------------------------------------
                [four_clean,~] = FouryComp(svals,t,u,mu);  
                ext_noise = randn(size(four_clean)) + 1i*randn(size(four_clean));
                four_noise = four_clean + noise_level*ext_noise;
                
                %----------------------------------------------------------
                % Now compute the sum of the estimated Fourier samples of
                % the past levels (1 to l -1) at locations 
                % s0(l) -  m : s0(l) + m - 1
                %----------------------------------------------------------
                if (l > 1)
                   for j=1:l-1 
                        pastFour = FouryCompSingle(svals,test(j,:),uest(j,:),mu(j),K); 
                        pastFourTot = pastFourTot + pastFour;
                   end
                end
                
                %----------------------------------------------------------
                % Subtract the effect of Fourier samples of past estimated
                % levels at locations s0(l) -  m : s0(l) + m - 1, followed
                % by deconvolution step
                %----------------------------------------------------------
                four_deconv = (four_noise - pastFourTot)./gbar;                
                [test(l,:),uest(l,:)] = MPestim(four_deconv,s0(l),m,K);
                                
                %-------------------------------------------------------
                % Plot the reconstructed spikes and the original ones
                %-------------------------------------------------------
                temp1 = max(abs(real(uest(l,:)))); % take real part of uest here since that is what we will plot below
                temp2 = max(abs(u(l,:)));
                temp = max(temp1, temp2);
                stem(test(l,:),real(uest(l,:)),'-Xr','MarkerSize', 8, 'LineWidth',2 );  % Plot test_i and real part of uest_i
                hold on
                stem(t(l,:) , u(l,:) ,'-Ob','MarkerSize', 8, 'LineWidth',2 );  % Plot t_i and real part of u_i
                hold off
                xlim([-0.2 1.2]) % Set x-axis limits
                ylim([(-1)*temp-1 temp + 1]) % Set y-axis limits
                title(strcat('l=',num2str(l),', K=',num2str(K),', run=',num2str(mc)))       
                
                %------------------------------
                % Compute statistics for plots
                %------------------------------
                t_locs_wrapp_err = zeros(K,1); % Stores error for t_i at index i
                t_locs_match_inds = zeros(K,1); % index i stores the index of test which matches with t_i
                rel_ampl_errs = zeros(K,1); % ith index stores relative error in amplitude for t_i w.r.t above matching
                
                [t_locs_wrapp_err,t_locs_match_inds] = MatchLocations(test(l,:)',t(l,:)'); % Compute matching
                
                max_loc_wrap_err(mc,K,l) = max(t_locs_wrapp_err); % compute max error in spike locations
                mean_loc_wrap_err(mc,K,l) = mean(t_locs_wrapp_err); % compute mean error in spike locations           
                rel_ampl_errs = abs(u(l,t_locs_match_inds) - uest(l,:))./(abs(u(l,:))); % Compute rel amplitude errors for the spikes
                max_ampl_err(mc,K,l) = max(rel_ampl_errs); % Compute the maximum rel amplitude error
                mean_ampl_err(mc,K,l) = mean(rel_ampl_errs); % Compute the average rel amplitude error
                
                %-------------------------------------------
                % Compute minimum separation between t_i's
                %------------------------------------------
                temp1 = sort(t(l,:));
                temp2 = [temp1-1,temp1,temp1+1];
                min_sep_t(mc,K,l) = min(abs(diff(temp2)));
                
                %---------------------------------------------
                % Display the computed statistics on screen
                %---------------------------------------------
                disp( ['mc = ', num2str(mc), ', l = ', num2str(l), ' min sep = ', num2str(min_sep_t(mc,K,l)), ...  
                       ', dw_max = ' , num2str(max_loc_wrap_err(mc,K,l)), ', dw_mean = ', num2str(mean_loc_wrap_err(mc,K,l)), ...
                       ', amp_err_max = ', num2str(max_ampl_err(mc,K,l)), ', amp_err_mean = ', num2str(mean_ampl_err(mc,K,l))] );
            end
        end
end

%-----------------
% Save workspace 
%-----------------
fileToSave = [ 'results/Results_', 'MC', num2str(MC), '_C',num2str(C),'_GminSep', num2str(global_min_sep), ... 
               '_K', num2str(Kmin) ,'_', num2str(Kmax), '_uminmax', num2str(umin) ,'_', num2str(umax), ...
               '_L', num2str(L), '_noise', num2str(noise_level), ...
               '_mumax', num2str(mu_largest), '.mat' ];
save(fileToSave);        




