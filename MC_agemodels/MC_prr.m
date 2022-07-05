% Written by Gavin Piccione and Terrence Blackburn, UCSC 2022.

% This code uses a monte carlo simulation to calculate an optimized
% age model for subglacial precipitate PRR50489 within age uncertainty bounds. 
% It then calculates the fit between the age model and a cliamte threshold.

% Outputs: variable fitprob calculated the probability of fit between
% sample record and threshold

% Figure 1 plots Monte Carlo outputs for precipitate record within age
% uncertainty bounds

% Figure 2 plots outputs for EPICA Dome C ice core record within age
% uncertainty bounds

% Figure 3, subplot 1 plots the two best fit curves from precipitate record
% and EPICA record, with threshold infomation overlaid

% Figure 3, subplot 2 plots the binary output describing fit to climate
% threshold 1 =  threshold fits; 0 = does not fit.

% Code takes ~60s to run on "normal" computer

%% Load data
clear 
close all

load('age_model_MC_optimize_PRR.mat')
load('EDC_fine.mat')

% List of variables loaded variables:
  % Ca_chron = precipitate spectra
  % CHRON_age = Bayesian age for precipitate
  % PRR50489_data = precipitate z positions; age model uncertainty bounds
  % EDC = δ18O data from EDC
  % t_old = EDC age model

%% Define climate threshold bounds

% define threshold on EDC
thresh = 0.5 ; 

%define Ca threshold (e.g. above certain value = calcite, below this value= opal)
Ca_thres = 35;

%% Input criteria for Monte Carlo 

% Define Monte Carlo method: 1 = calcuates age model without depth nodes
% 2 = calcuates age model with defined set of depth nodes set by z_inc.
method = 2; %  1 = no_nodes; 2 nodes;

% Number synthetic data sets:
k = 1e3;

% If using method 2;  create new coarse depth vector based on number of nodes. 
z_inc = 2000;

%% step 1 : input age age model + uncertainty bounds

v_l = length(PRR50489_data(:,1));

for i = 1:(v_l-1)
    del_z(i) = PRR50489_data(i+1,1)-PRR50489_data(i,1);
end

% set bottom of sample z to zero. 
PRR50489_data(1,4) = 0;

% use del z to reassign depths for entire sample
for i = 1:v_l-1
   PRR50489_data(i+1,4) =  PRR50489_data(i,4)+ del_z(1);
end

%% step 2: define 1e6 radom accumulation histories within that age model

% step a. Randomly select time at z = 0 (n=1) within permissable time range

% time range: 
a = PRR50489_data(1,2); %ka
b = PRR50489_data(1,3); % ka

for j=1:k
r = a + (b-a) .* rand(1,1);
syn(1,j) = r;
end

% step b. Randomly select accuulation rate between n=i & n=i+1, exclude negative slopes

if method ==1

 for j=1:k
     for v=2:v_l
       a = PRR50489_data(v,2); %ka
       b = PRR50489_data(v,3); % ka  
 r = a + (b-a) .* rand(1,1);
  if r > syn(v-1,j)
      r = syn(v-1,j);
  end
 syn(v,j) = r;
     end
 end
 
z = PRR50489_data(:,4);
 
elseif method ==2
    zdif =  PRR50489_data(end,4)- PRR50489_data(end,1);
    z_incb = zdif/(round(zdif/z_inc));
z = PRR50489_data(1,4):z_incb:PRR50489_data(end,4);
f = round(v_l/length(z));
    
z(end) = PRR50489_data(end,4);

for j=1:k
    n=1;
    for v=2:length(z)
        
      n = n+f;
      a = PRR50489_data(n,2); %ka
      b = PRR50489_data(n,3); % ka  
      
r = a + (b-a) .* rand(1,1);

syn(v,j) = r;

    end
end
end
zfine = PRR50489_data(:,4);

 for j=1:length(syn(1,:))
 syn_fine(:,j) = interp1(z,syn(:,j),zfine);
 end


%% Create a record of precipitate mineralogy using the optimized age model

  for g = 1:length(syn_fine(1,:))
     Ca_mc(:,g) = interp1(CHRONage,Ca_chron,syn_fine(:,g)); % 
  end

% Plot updated precipitate mineralogy versus time
figure()
hold on
for h=1:length(syn(1,:))
 plot(syn_fine(:,h),Ca_chron)
 
 end
xlabel('time (ka)')
 ylabel('Ca normalized')
 hold off

% syn(v,k) ages for v = depths and k = number of simulations
% z = course depths that correspond to syn.

%% step 3:  Calculate optimized EPICA Dome C record over the time interval pf precipitate formation

% Input δ18O and age data for EPICA Dome C ice core
EDC_O = (EDC(:,1));
EDC_t= (t_old(:,1));

% identify the start and end time of syn_fine (precipitate age model). 
 for h=1:length(syn(1,:))
dif_st = abs(syn_fine(end,h)- EDC_t); %calculate differences btwn WDC times and a given age models start time. 
[~,ind_st] = min(dif_st(:)); %index of closest value in WDC to the start time of the given age model
st_log(h) = ind_st;

dif_end = abs(syn_fine(1,h)- EDC_t); %calculate differences btwn WDC times and a given age models END time. 
[~,ind_end] = min(dif_end(:)); %index of closest value in WDC to the END time of the given age model
end_log(h) = ind_end;
 end
 
 for h=1:length(syn(1,:))
     ist = st_log(h);
     ien = end_log(h);
 EDC_new(:,h)= interp1(EDC_t(ist:ien),EDC_O(ist:ien),syn_fine(:,h)); % interpolate WAIS(t) at size of syn.

   EDC_new(isnan(EDC))=0;
  b =   EDC_new(:,h);
nearestfun = @(b) interp1(find(b),b(b~=0),(1:length(b))','nearest','extrap');

  EDC_new(:,h) = 0.5*(nearestfun(b) + flip(nearestfun(flip(b))));
%    for i =  1:length(WDC(:,h))
%       if WDC(i,h) ==0
%            WDC(i,h) = WD_O
%            
%        end
 end
 
 % Plot Monte Carlo outputs for EPICA Dome C record within age uncertainty
 % bounds
 figure()
 hold on
 for h=1:length(syn(1,:))
 plot(syn_fine(:,h),  EDC_new(:,h))
 
 end
 xlabel('time (ka)')
 ylabel('EDC')
 hold off

 %% step 4: Identify fit to threshold by vector points

sum_MC = zeros(1,length(syn_fine(1,:)));
M_c_n = zeros(length(syn_fine(:,1)),length(syn_fine(1,:)));

% Output binary vector describing fit to threshold
% In point in precipitate vecotr fits = 1, else = 0
     for h=1:length(syn_fine(1,:)) %100
 
        for i =1:length(syn_fine(:,1)) %1064

            EDCtest = EDC_new(i,h);
            Ca_test = Ca_mc(i,h);

          if  EDCtest > thresh &&  Ca_test  > Ca_thres
          M_c_n(i,h) = 1;
         
          elseif EDCtest < thresh &&  Ca_test  < Ca_thres
          M_c_n(i,h) = 1;
 
          elseif EDCtest > thresh &&  Ca_test  < Ca_thres
          M_c_n(i,h) = 0;
          
          elseif EDCtest < thresh &&  Ca_test  > Ca_thres
          M_c_n(i,h) = 0;
          end
 
        end
 
       sum_MC(h) = sum(M_c_n(:,h));
 
     end

%% step 5a Identify best fit path
% Best fit based on sum of total fit points

 fit_points = max(sum_MC); % number of points out of 736 that fit threshold values
 fitprob = (fit_points/length(syn_fine(:,1)))*100 % probability of fit to threshold
 best_row = find(sum_MC == fit_points);
 
for h=1:length(syn_fine(1,:))
    sEDC(:,h)= smooth(EDC_new(:,h),3);
end

%% Plot Ca spectra from best fit age model EDC
 figure(3)
subplot(2,1,1) 
 yyaxis right
 plot(syn_fine(:,best_row),Ca_mc(:,best_row), 'r-','DisplayName','MA113')
 yline(Ca_thres, 'r--','HandleVisibility','off')
 hold on
 yyaxis left
 plot(syn_fine(:,best_row),sEDC(:,best_row), 'b-','DisplayName','WDC')
 yline(thresh, 'b--','HandleVisibility','off')
 title('Best Fit MC Age Model')
 xlabel('Time (years)')
 legend

 % Plot fit points
 subplot(2, 1, 2)
 plot(syn_fine(:,best_row),M_c_n(:,best_row))
 title('Fit to Climate Threshold')
 xlabel('Time (years)')
 ylim([-0.01 1.1])

