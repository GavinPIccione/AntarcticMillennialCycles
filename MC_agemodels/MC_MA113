% Written by Gavin Piccione and Terrence Blackburn, UCSC 2022.

% This code uses a monte carlo simulation to calculate an optimized
% age model for subglacial precipitate MA113 within age uncertainty bounds. 
% It then calculates the fit between the age model and a cliamte threshold.

% Outputs: variable fitprob calculated the probability of fit between
% sample record and threshold

% Figure 1 plots Monte Carlo outputs for precipitate record within age
% uncertainty bounds

% Figure 2 plots outputs for WAIS Divide ice core record within age
% uncertainty bounds

% Figure 3, subplot 1 plots the two best fit curves from precipitate record
% and EPICA record, with threshold infomation overlaid

% Figure 3, subplot 2 plots the binary output describing fit to climate
% threshold 1 =  threshold fits; 0 = does not fit.

%% Load data
clear all
close all

load('age_model_MC_optimize_MA113.mat')

%% Define climate threshold bounds

% define threshold on WDC
thresh = -39.6 ; 
%define Ca threshold (e.g. above certain value = calcite, below this value= opal)
Ca_thres = -5;

% List of variables loaded variables:
  % Ca_chron = precipitate spectra
  % chron_age = Bayesian age for precipitate
  % MA113_age)model = MA113 age model uncertainty bounds
  % WDC_fine = δ18O data from WDC
  % WDC_age = WDC age model

%% Input criteria for Monte Carlo 

strat_sub = 1e-4; 

% Define Monte Carlo method: 1 = calcuates age model without depth nodes
% 2 = calcuates age model with defined set of depth nodes set by z_inc.
method = 2; %  1 = no_nodes; 2 nodes;

% number synthetic data sets:
k = 1e3;

% If using method 2;  create new coarse depth vector based on number of nodes. default is 50
z_inc = 2000;


%% step 1 : input age age model + bounds

v_l = length(MA113_age_model(:,1));

for i = 1:(v_l-1)
    del_z(i) = MA113_age_model(i+1,1)-MA113_age_model(i,1);
end

% set bottom of sample z to zero. 
MA113_age_model(1,4) = 0;

% use del z to reassign depths for entire sample
for i = 1:v_l-1
    MA113_age_model(i+1,4) =  MA113_age_model(i,4)+ del_z(1);
end

%% step 2: define 1e6 radom accumulation histories within that age model

% step a. Randomly select tim at z = 0 (n=1) within permissable time range

%time range: 
a = MA113_age_model(1,2); %ka
b = MA113_age_model(1,3); % ka

for j=1:k
r = a + (b-a) .* rand(1,1);
syn(1,j) = r;
end

% step b. Randomly select accuulation rate between n=i & n=i+1, exclude negative slopes

if method ==1

 for j=1:k
     for v=2:v_l
       a = MA113_age_model(v,2); %ka
       b = MA113_age_model(v,3); % ka  
 r = a + (b-a) .* rand(1,1);
  if r > syn(v-1,j)
      r = syn(v-1,j);
  end
 syn(v,j) = r;
     end
 end

 
z = MA113_age_model(:,4);
 
elseif method ==2
    zdif =  MA113_age_model(end,4)-MA113_age_model(end,1);
    z_incb = zdif/(round(zdif/z_inc));
z = MA113_age_model(1,4):z_incb:MA113_age_model(end,4);
f = round(v_l/length(z));
    
z(end+1) = MA113_age_model(end,4);

for j=1:k
    n=1;
    for v=2:length(z)
        
      n = n+f;
      a = MA113_age_model(n,2); %ka
      b = MA113_age_model(n,3); % ka  
      
r = a + (b-a) .* rand(1,1);

if r > syn(v-1,j)
    r = syn(v-1,j)-strat_sub;
end

syn(v,j) = r;

    end
end
end
 zfine = MA113_age_model(:,4);

 for j=1:length(syn(1,:))
 syn_fine(:,j) = interp1(z,syn(:,j),zfine);
 end



%% Create a record of precipitate mineralogy using the optimized age model

  for g = 1:length(syn_fine(1,:))
     Ca_mc(:,g) = interp1(chron_age,Ca_chron,syn_fine(:,g)); % 
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

%% step 3:  Calculate optimized WAID Divide Ice core record over the time interval of precipitate formation

% Input δ18O and age data for WAIS Divide ice core
 WD_O = (WD_fine(:,1));
 WD_t= (WD_fine(:,2));
 
 % idetify the start and end time of syn_fine (precipitate age model)
 for h=1:length(syn(1,:))
dif_st = abs(syn_fine(end,h)-WD_t); %calculate differences btwn WDC times and a given age models start time. 
[~,ind_st] = min(dif_st(:)); %index of closest value in WDC to the start time of the given age model
st_log(h) = ind_st;

dif_end = abs(syn_fine(1,h)-WD_t); %calculate differences btwn WDC times and a given age models END time. 
[~,ind_end] = min(dif_end(:)); %index of closest value in WDC to the END time of the given age model
end_log(h) = ind_end;
 end
 

 for h=1:length(syn(1,:))
     ist = st_log(h);
     ien = end_log(h);
 WDC(:,h)= interp1(WD_t(ist:ien),WD_O(ist:ien),syn_fine(:,h)); % interpolate WAIS(t) at size of syn.

  WDC(isnan(WDC))=0;
  b =  WDC(:,h);
nearestfun = @(b) interp1(find(b),b(b~=0),(1:length(b))','nearest','extrap');

 WDC(:,h) = 0.5*(nearestfun(b) + flip(nearestfun(flip(b))));
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
 plot(syn_fine(:,h), WDC(:,h))
 
 end
 xlabel('time (ka)')
 ylabel('WDC')
 hold off

%% step 4: Identify fit to threshold by vector points

% if then

sum_MC = zeros(1,length(syn_fine(1,:)));
M_c_n = zeros(length(syn_fine(:,1)),length(syn_fine(1,:)));

% Output binary vector describing fit to threshold
% In point in precipitate vecotr fits = 1, else = 0
     for h=1:length(syn_fine(1,:)) %100
 
        for i =1:length(syn_fine(:,1)) %1064

            WDCtest = WDC(i,h);
            Ca_test = Ca_mc(i,h);

          if  WDCtest > thresh &&  Ca_test  > Ca_thres
          M_c_n(i,h) = 1;
         
          elseif WDCtest < thresh &&  Ca_test  < Ca_thres
          M_c_n(i,h) = 1;
 
          elseif WDCtest > thresh &&  Ca_test  < Ca_thres
          M_c_n(i,h) = 0;
          
          elseif WDCtest < thresh &&  Ca_test  > Ca_thres
          M_c_n(i,h) = 0;
          end
 
        end
 
       sum_MC(h) = sum(M_c_n(:,h));
 
     end

%% step 5a Identify best fit path. WAIS

 fitprob = max(sum_MC);% number of points out of 1064 that fit threshold values
 percprob = ( fitprob/length(syn_fine(:,1)))*100 % probability of fit to threshold
 best_row = find(sum_MC ==  fitprob);

for h=1:length(syn_fine(1,:))
    sWDC(:,h)= smooth(WDC(:,h),3);
end

%% Plot Ca spectra from best fit age model WAIS
 figure()
 subplot(2,1,1) 
 yyaxis right
 plot(syn_fine(:,best_row),Ca_mc(:,best_row), 'r-','DisplayName','MA113')
 yline(Ca_thres, 'r--','HandleVisibility','off')
 hold on
 yyaxis left
 plot(syn_fine(:,best_row),sWDC(:,best_row), 'b-','DisplayName','WDC')
 yline(thresh, 'b--','HandleVisibility','off')
 title('Best Fit MC Age Model')
 xlabel('Time (years)')
 legend
% 
 %Plot rolling Fit to threshold
 subplot(2, 1, 2)
 plot(syn_fine(:,best_row),M_c_n(:,best_row))
 title('Fit to Climate Threshold')
 xlabel('Time (years)')
 ylim([-0.01 1.01])
