% radioactive_labeling_P32.m 
% Elizabeth S Van Itallie - 4/20/2021

% (This code requires MATLAB SIMULINK.)
%
% Use this code to determine the amount of P32 labeled nucleotide and 
% unlabeled nucleotide to add to a reaction of specifed volume to achieve
% a specified final concentration of radioacitivity and total (radioactive
% + not) nucleotide concentration taking into account decay since "as of
% day".

% This code can be used for radioactive labeling reactions using
% alpha-UTP to make antisense probes for TaPR and for those using
% gamma-ATP to make ladders for TaPR. 

% Example variables for each that correspond to Chapter 4 of 
% Elizabeth Van Itallie's PhD Thesis are below.

% ~~~~~~~ variable parameters for ladder labeling with ATP- y32P ~~~~~ 

nucl = 'ATP';
T_conc_nucl = 0.1; % uM;  total concentration ATP
frac_hot_use = 0.15; % desired fraction of ATP that is labelleld 

dil_factor = 10; % fold dilution for the radioactive ATP because without 
% dilution the volume would be too small 
T_vol = 10; % reaction volume
conc_cold_nucl = 1; % uM stock concentration of unlabeled ATP


S_sa = 3000; % Ci/mmol starting specific activity in Ci/mmol
S_hot_conc = 10; % mCi/mL starting hot concentration in mCi/mL 

initial_cond = [1 1]; % starting point for optimization

% ~~~~~~~~~ 
% 
% nucl = 'UTP';
% T_conc_nucl = 2.5; % the desired total concentration of UTP is uM 
% frac_hot_use = 0.1757; % the desired final specific activity of hot UTP
% 
% dil_factor = 1;
% T_vol = 10; % total volume of transcription reaction 
% conc_cold_nucl = 30; % concentration of cold UTP added in uM 
% 
% S_sa = 6000; % Ci/mmol starting specific activity in Ci/mmol
% S_hot_conc = 40; % mCi/mL starting hot concentration in mCi/mL 
% 
% initial_cond = [1 1]; % starting point for optimization

% ~~~ Days elapsed sinse "as of data" variable ~~~~~~~

day_elapsed = 11; % number of days elapsed since "as of date" for the 
% radioactive ATP

% ~~~~~~~~~~~ constants for P-32 ~~~~~~~~~~~ 

M_sa = 9120; % Ci/mmol theoretical maximum specific activity in Ci/mmol 
r = 0.0485; %-log(1/2)/14.3  

% ~~~~~~~~~~ Calculate inputs from variables and constants ~~~~~~

S_mc = S_hot_conc/S_sa*1000; % umol/L starting concentration in umol/L (uM)
E_mc = S_mc/dil_factor; 
Hot_frac = S_sa/M_sa; 
Carrier_frac = 1 - Hot_frac;
Hot_conc_start = Hot_frac*E_mc;
Carrier_conc = Carrier_frac*E_mc; % this is 2.28 uM in stock 
Hot_conc_use = Hot_conc_start*exp(-r*day_elapsed); 

% set up the constants for solving the system of linear equations - 
% pass the constants and initial values to fsolve and then get back the
% result 


c = [frac_hot_use Hot_conc_use Carrier_conc conc_cold_nucl T_vol T_conc_nucl]; 
x = fsolve(@(x) myfun(x,c), initial_cond); %eve

x_r = round(x,2,'decimals');

% display the answer 
disp(['Final conditions for labelling reactions: [',nucl,'] = ' num2str(T_conc_nucl) ' uM; Fraction P-32:',nucl,' = ' num2str(frac_hot_use)])
disp(['Volume (uL) of DILUTED by ' num2str(dil_factor) ' fold, hot P-32:',nucl,' to use = ' num2str(x_r(1)) ', Volume (uL) cold ' num2str(conc_cold_nucl) ' uM ',nucl,' to use = ' num2str(x_r(2))])


function F = myfun(x,c)

    F = [c(1)*(x(1)*(c(2)+c(3))+x(2)*c(4))-x(1)*c(2);
        (c(5)*c(6))-(x(1)*(c(2)+c(3))+x(2)*c(4))];
  
end