% probes_labeling.m 
% Elizabeth S Van Itallie - 4/20/2021

% Use this code to determine the amount of P32 labeled nucleotide and 
% unlabeled nucleotide to add to reaction of specifed volume to achieve
% a specified final concentration of radioacitivity and total (radioactive
% + not) nucleotide concentration taking into account decay since "as of
% day".

% these are the experimental variables 

day_elapsed = 34; % number of days relative to the "as of day" 
T_vol = 10; % total volume of transcription reaction 
conc_cold_utp = 30; % concentration of cold UTP added in uM 

frac_hot_use = 0.1757; % the desired final specific activity of hot UTP
T_conc_utp = 2.5; % the desired total concentration of UTP is uM 

S_sa = 6000; % Ci/mmol starting specific activity in Ci/mmol
S_hot_conc = 40; % mCi/mL starting hot concentration in mCi/mL 


% these are constants 
M_sa = 9120; % Ci/mmol theoretical maximum specific activity in Ci/mmol 
r = 0.0485; % -log(1/2)/14.3 


% these are calculation for the inputs in the equation solver

S_mc = S_hot_conc/S_sa*1000; % umol/L starting concentration in umol/L (uM)
Hot_frac = S_sa/M_sa; 
Carrier_frac = 1 - Hot_frac;
Hot_conc_start = Hot_frac*S_mc;
Carrier_conc = Carrier_frac*S_mc;
Hot_conc_use = Hot_conc_start*exp(-r*day_elapsed); 

% set up the constants for solving the system of linear equations - c
% pass the constants and initial values to fsolve and then get back the
% result 

c = [frac_hot_use Hot_conc_use Carrier_conc conc_cold_utp T_vol T_conc_utp]; 
x = fsolve(@(x) myfun(x,c),[1; 0.5]); %eve

x_r = round(x,2,'decimals');

disp(['[UTP] = ' num2str(T_conc_utp) ' uM; Fraction P-32:UTP = ' num2str(frac_hot_use)])
disp(['Volume (uL) labelled P-32:UTP to use = ' num2str(x_r(1)) ', Volume (uL) cold ' num2str(conc_cold_utp) ' uM UTP to use = ' num2str(x_r(2))])

myfun([x_r],c)

function F = myfun(x,c)

    F = [c(1)*(x(1)*(c(2)+c(3))+x(2)*c(4))-x(1)*c(2);   % this is the fraction hot use
        (c(5)*c(6))-(x(1)*(c(2)+c(3))+x(2)*c(4))];      % this is total concentration utp 

end
 
