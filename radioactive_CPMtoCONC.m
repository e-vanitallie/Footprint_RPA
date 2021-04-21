% radioactive_CPMtoCONC.m 
% Elizabeth S Van Itallie - 4/20/2021

% 1. First determine the linear relationship between nanomoles of P32 and
%   measured CPMs using the scintillation counter that is used to measure
%   the probes 
%
% 2. Use the standard curve to determine the number of CPMs that should be
%   in a hybridization reaction of specified volume in order to have a 
%   desired concentration of labeled probe. 


%% Part 1 - Establish the Standard Curve 

% experimental design 
initial_dil = 100;
standard_vols = [0 1 2 5 1 2 5 1 5 10]; % in uL 
standard_dils = [100 100 100 100 10 10 10 1 1 1]; % dilution factors

% calculated parmeters 
dayof_conc = 3.0974; % uM ATP
dayof_SA = 2533.8558; % specific activity

% measured data
measured_cpms = [41.50 40.50 38.01; ...
        1701.44 1702.56  1728.19; ...
        3637.25 3585.98  3554.21; ...
        9127.89 9090.01 9113.66; ...
        17500.61 17444.79 17519.60; ...
        36422.57 36151.32 36104.82; ...
        93246.67 93753.58 93385.31; ...
        189766.8 189490.5 189314.7; ...
        944720.7 944485.6 943910.8; ...
        1890926 1887062 1889304];
    
% ~ Fit a line to find the relationship between the measured CPM values and
% the known P32 nanomoles ~~~~~~~~~~

max_SA = 9120;
frac_hot = dayof_SA/max_SA;
dayof_hot_conc = frac_hot*dayof_conc;
initial_dil_conc = dayof_hot_conc/initial_dil;

num_standards = length(standard_dils);

standard_vals = repmat(initial_dil_conc,1,num_standards)./standard_dils.*standard_vols; % result is in micro-moles
standard_vals_nano = standard_vals*1e3; % result in nano-moles

measured_cpms_T = measured_cpms'; 

standard_vals_MAT = repmat(standard_vals_nano,3,1);
 
P = polyfit(standard_vals_MAT(:,1:end),measured_cpms_T(:,1:end),1);
y_out = polyval(P,standard_vals_nano);
y_out_MAT = repmat(y_out,3,1);

Y_delta_mean = mean(abs(measured_cpms_T - y_out_MAT))./y_out;

figure(1); hold on; set(gcf,'color','white')
p2 = plot(log(reshape(standard_vals_MAT(:,1:end),1,30)),log(reshape(measured_cpms_T(:,1:end),1,30)),'ro','MarkerFaceColor','r');    
p3 = plot(log(standard_vals_nano),log(P(1).*standard_vals_nano+repmat(P(2),1,num_standards)),'b-');
xlabel('Nanomoles P32, Natural Log')
ylabel('CPM measured w/ scintillation counter, Natural Log')

legend([p2,p3],{'Measured Standards (used for fit)',...
    ['CPMs = A * nanomoles P32 + B; A = ' num2str(round(P(1),2)) ', B = ' num2str(round(P(2)),2)]})  

set(gca,'FontSize',20)
title('Standard Curve for Relating CPMs to [P-32]')

%% Part 2 - Use the standard curve to determine the volume of probe to use
% for the TaPR Part 2 hybridization reactions

% ~~~~~~~~~ First determine the average number of UTPs per probe ~~~~~~~~ 

% ~~~~~~ Experimental Variables ~~~~~~~~~~ 

amount_probe_measured = 1; % volume probe added to scintilliation vial in uL 
exp_conc = [7e-6; 7e-6]; % M desired concentration for the final reaction
hybrid_vol = 20; % hybridization reaction volume 
frac_hot_UTP =  0.1757; % fraction P32-UTP/total UTP

% ~~~~ Information About The Probes ~~~~~~  

probe_names = 'Ctnnb1.L 3prime-Sp1, Ctnnb1.S 3prime-Sp1';
probe_f_1 = fastaread('Ctnnb1L-3p1-Sp-PCR-rna-probe.fasta');
probe_f_2 = fastaread('Ctnnb1S-3p1-Sp-PCR-rna-probe.fasta');

% lengths of the probes
probe_lengths = [length(probe_f_1.Sequence); length(probe_f_2.Sequence)];

% find the reverse complement of the rna probe sequence
probe_rc_1 = seqrcomplement(probe_f_1.Sequence);
probe_rc_2 = seqrcomplement(probe_f_2.Sequence);

% find the number of Ts in the reverse complement which will be Us in RNA 
num_Us_1 = length(strfind(probe_rc_1,'T'));
num_Us_2 = length(strfind(probe_rc_2,'T'));
num_Us = [num_Us_1; num_Us_2];

frac_Us = num_Us./probe_lengths;

% ~~~~~~~~~ Then determine the amount of probe in CPM to add to the 
% hybridization reaction for the desired final concentration. 

% how many expected hot UTP per probe
P32_probe_conv = probe_lengths.*frac_Us*frac_hot_UTP;  

% what is the desired concentration of P32
P32_forhyb = (exp_conc*1000*hybrid_vol).*P32_probe_conv;

% P constants come from the standard curve determined above 
forhyb_CPM = P(1)*P32_forhyb+P(2)*ones(size(P32_forhyb)); 

disp(['Add ',num2str(round(forhyb_CPM,0)'), ' CPMs of ',probe_names])
disp(['For a ',num2str(hybrid_vol),'uL hybridization reaction with ',num2str(exp_conc'), 'M of labelled probes.'])

