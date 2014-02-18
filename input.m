S=[ 40; 40; 60]
R=[ 1.00000    0.0000      0.0000;
    0.00000    1.0000      0.0000;
    0.00000    0.0000      1.4550]

alat = 7.115;
R *= alat;
coreflag = 1;
atoms = [0.75 0.25 0.00; 0.25 0.75 0.00]*alat;

%# freq. grid for BGW epsilon calculation:
init_frequency = 0.0
delta_frequency = 10.0
delta_frequency_step = 0.5
frequency_low_cutoff = 2.0 
frequency_high_cutoff = 2.1