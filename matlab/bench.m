% Run script to demonstrate the properties of models with consumption and savings
% as well as polyline class. To be used as starting point for working with this models code.
% Written by Fedor Iskhakov, Australian National University, 2016

close all
clear
clear classes % update classes in the memory
addpath('utils');

fprintf('Hi! This is Matlab version %s running on my laptop\n',version())

%% Nice simulation graphics using retirement model
m5=model_retirement;
m5.ngridm=500;
m5.df=1/(1+m5.r); %flat consumption hopefully
m5.sigma=0.35;
m5.lambda=0.000002; 
m5.nsims=50;
m5.init=[5 20];
tic
m5.solve_dcegm;
t=toc;
fprintf('Retirement model solved with\n %d asset points\n %d periods at \n %.7f lambda  \n %.3f sigma \n in %s\n',m5.ngridm,m5.Tbar,m5.lambda,m5.sigma,ht(t));
% save output
m5.policy.write('output/policy')
m5.value.write('output/value')
fprintf('wrote policy and value function to ascii in output/ . exiting matlab.')
% 
% 
% m5.plot('policy');
% m5.plot('value');
% m5.plot('prob_work');
% m5.sim;
% m5.plot('sim');
% fprintf('Simulation plots for retirement model produced\n')
% 
% 
