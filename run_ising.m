%clear all

%% default parameters
N=100;
s0=0;
seed=1;
tmax=50*N^2;
tlog=N;

% critical temperature
Tc = 2/log(1+sqrt(2));

if 1
%% run the model
tic;
[M,Mvar, Ac]=ising (N, s0, seed, Tc, tmax, tlog, 10);
toc;

%% compute autocorrelation
figure
plot(Ac)
xlim([1 10000]);
xlabel('time')
ylabel('Autocorrelation')
title('lag=10000')

%% plot
figure
t=0:tlog:tmax;
shadedErrorBar(1:length(M),M,Mvar)
xlabel('Time t')
ylabel('Magnetization M')
ylim([0 0.8])
box off
end

if 0
%% run the sweep
% +/- interval around Tc
interval=.35; % percent
min_T=Tc-(Tc*interval);
max_T=Tc+(Tc*interval);
step_T=.05*Tc;

steps_T=min_T:step_T:max_T;
Ms=zeros(1,length(steps_T));
Mb={};

for i=1:length(steps_T)
    M=ising (N, 0, seed, steps_T(i), tmax, tlog);
    %take the ergodic average of last time series values
    Mb{i}=M;
    Ms(i)=mean(M(end-2*N^2:end));
    Mr(j)=2*std(M(end-2*N^2:end))/N;
    i
end

%% plot temperature vs magnetization
plot(steps_T,Ms)
xlabel('Temperature')
ylabel('Magnetization M')
box off
vline(Tc,'r','Tc')
end