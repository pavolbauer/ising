function [Mout, Mvar, Ac] = ising (N, s0, seed, T, tmax, tlog, lag)
%ISING - compute the 2d square lattice Ising model using Metropolis
% Monte carlo (importance sampling.)
%
% Parameters:
% N - Side length of the lattice
% s0 - initial statate (-1=all -1, 1=all +1, 0 = random)
% seed - initial seed for RN generator
% T - temperature
% tmax - max. amount of timesteps (before the equibrilium)
% tlog - time steps until update of output (0=no log/plot)
% lag - the lag used for autocorrelation, in terms of multiples of tlog 
%       (0=no autocorrelation monitoring)
%
% Output:
% Mout - Order parameter sampled at t=0:tlog:tmax
% Mvar - Variance of the order parameter sampled at t=0:tlog:tmax
% Ac - Autocorrelation time at t=0:tlog:tmax
%
% P.Bauer, May 2016
% parts of this code were inspired by the tutorial found at;
% http://www.physics.usyd.edu.au/~sflammia/Courses/StatMech2014/
% advanced/4/a4.html

%initialize random seed
rng(seed);

%initialize state matrix
if(s0==0)
    S=randi(2,N,N);
    S=S-1;
    S(S==0)=S(S==0)-1;
else
    S=ones(N,N).*s0;
end

% log magnetization & autocorrelation over time
Mout=zeros(1,length(0:tlog:tmax));
Mvar=zeros(1,length(0:tlog:tmax));
Ac=zeros(1,length(0:tlog:tmax));
logi=1;

for  t=0:tmax
   %1) pick random spin
   i1= randi(N,1);
   i2= randi(N,1);
   
   %2) compute the change in Energy in flipping this spin
   dE = 2*S(i1,i2)*neighbsum(S,i1,i2);

   %3) compute flipping probability
   P = exp(-dE/T);
   
   %4) decide if a transition will occur
   if rand <= P
    S(i1,i2) = -1*S(i1,i2); 
   end
    
   %5) log magnetization & autocorrelation
   if tlog && mod(t,tlog)==0
     Mout(logi)=abs(sum(sum(S))/N^2);
     Mvar(logi)=mean(mean(S-Mout(logi)));
     
     if lag && lag<logi %compute autocorrelation at log steps
        mean_ac=mean(Mout(logi-lag:logi));
        Ac(logi-lag)=(Mout(logi-lag)-mean_ac)*(Mout(logi)-mean_ac);
        Ac(logi-lag)=Ac(logi-lag)/var(Mout(logi-lag:logi));
     end
     
     logi=logi+1;
   end     
end

% in case no output is needed
if(~tlog)
    Mout=sum(sum(S))/N;
end    

end

function ns=neighbsum(S,x,y)
% returns the # of neighbors in +1 state
    N=size(S,1);
    ns =S(mod((x-1-1), N)+1, y) + ...
        S(mod((x-1+1), N)+1, y) + ...
        S(x, mod((y-1-1), N)+1) + ...
        S(x, mod((y-1+1), N)+1) ;   
end

