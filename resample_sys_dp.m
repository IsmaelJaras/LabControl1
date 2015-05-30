%% Resampling function
function [xk, wk, idx] = resample_sys_dp(xk, wk)

Ns = length(wk);  % Ns = number of particles

% this is performing latin hypercube sampling on wk
edges = min([0 cumsum(wk)'],1); % protect against accumulated round-off
edges(end) = 1;                 % get the upper edge exact
u1 = rand/Ns;
% this works like the inverse of the empirical distribution and returns
% the interval where the sample is to be found
[~, idx] = histc(u1:1/Ns:1, edges);

xk = xk(idx,:);
wk = 1/Ns*ones(Ns,1);

return;  % bye, bye!!!
