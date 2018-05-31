function [U,D] = Helmholtz(m)
% define the grid of equally sized elemnts on the 2D domain
nZ = 101; % number of grid elements in Z (vertical)
nX = 101; % number of grid elements in X (horizontal)
h  = [10 10]; % width of each grid element in Z and X directions respectively
n  = [nZ nX];
n_rs = [51,51]; % [number of sensors, number of receivers]
%z  = [0:nD(1)-1]'*h(1);
%x  = [0:nD(2)-1]'*h(2);

% frequency
f  =  25;  % 25 

% sampling operators
Ps = getP(n,[1:2:101],2);
Pr = getP(n,[1:2:101],100);

% get Helmholtz matrix and source functions
A  = getA(f,m,h,n);
Q  = speye(n_rs(1)); % square domain: n(1) = n(2)

% compute wavefield and data
U  = A\Ps'*Q;
D  = Pr*U;
return