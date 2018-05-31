%% Heaviside function
function out = Heaviside(x)
eps_phi=0.1;
x=x-eps_phi;
out=zeros(size(x));
I= x>0;
out(I)=1;
return

% function out = Heaviside(x)
% zero_one = 0.5*(1+sign(x));
% out=zeros(size(x));
% I= zero_one>0.5;
% out(I)=1;
% return