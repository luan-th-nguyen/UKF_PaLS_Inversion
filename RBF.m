%% Radial basis function
function out=RBF(r)
out=(max(0,(1-r)).^4).*(4*r+1);
return