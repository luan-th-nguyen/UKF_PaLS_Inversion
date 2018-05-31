function [xPts, wPts, nPts] = scaledSymmetricSigmaPoints(x,P,lambda)
% This compact function returns the scaled symmetric sigma point distribution, 
% written by Rudolph van der Merwe  

% Number of sigma points and scaling terms
n    = size(x(:),1);
nPts = 2*n+1;            % we're using the symmetric SUT

% Allocate space

wPts=zeros(1,nPts);
xPts=zeros(n,nPts);

% Calculate matrix square root of weighted covariance matrix
[Psqrtm,pp]=chol((n+lambda)*P);  

% Array of the sigma points
xPts=[zeros(size(P,1),1) Psqrtm -Psqrtm];

% Add mean back in
xPts = xPts + repmat(x,1,nPts);  

% Array of the weights for each sigma point
wPts=[lambda 0.5*ones(1,nPts-1)]/(n+lambda);

