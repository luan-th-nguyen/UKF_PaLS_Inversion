function [xEst,PEst,innov,xSigmaPts]=ukf(xEst,PEst,z,D_b,Q,R)
% UKF for function miminization
% xEst: current estimate
% PEst: current estimation error covariance matrix
% z: measured data
% D_b: data from the fixed background model
% Q, R: covariance matrices for the psuedo state transition noise and  measurement noise
states       = size(xEst(:),1);
observations = size(z(:),1);

% Calculate the sigma points and there corresponding weights using the Scaled Unscented Transformation
%lambda = 1;
lambda = 1e-3;
[xSigmaPts, wSigmaPts, nsp] = scaledSymmetricSigmaPoints(xEst, PEst, lambda); 

% Work out the projected sigma points and their means
xPredSigmaPts = zeros(states,nsp);
zPredSigmaPts = zeros(observations,nsp);
parfor i=1:nsp
    xPredSigmaPts(:,i) = xSigmaPts(:,i);
    mi = map_materials_8bumps(xPredSigmaPts(:,i));
    [~,Di] = helmholtz(mi);
    %yi = real(Di(:)-D_b(:));
    yi = real(Di-D_b);
    zPredSigmaPts(:,i) = yi(:);
end

% Calculate the means 
xPred = zeros(states,1);
zPred = zeros(observations,1);
for i=1:nsp
    xPred = xPred + wSigmaPts(i)*xPredSigmaPts(:,i);
    zPred = zPred + wSigmaPts(i)*zPredSigmaPts(:,i);
end


% Work out the covariances and the cross correlations
PPred = zeros(states,states);
PxzPred = zeros(states,observations);
S = zeros(observations,observations);
for i=1:nsp
    PPred   = PPred + wSigmaPts(i)*(xPredSigmaPts(:,i) - xPred)*(xPredSigmaPts(:,i) - xPred)';    
    PxzPred  = PxzPred + wSigmaPts(i)*(xPredSigmaPts(:,i) - xPred)*(zPredSigmaPts(:,i) - zPred)';
    S = S + wSigmaPts(i)*(zPredSigmaPts(:,i) - zPred)*(zPredSigmaPts(:,i) - zPred)'; 
end
PPred = PPred + Q;
S = S + R;

% Calculate Kalman gain
K  = PxzPred / S;

% Calculate Innovation
innov = z - zPred;

% Update mean
xEst = xPred + K*innov;

% Update covariance
PEst = PPred - K*S*K';