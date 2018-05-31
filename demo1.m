close all;
clearvars;

%% INITIALIZATION:
% ==============
T = 16;   % Number of iterative steps.
nZ = 101; % number of grid elements in Z (vertical)
nX = 101; % number of grid elements in X (horizontal)
nD  = [nZ nX];
h  = [10 10]; % width of each grid element in Z and X directions respectively
z  = [0:nD(1)-1]'*h(1);
x  = [0:nD(2)-1]'*h(2);

%% SYNTHETIC DISTURBANCE SCENARIO:
v_b = 2000*ones(nD);    % background velocity
m_b = 1./(v_b(:)).^2;
[U_b,D_b] = helmholtz(m_b);
dv = zeros(nD);
dv(41:61,41:61) = 300;	% 1 rectangular disturbance
v_true = v_b+dv;
m_true = 1./(v_true(:)).^2;
[U_true,D_true] = helmholtz(m_true);

% plot true model and its outputs
figure(10);imagesc(x,z,reshape(v_true,nD)); colorbar; hold on
title('True model');
axis equal; xlabel('x [m]','fontsize',16);ylabel('z [m]','fontsize',16);
x_s = ones(1,51)*h(1);
z_s = [1:2:101]*h(2);
scatter(x_s,z_s,50,'xr');
x_r = ones(1,51)*100*h(1);
z_r = [1:2:101]*h(2);
scatter(x_r,z_r,50,'vr');
legend('Source','Receiver');
set(gca,'fontsize',16);
print(10,'-dpng',['Helmholtz_model1_true']);

figure(20);imagesc(x,z,reshape(real(U_true(:,26)-U_b(:,26)),nD));
title('Full scattered wavefield');
xlim([0 101]);
axis equal; xlabel('x [m]','fontsize',16);ylabel('z [m]','fontsize',15);
set(gca,'fontsize',16);
print(20,'-dpng',['Helmholtz_model1_true_wavefield']);

y = real(D_true(:)-D_b(:)); % pure scattered measurements
r = length(y); % number of measurements
rms_std = rms(y);
y_randn = 0.03*rms_std.*randn(r,1); %3% rms noise
yTrue = y + y_randn;
y_s50 = real(D_true(:,26)-D_b(:,26)); % the 26th source
yTrue_all = reshape(yTrue,[51 51]);
yTrue_s50 = yTrue_all(:,26);

%% LEVEL-SET CONFIGURATIONS
nb = 8;                         % number of bumps nb x nb
alpha = (10.*(+1).^[1:nb^2])';  % free parameters
beta = 20.0*ones(nb^2,1);       % free parameters
vi = 20.0;                      % object's material
ls_pars = [alpha;beta;vi];
m0 = map_materials_8bumps(ls_pars);
[U0,D0] = helmholtz(m0);
y0 = real(D0(:)-D_b(:));
n = length(ls_pars);            % number of variables 

% upper and lower bounds of the free parameters
alphaL = (-600.*(+1).^[1:nb^2])';
alphaU = (600.*(+1).^[1:nb^2])';
betaL = 0.0*ones(nb^2,1);
betaU = 600.0*ones(nb^2,1);
vL = [1700.0];
vU = [2600.0];
ls_parsL = [alphaL;betaL;vL];
ls_parsU = [alphaU;betaU;vU];

%% UKF CONFIGURATIONS
% covariance matrices
R = (4.0e-1)^2*eye(r);     
P0 = diag(((ls_parsU-ls_parsL)/30).^2);
Q = 0.1^2*P0;

xEst = ls_pars;                 % initial estimate
PEst = P0;          
yEst = y0;

mu_ukf = zeros(n,T);            % to store in arrays
mu_ukf(:,1) = ls_pars;
J_ukf = zeros(1,T);

innov = (yTrue - yEst);
J_ukf(1) = innov'/R*innov;      % to trace data misfit

% plot the initial material map
v0 = sqrt(1.0./m0(:));
m0_bound = map_materials_8bumps(mu_ukf(:,1)+3*sqrt(diag(P0)));
v0_bound = sqrt(1.0./m0_bound(:));
figure(1);imagesc(x,z,reshape(v0,nD),[min(v0(:)) max(v0(:))]); colorbar; hold on;
figure(1); h=imagesc(x,z,reshape(v0_bound,nD)); set(h, 'AlphaData', 0.2);hold off; 
title('Initial');
axis equal; xlabel('x [m]','fontsize',16); ylabel('z [m]','fontsize',16);
set(gca,'fontsize',16);
print(1,'-dpng',['Helmholtz_model1_initial']);

figure(2); scatter(1,J_ukf(1),16,'b','filled'); hold on % plot J_ukf
xlabel('Iteration')
ylabel('Data misfit')
pause(.1)

%% UKF ITERATIONS
for t=2:T+1
    fprintf('UKF : iter = %i / %i \n',t,T);

    [mu_ukf(:,t),PEst,innov,sigPts]=ukf(mu_ukf(:,t-1),PEst,yTrue,D_b,Q,R);
        
    m = map_materials_8bumps(mu_ukf(:,t));
    v = sqrt(1.0./m(:));
    m_bound = map_materials_8bumps(mu_ukf(:,t)+3*sqrt(diag(PEst)));
    v_bound = sqrt(1.0./m_bound(:));
    
    figure(1); imagesc(x,z,reshape(v,nD),[min(v(:)) max(v(:))]); colorbar; hold on;
    figure(1); h=imagesc(x,z,reshape(v_bound,nD)); set(h, 'AlphaData', 0.2); hold off;
    title(['Iter. ',num2str(t-1)])
    axis equal; xlabel('x [m]','fontsize',16); ylabel('z [m]','fontsize',16);
    set(gca,'fontsize',16);
    %print(1,'-dpng',['Helmholtz_model1_reconstructed_material_map_iter_',num2str(t-1)]);
    
    J_ukf(t) = innov'/R*innov;
    figure(2); scatter(t,J_ukf(t),16,'b','filled'); % plot J_ukf
    pause(.1)
end;

%% OUTPUTS
% plot comparative data for one exemplary shot
mu_est = mu_ukf(:,end);
m_est = map_materials_8bumps(mu_est);
[U_est,D_est] = helmholtz(m_est);
yEst_s50 = real(D_est(:,26)-D_b(:,26));
figure(3); plot(z(1:2:101),real(D0(:,26)),'--k'); hold on
plot(z(1:2:101),yTrue_s50,'-.b'); hold on
plot(z(1:2:101),yEst_s50,'-r'); 
xlim([0 1000]);
legend('Initial', '3% RMS synthetic','Reconstructed')
xlabel('z [m]','fontsize',16);ylabel('Real(Magnitude)','fontsize',16);
set(gca,'fontsize',16);

% save plots to png files
print(1,'-dpng',['Helmholtz_model1_reconstructed_material_map']);
print(2,'-dpng',['Helmholtz_model1_reconstructed_objective']);
print(3,'-dpng',['Helmholtz_model1_reconstructed_observables']);

% view full data residuals
yEst = real(D_est(:)-D_b(:));
yEst_all = reshape(yEst,[51 51]);
y0_all = reshape(y0, [51 51]);
fig4 = figure('rend','painters','pos',[10 10 800 400]); subplot(1,2,1); imagesc(full(y0_all - yTrue_all)', [min(min(y0_all - yTrue_all)) max(max(y0_all - yTrue_all))]); hold on;
colormap gray
set(gca,'YDir','normal', 'fontsize', 16)
pbaspect([1 1 1])
xlabel('Receiver number', 'fontsize',16); ylabel('Source number', 'fontsize',16); title('a) h(m_0) - d^{obs}', 'fontsize',16);
subplot(1,2,2); imagesc(full(yEst_all - yTrue_all)',[min(min(y0_all - yTrue_all)) max(max(y0_all - yTrue_all))] ); hold on;
colormap gray
set(gca,'YDir','normal', 'fontsize', 16)
pbaspect([1 1 1])
xlabel('Receiver number', 'fontsize',16); ylabel('Source number', 'fontsize',16); title('b) h(m_{16}) - d^{obs}', 'fontsize',16);
print(fig4,'-dpng',['Helmholtz_model1_reconstructed_observables_all']);
