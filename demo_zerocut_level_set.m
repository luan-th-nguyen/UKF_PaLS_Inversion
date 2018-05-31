clear all; close all;

% define the grid of equally sized elemnts on the 2D domain
nZ = 101; % number of grid elements in Z (vertical)
nX = 101; % number of grid elements in X (horizontal)
n  = [nZ nX];
h  = [10 10]; % width of each grid element in Z and X directions respectively
z  = [0:n(1)-1]'*h(1);
x  = [0:n(2)-1]'*h(2);

% Radial basis parametric level set 
nb = 8; % number of bumps in each direction (the same for Z and X)
nCentersZ = [30,36,42,48,54,60,66,72];
nCentersX = [30,36,42,48,54,60,66,72];
coordsCentersZ = zeros(nb^2,1);
coordsCentersX = zeros(nb^2,1);
for i=1:nb
    for j=1:nb
        coordsCentersZ((i-1)*nb+j) = z(nCentersZ(i));
        coordsCentersX((i-1)*nb+j) = x(nCentersX(j));
    end
end

%% sample model 1
alpha = (1.*(-1).^[1:nb^2])'; % free parameters
beta = 80.0*ones(nb^2,1); % free parameters

%% sample model 2
% beta = 80.0*ones(nb^2,1); % free parameters
% alpha = -1*ones(nb^2,1);
% for i = 3:6
%     for j = 3:6
%         alpha(nb*(i-1) + j) = 1;
%         beta(nb*(i-1) + j) = 80;
%     end
% end

vi = 20.0; % object's material

[Ztemp, Xtemp]=meshgrid(z,x);
Z = Ztemp(:);
X = Xtemp(:);
zMat = repmat(Z,1,nb^2);
xMat = repmat(X,1,nb^2);

centersZ = ones(size(zMat))*diag(coordsCentersZ);
centersX = ones(size(xMat))*diag(coordsCentersX);
nu = 0.01;
r = sqrt(  ((xMat-centersX).^2 + (zMat-centersZ).^2 + nu^2)* diag(1./beta.^2)   );
%r = sqrt(  ((xMat-centersX).^2 + (zMat-centersZ).^2 + nu^2)* diag(1./beta)   );
A = RBF(r);
%phi = A*alpha -.4;
phi = A*alpha;
v = vi*Heaviside(phi) + 2000;
m = 1./(v(:)).^2;

% ploting
figure(1);
s = surf(Xtemp,Ztemp,reshape(phi,n),'FaceAlpha',0.3); hold on
s.EdgeColor = 'none';
contour3(Xtemp,Ztemp,reshape(phi,n),'LevelList',0.01,'Color','red')
xlabel('x' ,'fontsize',16);ylabel('z','fontsize',15);
figure(2);
pcolor(Xtemp,Ztemp,reshape(v,n)); hold on
for i = 1:nb
    for j = 1:nb
        if (sign(alpha(nb*(i-1)+j)) == 1)
            scatter(coordsCentersX(nb*(i-1)+j),coordsCentersZ(nb*(i-1)+j),18,'or','filled')
        else
            scatter(coordsCentersX(nb*(i-1)+j),coordsCentersZ(nb*(i-1)+j),18,'ok','filled')
        end
    end
end
axis equal; xlabel('x','fontsize',16);ylabel('z','fontsize',15);
%print(1,'-dpng',['ls_function2']);
%print(2,'-dpng',['ls_zerolevel2']);