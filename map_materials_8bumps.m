%% Map material properties according to level set parameters
function m = map_materials_8bumps(ls_pars)
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
alpha = ls_pars(1:nb^2);
beta = ls_pars(nb^2+1:end-1);
%vo = ls_pars(end-1); % background material
vi = ls_pars(end);   % object's material

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
%v = 300*Heaviside(phi) + 2000;
v = vi*Heaviside(phi) + 2000;
m = 1./(v(:)).^2;
return