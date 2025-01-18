%%solve problem Keller-Segel system in 3D
%Dec-MC-BCFD method
%=============================================================
clc,clear

%%%%%%%%%%%%%%%%some input
M = 40;
hr = 1/M;
kr = hr;
lr = hr;
xL = 0;
xR = 1;
yL = 0;
yR = 1;
zL = 0;
zR = 1;
T = 6e-05;
non=1;

lambda = 1;
alpha = 1;
%%%%%%%%%%%%%%%temporal grid information
tau = 1e-06;
t = 0 : tau : T;    

%%%%%%%%%%%%%%%%spatial grid information
x = xL : hr : xR;
y = yL : kr : yR;
z = zL : lr : zR;
N = length(t);
Nx = length(x) - 1;
Ny = length(y) - 1;
Nz = length(z) - 1;

[x,y,z,hf,kf,lf,h,k,l,xf,yf,zf] = Grid_mid_3D(Nx,Ny,Nz,x,y,z,non);

%%%%%%%%%%%%%Initial condtion
[c_old,n_old] = Initial_solution_3D(xf,yf,zf);
c00 = c_old;
n00 = n_old;

c_old = c_old(:);
n_old = n_old(:);
%%%%%%%%%%%%%%%laplace operator aproximation
[Ddx,dDx] = Laplace(h,hf,Nx);
[Ddy,dDy] = Laplace(k,kf,Ny);
[Ddz,dDz] = Laplace(l,lf,Nz);
%%%%%%%%%%%%%%gradient operator aproximation
Dx = Gradient(h,Nx);
Dy = Gradient(k,Ny);
Dz = Gradient(l,Nz);
dx = Gradient(hf,Nx-1);
dy = Gradient(kf,Ny-1);
dz = Gradient(lf,Nz-1);
%%%%%%%%%%%%Interpolation operator
Interp_x = Interpolation(h,hf,Nx-1);
Interp_y = Interpolation(k,kf,Ny-1);
Interp_z = Interpolation(l,lf,Nz-1);
%%%%%%%%%%%%%%%%%%%%%%%%%
dx_n  = zeros(Nx,Nx-1);
dy_n  = zeros(Ny,Ny-1);
dz_n  = zeros(Nz,Nz-1);
dx_n(2:Nx-1,:)  = dx;
dy_n(2:Ny-1,:)  = dy;
dz_n(2:Nz-1,:)  = dz;
dx_n(1,1) = 1/hf(1);
dx_n(Nx,Nx-1) = -1/hf(Nx);
dy_n(1,1) = 1/kf(1);
dy_n(Ny,Ny-1) = -1/kf(Ny);
dz_n(1,1) = 1/lf(1);
dz_n(Nz,Nz-1) = -1/lf(Nz);
%%%%%%%%%%%%%%%%%%%
Intn_x = Interpolation(h,hf,Nx);
Intn_y = Interpolation(k,kf,Ny);
Intn_z = Interpolation(l,lf,Nz);

Iy = eye(Ny);
Ix = eye(Nx);
Iz = eye(Nz);
%%%%%%%%%Step 1: compute n %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LHS_n = 1/tau*kron( Iz,kron( Iy,Ix ) ) ...
        - kron( Iz,kron( Iy,dx_n*Dx ) )/2 - kron( Iz,kron( dy_n*Dy,Ix ) )/2 - kron( dz_n*Dz,kron( Iy,Ix ) )/2 ...
         + lambda*( kron( Iz,kron( Iy,dx_n ) )*( (kron( Iz,kron( Iy,Dx ) )*c_old).*kron( Iz,kron( Iy,Intn_x ) ) )...
                    + kron( Iz,kron( dy_n,Ix ) )*( (kron( Iz,kron( Dy,Ix ) )*c_old).*kron( Iz,kron( Intn_y,Ix ) ) ) ...
                    + kron( dz_n,kron( Iy,Ix ) )*( (kron( Dz,kron( Iy,Ix ) )*c_old).*kron( Intn_z,kron( Iy,Ix ) ) ));
RHS_n = 1/tau*n_old +  ( kron( Iz,kron( Iy,dx_n*Dx ) ) + kron( Iz,kron( dy_n*Dy,Ix ) ) + kron( dz_n*Dz,kron( Iy,Ix ) ) )*n_old/2;
n_new = LHS_n\RHS_n;
%%%%%%%%%Step 2: compute c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LHS_c = (1/tau+1/2)*kron( Iz,kron( Iy,Ix ) ) ... 
        - ( kron( Iz,kron( Iy,dx_n*Dx ) ) + kron( Iz,kron( dy_n*Dy,Ix ) ) + kron( dz_n*Dz,kron( Iy,Ix ) ) )/2;
RHS_c = (1/tau-1/2)*c_old  + ( kron( Iz,kron( Iy,dx_n*Dx ) ) + kron( Iz,kron( dy_n*Dy,Ix ) ) + kron( dz_n*Dz,kron( Iy,Ix ) ) )*c_old/2 + (n_new+n_old)/2;
c_new = LHS_c\RHS_c;
%%%%%%%%%%%%%%%%%%%%generate variables%%%%%%%%%%%%%%%%%%%%%%%
c_old = c_new;
n_oold = n_old;
n_old = n_new;

for m = 3:N
    Time = t(m);

%%%%%%%%%Step 1: compute c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LHS_c = (1/tau+alpha/2)*kron( Iz,kron( Iy,Ix ) ) ... 
        - ( kron( Iz,kron( Iy,dx_n*Dx ) ) + kron( Iz,kron( dy_n*Dy,Ix ) ) + kron( dz_n*Dz,kron( Iy,Ix ) ) )/2;
RHS_c = (1/tau-alpha/2)*c_old  + ( kron( Iz,kron( Iy,dx_n*Dx ) ) + kron( Iz,kron( dy_n*Dy,Ix ) ) + kron( dz_n*Dz,kron( Iy,Ix ) ) )*c_old/2 ...
        + (3/2*n_old - 1/2*n_oold);
c_new = LHS_c\RHS_c;
%%%%%%%%%Step 2: compute n %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LHS_n = 1/tau*kron( Iz,kron( Iy,Ix ) ) ...
         - kron( Iz,kron( Iy,dx_n*Dx ) )/2 - kron( Iz,kron( dy_n*Dy,Ix ) )/2 - kron( dz_n*Dz,kron( Iy,Ix ) )/2 ...
         + lambda/2*( kron( Iz,kron( Iy,dx_n ) )*( (kron( Iz,kron( Iy,Dx ) )*c_old).*kron( Iz,kron( Iy,Intn_x ) ) )...
                    + kron( Iz,kron( dy_n,Ix ) )*( (kron( Iz,kron( Dy,Ix ) )*c_old).*kron( Iz,kron( Intn_y,Ix ) ) ) ...
                    + kron( dz_n,kron( Iy,Ix ) )*( (kron( Dz,kron( Iy,Ix ) )*c_old).*kron( Intn_z,kron( Iy,Ix ) ) ) );
RHS_n = 1/tau*n_old + ( kron( Iz,kron( Iy,dx_n*Dx ) ) + kron( Iz,kron( dy_n*Dy,Ix ) ) + kron( dz_n*Dz,kron( Iy,Ix ) ) )*n_old/2 ...
        - lambda/2*( kron( Iz,kron( Iy,dx_n ) )*( (kron( Iz,kron( Iy,Dx ) )*c_old).*kron( Iz,kron( Iy,Intn_x ) ) )...
                    + kron( Iz,kron( dy_n,Ix ) )*( (kron( Iz,kron( Dy,Ix ) )*c_old).*kron( Iz,kron( Intn_y,Ix ) ) ) ...
                    + kron( dz_n,kron( Iy,Ix ) )*( (kron( Dz,kron( Iy,Ix ) )*c_old).*kron( Intn_z,kron( Iy,Ix ) ) ) )*n_old;
n_new = LHS_n\RHS_n;
%%%%%%%%%%%%%%%%%%%%generate variables%%%%%%%%%%%%%%%%%%%%%%%
c_old = c_new;
n_oold = n_old;
n_old = n_new;

end

c = reshape(c_new,[Nx,Ny,Nz]);
n = reshape(n_new,[Nx,Ny,Nz]);
