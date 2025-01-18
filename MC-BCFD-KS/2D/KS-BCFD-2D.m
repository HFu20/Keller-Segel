%%solve problem Keller-Segel system in 2D
%Dec-MC-BCFD method
%=============================================================
clc,clear

%%%%%%%%%%%%%%%%some input
M = 20;
hr = 1/M;
% hr = 0.025;
kr = hr;
xL = 0;
xR = 1;
yL = 0;
yR = 1;
T = 6e-5;
non=1;   %1-non-uniform grid, 0-uniform grid

lambda = 1;
alpha = 1;
%%%%%%%%%%%%%%%temporal grid information
tau = 1e-06;
t = 0 : tau : T;    

%%%%%%%%%%%%%%%%spatial grid information
x = xL : hr : xR;
y = yL : kr : yR;
N = length(t);
Nx = length(x) - 1;
Ny = length(y) - 1;
%%
% [x,y,hf,kf,h,k,xf,yf] = Grid(hr,kr,Nx,Ny,x,y,non); %random perturbation grids
% [x,y,hf,kf,h,k,xf,yf] = Grid_cor(hr,kr,Nx,Ny,x,y,non); %cornor refinement grids
[x,y,hf,kf,h,k,xf,yf] = Grid_mid(Nx,Ny,x,y,non);  %middle refinement grids

c_end = zeros(Nx,Ny,N);
n_end = zeros(Nx,Ny,N);
%%%%%%%%%%%%%Initial condtion
[c_old,n_old] = Initial_solution(xf,yf);
c_end(:,:,1) = c_old;
n_end(:,:,1) = n_old;

c_old = c_old(:);
n_old = n_old(:);
%%%%%%%%%%%%%%%laplace operator aproximation
[Ddx,dDx] = Laplace(h,hf,Nx);
[Ddy,dDy] = Laplace(k,kf,Ny);

%%%%%%%%%%%%%%gradient operator aproximation
Dx = gradient(h,Nx);
Dy = gradient(k,Ny);
dx = gradient(hf,Nx-1);
dy = gradient(kf,Ny-1);

%%%%%%%%%%%%Interpolation operator
Interp_x = Interpolation(h,hf,Nx-1);
Interp_y = Interpolation(k,kf,Ny-1);

%%%%%%%%%%%%%%%%%%%%%%%%%
dx_n  = zeros(Nx,Nx-1);
dy_n  = zeros(Ny,Ny-1);
dx_n(2:Nx-1,:)  = dx;
dy_n(2:Ny-1,:)  = dy;
dx_n(1,1) = 1/hf(1);
dx_n(Nx,Nx-1) = -1/hf(Nx);
dy_n(1,1) = 1/kf(1);
dy_n(Ny,Ny-1) = -1/kf(Ny);
%%%%%%%%%%%%%%%%%%%
Intn_x = Interpolation(h,hf,Nx);
Intn_y = Interpolation(k,kf,Ny);
%%
Iy = eye(Ny);
Ix = eye(Nx);
%%%%%%%%%Step 1: compute n %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LHS_n = 1/tau*kron( Iy,Ix ) - kron( Iy,dx_n*Dx )/2 - kron( dy_n*Dy,Ix )/2 ...
         + lambda*( kron( Iy,dx_n )*( (kron( Iy,Dx )*c_old).*kron( Iy,Intn_x ) )...
                    + kron( dy_n,Ix )*( (kron( Dy,Ix )*c_old).*kron( Intn_y,Ix ) )   );
RHS_n = 1/tau*n_old +  ( kron( Iy,dx_n*Dx ) + kron( dy_n*Dy,Ix ) )*n_old/2;
n_new = LHS_n\RHS_n;
%%%%%%%%%Step 2: compute c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LHS_c = (1/tau+1/2)*kron( Iy,Ix ) ... 
        - ( kron( Iy,dx_n*Dx ) + kron( dy_n*Dy,Ix ) )/2;
RHS_c = (1/tau-1/2)*c_old  + ( kron( Iy,dx_n*Dx ) + kron( dy_n*Dy,Ix ) )*c_old/2 + (n_new+n_old)/2;
c_new = LHS_c\RHS_c;
%%%%%%%%%%%%%%%%%%%%generate variables%%%%%%%%%%%%%%%%%%%%%%%
c_old = c_new;
n_oold = n_old;
n_old = n_new;

for m = 3:N
    Time = t(m);
%%%%%%%%%Step 1: compute c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LHS_c = (1/tau+alpha/2)*kron( Iy,Ix ) ... 
        - ( kron( Iy,dx_n*Dx ) + kron( dy_n*Dy,Ix ) )/2;
RHS_c = (1/tau-alpha/2)*c_old  + ( kron( Iy,dx_n*Dx ) + kron( dy_n*Dy,Ix ) )*c_old/2 + (3/2*n_old - 1/2*n_oold);% + fc1;
c_new = LHS_c\RHS_c;
%%%%%%%%%Step 2: compute n %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LHS_n = 1/tau*kron( Iy,Ix ) - kron( Iy,dx_n*Dx )/2 - kron( dy_n*Dy,Ix )/2 ...
         + lambda/2*( kron( Iy,dx_n )*( (kron( Iy,Dx )*c_new).*kron( Iy,Intn_x ) )...
                    + kron( dy_n,Ix )*( (kron( Dy,Ix )*c_new).*kron( Intn_y,Ix ) )   );
RHS_n = 1/tau*n_old +  ( kron( Iy,dx_n*Dx ) + kron( dy_n*Dy,Ix ) )*n_old/2 ...
        - lambda/2*( kron( Iy,dx_n )*( (kron( Iy,Dx )*c_new).*kron( Iy,Intn_x ) )...
                    + kron( dy_n,Ix )*( (kron( Dy,Ix )*c_new).*kron( Intn_y,Ix ) )   ) *n_old;% + fn2;
n_new = LHS_n\RHS_n;
%%%%%%%%%%%%%%%%%%%%generate variables%%%%%%%%%%%%%%%%%%%%%%%
c_old = c_new;
n_oold = n_old;
n_old = n_new;

end

c = reshape(c_new,[Nx,Ny]);
n = reshape(n_new,[Nx,Ny]);
