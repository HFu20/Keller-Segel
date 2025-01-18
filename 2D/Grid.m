function [x,y,hf,kf,h,k,xf,yf]=Grid(hr,kr,Nx,Ny,x,y,non)
%random perturbation grids
beta = 0.5;

if non == 1
for i = 2:Nx          
    x(i) = x(i) + beta*hr*(-1+2*rand);
end
for i = 2:Ny          
    y(i) = y(i) + beta*kr*(-1+2*rand);
end
end

hf = zeros(1,Nx);       %
kf = zeros(1,Ny); 
for i = 1:Nx
    hf(i) = x(i+1) - x(i);
end
for i = 1:Ny
    kf(i) = y(i+1) - y(i);
end

xf = zeros(1,Nx);
yf = zeros(1,Ny);
for i = 1:Nx
    xf(i) = ( x(i) + x(i+1) )/2;     
end
for i = 1:Ny
    yf(i) = ( y(i) + y(i+1) )/2;     
end

h = zeros(1,Nx+1);       
k = zeros(1,Ny+1); 
for i = 1:Nx-1
    h(i+1) = xf(i+1) - xf(i);
end
h(1) = xf(1) - x(1);
h(Nx+1) = x(Nx+1) - xf(Nx);
for i = 1:Ny-1
    k(i+1) = yf(i+1) - yf(i);
end
k(1) = yf(1) - y(1);
k(Ny+1) = y(Ny+1) - yf(Ny);
end
