function [x,y,z,hf,kf,lf,h,k,l,xf,yf,zf] = Grid_mid_3D(Nx,Ny,Nz,x,y,z,non)

if non == 1
% for i = 2:Nx          
%     x(i) = x(i) + 0.1*hr*(-1+2*rand);
% end
% for i = 2:Ny          
%     y(i) = y(i) + 0.1*kr*(-1+2*rand);
% end
% end

%%%% mid mesh encryption
for i = 1:Nx/2
x((Nx)/2+i) = 1/2 + (i-1)^(2)/(2*(Nx/2+1)^(2));
end
for i = 1:Nx/2-1
x((Nx)/2-i+1) = 1/2 - i^(2)/(2*(Nx/2+1)^(2));
end

for i = 1:Ny/2
y((Ny)/2+i) = 1/2 + (i-1)^(2)/(2*(Ny/2+1)^(2));
end
for i = 1:Nx/2-1
y((Ny)/2-i+1) = 1/2 - i^(2)/(2*(Ny/2+1)^(2));
end

for i = 1:Nz/2
z((Nz)/2+i) = 1/2 + (i-1)^(2)/(2*(Nz/2+1)^(2));
end
for i = 1:Nx/2-1
z((Nz)/2-i+1) = 1/2 - i^(2)/(2*(Nz/2+1)^(2));
end
end

hf = zeros(1,Nx);       %
kf = zeros(1,Ny); 
lf = zeros(1,Nz); 
for i = 1:Nx
    hf(i) = x(i+1) - x(i);
end
for i = 1:Ny
    kf(i) = y(i+1) - y(i);
end
for i = 1:Nz
    lf(i) = z(i+1) - z(i);
end

xf = zeros(1,Nx);
yf = zeros(1,Ny);
zf = zeros(1,Nz);
for i = 1:Nx
    xf(i) = ( x(i) + x(i+1) )/2;     
end
for i = 1:Ny
    yf(i) = ( y(i) + y(i+1) )/2;     
end
for i = 1:Ny
    zf(i) = ( z(i) + z(i+1) )/2;     
end

h = zeros(1,Nx+1);       %半网格空间步长
k = zeros(1,Ny+1); 
l = zeros(1,Nz+1);
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
for i = 1:Nz-1
    l(i+1) = zf(i+1) - zf(i);
end
l(1) = zf(1) - z(1);
l(Ny+1) = z(Ny+1) - zf(Ny);
end
