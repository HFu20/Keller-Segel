function  A = Interpolation(a,b,M)

A = zeros(M-1,M);

A(1:M:end-1) = b(2:M)./(2*a(2:M));

for i = 1:M-1
    A(i,i+1) = b(i)./(2*a(i+1));
end
end
