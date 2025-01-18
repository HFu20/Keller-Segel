function  A = Gradient(a,M)

A = zeros(M-1,M);

A(1:M:end-1) = -1./a(2:M);

for i = 1:M-1
    A(i,i+1) = 1./a(i+1);
end

end
