function  [Dd,dD ]= Laplace(a,b,M)

Dd = diag( -1./a(2:M).*( 1./b(2:M) + 1./b(1:M-1) ) )...
    + diag( 1./(a(2:M-1).*b(2:M-1)) , 1 )...
    + diag( 1./(a(3:M).*b(2:M-1)) , -1 );

dD = diag( -1./b(1:M).*( 1./a(2:M+1) + 1./a(1:M) ) )...
    + diag( 1./(b(1:M-1).*a(2:M)) , 1 )...
    + diag( 1./(b(2:M).*a(2:M)) , -1 );

end
