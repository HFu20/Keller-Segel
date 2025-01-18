function  [c,n] = Initial_solution(xf,yf)
%choose the initial value below

%no blow up less than 8*pi, initial value for Example 4.2
n = 50*exp( -5*((xf-0.5).^2'+ (yf-0.5).^2) );
c = 25*exp( -2.5*((xf-0.5).^2'+ (yf-0.5).^2) );

%blow up great than 8*pi, initial value for Example 4.3
n = 130*exp( -15*((xf).^2'+ (yf).^2) );
c = 13*exp( -2*((xf).^2'+ (yf).^2) );

%blow up at mid, initial value for Example 4.4
n = 1000*exp( -100*((xf-0.5).^2' + (yf-0.5).^2) );
c = 500*exp( -50*((xf-0.5).^2' + (yf-0.5).^2) );

%blow up at corner, initial value for Example 4.5
n = 1000*exp( -100*((xf-0.15).^2'+ (yf-0.15).^2) );
c = zeros(length(xf),length(yf));

end
