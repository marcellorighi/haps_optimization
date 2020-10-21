
clear; 
// 
rho = 0.08; 
mu  = 1.5e-5; 
nu  = mu/rho; 
u   = 18; 
L   = 750;  
b   = 28.; 
a   = %pi/8*b; 
G0  = 4 * L / (%pi * rho * u * b); 

// 
n = 201; 
s = ([0:1:n] + 0.5) * b / (n+1) - b/2; // +  linspace(-b/2,b/2,n+1); 
ds = b / (n+1); 
G = G0 * sqrt( 1 - (2*s/b).^2 ); 
dGds = (1/2) * G0 * (-2*(2*s/b) * 2/b) ./ sqrt( 1 - (2*s/b).^2 ) ;

//
ny = n+1; 
// yy = linspace(-b/2,b/2,ny); 

eps=0.01;

f2=scf(2); clf; 

for j = 1:ny; 
    
    y = s(j); 
    w_s = dGds .* (-y + s) ./((-y + s).^2 + eps) * 1/( 4*%pi );
    w(j)= sum(w_s) * ds; 

    ddi(j) = G(j) * rho * w(j);

    plot(s,w_s);
    
end

f1=scf(1); clf; 
plot(s',w); 
xgrid; 


