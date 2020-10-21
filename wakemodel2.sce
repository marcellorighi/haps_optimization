
clear; 
// 
rho = 1.00; 
mu  = 1.5e-5; 
nu  = mu/rho; 
u   = 100; 
L   = 502649;  
b   = 40.; 
a   = %pi/8*b; 
G0  = 4 * L / (%pi * rho * u * b); 

x = 11*b;
z = 0.0;
ny = 50; 
yy = linspace(0.0*b,2*b,ny); 

for j = 1:ny; 
    
y   = yy(j);    
w(j)= - G0/(4*%pi*x) * ( (y+a)./ sqrt(x^2 + (y+a).^2) - (y-a)./ sqrt(x^2 + (y-a).^2) ) - G0/(2*%pi) * a/(y^2 -a^2) * (1 - x/sqrt(x^2 + (y+a)^2)); 

end

f1 = scf(1); 
clf; 
plot(yy,w,'o-'); 
xgrid; 



