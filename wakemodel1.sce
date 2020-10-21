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

// 
n = 201; 
s = ([0:1:n] + 0.5) * b / (n+1) - b/2; // +  linspace(-b/2,b/2,n+1); 
ds = b / (n+1); 
G = G0 * sqrt( 1 - (2*s/b).^2 ); 
dGds = (1/2) * G0 * (-2*(2*s/b) * 2/b) ./ sqrt( 1 - (2*s/b).^2 ) ;

// 
x = -1*b;

z = 0.0;

tau = -x / u;
// rc     = 20.0; //2.24 * sqrt(mu * tau_mu);
rc  = 0.52 * sqrt(tau);

ny = 250; 
yy = linspace(3*b/2,0*b/2,ny); 

//f2=scf(2); 
clf; 


for j = 1:ny; 
    
y = yy(j); 
w_s = dGds .* (y - s) ./ ((y - s).^2 + z^2 + rc^2) .* (1 - x * ones(s)./sqrt(x^2 + (y - s).^2 + z^2) ) * 1/( 4*%pi );
w(j)= sum(w_s) * ds; 

//plot(s,w_s); 

end




f1 = scf(1); 
clf; 
plot(yy/b,w/u,'o-'); 
xgrid; 








