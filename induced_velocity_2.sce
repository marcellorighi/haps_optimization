
clear; 
// 
rho = 0.08; 
mu  = 1.5e-5; 
nu  = mu/rho; 
u   = 18; 

q = 0.5 * rho * u * u; 
L   = 750;  
ar  = 28; 
cl_avg = 0.5; 
area = L / (q * cl_avg); 
b    = sqrt(ar*area); 
chord0 = area / (%pi/4 * b); 

// b   = 28.; 
// a   = %pi/8*b; 

// leading A/C
x = -10*b; 
z = 0; 
tau = -x / u;
rc  = 0.52 * sqrt(tau);

G0  = 4 * L / (%pi * rho * u * b); 






cla = 2*%pi; 
// alpha0 = cl_avg / cla; 

// twist, Fourier series 

A1 = 0.05; 
A2 = 0.02; 
A3 = 0.01; 


// 
n = 201; 
s = ([0:1:n] + 0.5) * b / (n+1) - b/2; // +  linspace(-b/2,b/2,n+1); 
ds = b / (n+1); 
chord = chord0 * sqrt( 1 - (2*s/b).^2 );

twist = A1 * cos(%pi * 2*s/b) + A2 * cos(%pi/2 * 2*s/b) + A3 * cos(%pi/4 * 2*s/b); 
alpha0 = L/(q*area*cla) - sum(twist.*chord)*ds/area; 

alpha = alpha0 + twist; 
cl    = cla * alpha; 
G     = cl .* chord * u/2; 
dGds_aux  = ([G G($) G($)] - [G(1) G(1) G]) / ds; 
dGds = dGds_aux(2:$-1); 
dGds(1) = dGds(2); 
dGds($) = dGds($-1); 


// Leading A/C
G_LA = G0 * sqrt( 1 - (2*s/b).^2 ); 
dGds_LA = (1/2) * G0 * (-2*(2*s/b) * 2/b) ./ sqrt( 1 - (2*s/b).^2 ) ;



// G = G0 * sqrt( 1 - (2*s/b).^2 ); 
// dGds = (1/2) * G0 * (-2*(2*s/b) * 2/b) ./ sqrt( 1 - (2*s/b).^2 ) ;


// induced velocity

eps=0.01;

f2=scf(2); clf; 

for j = 1:n+1; 
    
    y = s(j); 
    y_LA = s(j) + b*0.88; 

// induced velocity form LA

    w_s_LA = dGds_LA .* (y_LA - s) ./ ((y_LA - s).^2 + z^2 + rc^2) .* (1 - x * ones(s)./sqrt(x^2 + (y_LA - s).^2 + z^2) ) * 1/( 4*%pi );

    w_LA(j)= sum(w_s_LA) * ds; 
    
// own induced velocity     

    w_s = dGds .* (-y + s) ./((-y + s).^2 + eps) * 1/( 4*%pi );
    w(j)= sum(w_s) * ds; 

    ddi_FORM(j) = G(j) * rho * (w(j) - w_LA(j));
    ddi(j) = G(j) * rho * (w(j) );

    plot(s,w_s);
    
end

f1=scf(1); clf; 
plot(s',w,s',w_LA,'r-'); 
legend('$\Large w_{own}$','$\Large w^{LA}$');
xgrid; 

f3=scf(3); clf; 
plot(s',ddi); 
plot(s',ddi_FORM,'r-'); 
legend('$\Large D_i$','$\Large D_i^F$');
xgrid; 

f4=scf(4); clf; 
plot(s,G,s,G_LA,'r-');
legend('$\Large \Gamma$','$\Large \Gamma^{LA}$');
xgrid; 



