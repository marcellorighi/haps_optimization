

clear; 

function dfds = fdd(f,s)
    
    n = length(f); 
    
//    ff = f(2:$-1); 
    fm = f(1:$-2); 
    fp = f(3:$);
     
//    ss = s(2:$-1); 
    sm = s(1:$-2); 
    sp = s(3:$); 

    dfds = (fp - fm) ./ (sp - sm);
    
/*    for i=2:n-1 
        dfds(1,i) = (f(i+1) - f(i-1)) / (s(i+1) - s(i-1)); 
    end
*/
    
    dfds(1,1) = (f(2) - f(1)) / (s(2) - s(1));
    dfds(1,n) = (f(n) - f(n-1)) / (s(n) - s(n-1));
    
    dfds = max(min(dfds,dfds_max),-dfds_max);

endfunction

function intfds=ti(f,s)
    
    n = length(f); 
    intfds = 0.; 
  
    ff = f(1:$-1); 
    fp = f(2:$); 
    ss = s(1:$-1); 
    sp = s(2:$); 

    intfds = 0.5 * (ff + fp) * (sp - ss)'; 

/*    for i=1:n-1 
        intfds = intfds + 0.5 * (f(i) + f(i+1)) * (s(i+1) - s(i)); 
    end
*/
    
endfunction

function w=wakevelocity_old(n,b,s,G,u,x0,y0,z0)

// bL = b lead A/C

// x, y, z = non dimensional distance between two A/C 

bL = b; 

x1 = x0 * b; 
y1 = y0 * b; 
z1 = z0 * b; 

nL = length(G) - 1; 

dGds = fdd(G,s); 


tau = -x1 / u;
rc  = 0.52 * sqrt(tau);

for j = 1:n+1; 
    
    y = s(j) + y1; 

    w_s = dGds .* (- y + s) ./ ((y - s).^2 + z1^2 + rc^2) .* (1 - x1 *ones(s)./sqrt(x1^2 + (y - s).^2 + z1^2) ) * 1/( 4*%pi );

//    w(j)= (sum(w_s) - 0.5*(w_s(1) + w_s($))) * dsL; 

    w(j) = ti(w_s,s);

end
    
endfunction;

function w=wakevelocity(n,b,s,G,dGds,u,x0,y0,z0)

x1 = x0 * b; 
y1 = y0 * b; 
z1 = z0 * b; 

tau = -x1 / u;
rc  = 0.52 * sqrt(tau);

for j = 1:n; 
    
    y = s(j) + y1; 
    w_s = dGds .* (- y + s) ./ ((y - s).^2 + z1^2 + rc^2) .* (1 - x1 *ones(s)./sqrt(x1^2 + (y - s).^2 + z1^2) ) * 1/( 4*%pi );
    w(j) = ti(w_s,s);

end
    
endfunction;

function [G,w] = circulation(A)

    G = 2*u*b * A'*theta_sintable; 
    w = u * (A'*theta_nsintable)./(sin(theta) + eps);

endfunction; 

n=250; 
rho = 1.225; 
b=20; 
u = 20; 
ar =20; 
dfds_max = 10000; 

x0=-b; 
y0=0.0*b; 
z0=0; 


ne=10; 
A = zeros(ne,1); 
A(1) = 1; 

nh = [1:ne]';

theta = linspace(0,%pi,n+1)+%pi/(n+1)*0.5;
theta = theta(1:$-1); 

theta_sintable = sin(nh*theta); 
theta_costable = cos(nh*theta);
theta_nsintable = diag([1:ne])*theta_sintable; 
theta_ncostable = diag([1:ne])*theta_costable; 

A1 = A(1); A2=A(2); A3=A(3); 

 
s     = -(b/2)*cos(theta); 

l = 2*rho*u*u*b *(A1 * 1 * sin(theta) +  A2 * 1 * sin(2*theta) + A3 * 1 * sin(3*theta) ); 

ui = u * (A1 * 1 * sin(theta) +  A2 * 2 * sin(2*theta) + A3 * 3 * sin(3*theta) ) ./ sin(theta);
G = 2*u*b *(A1 * sin(theta) +  A2 * sin(2*theta) + A3 * sin(3*theta) ); 
dGds = 4*u * (A1 * 1 * cos(theta) +  A2 * 2 * cos(2*theta) + A3 * 3 * cos(3*theta) ) ./ sin(theta);

dgds = fdd(G,s); 

w=wakevelocity(n,b,s,G,dGds,u,x0,y0,z0);

w_old = wakevelocity_old(n-1,b,s,G,u,x0,y0,z0);


[Gtest,wtest] = circulation(A); 






