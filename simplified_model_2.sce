

clear; 

function dfds = fdd(f,s)
    
    n = length(f); 
    
//    ff = f(2:$-1); 
    fm = f(1:$-2); 
    fp = f(3:$);
     
//    ss = s(2:$-1); 
    sm = s(1:$-2); 
    sp = s(3:$); 

// disp([size(f) size(s)] )
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

function w=wakevelocity(G,dGds)

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

function [G,dGds,w] = circulation(A,u,b)

    G = 2*u*b * A'*theta_sintable; 
    dGds = 4*u * (A'*theta_ncostable)./(sin(theta) + eps);
    w = u * (A'*theta_nsintable)./(sin(theta) + eps);

endfunction; 


function f = induced_drag(x)

    lateral_distance = x(1); 
    long_distance    = x(2);
    ar               = x(3); 
    cl_avg           = x(4); 
    u                = x(5); 
    A                = x(6:(ne+6-1));

 // chord based on elliptic planform

    q    = 0.5 * rho * u * u;
    area = L / (q * cl_avg); 
    chord0 = ( area * 4) / (%pi * b);
    chord = chord0 * sqrt( 1 - (2*s/b).^2 ); 
    re    = rho * chord * u /nu; 

    [G,dGds,w] = circulation(A,u,b); 

// induced power 
    fi = -(rho * G' .* (ww )); 
    pwr_req = ti(fi',s) * u;     
    f(1,1) = pwr_req; // * 12 / power_density; 
    

// wing root bending moment 

    fi = rho * u * G .* abs(s) * 0.5;
    f(1,2) =  ti(fi,s);

endfunction; 



// numerics
eps=1e-13; 
dfds_max = 10000; 

// flow 
rho = 1.225; 
u = 20; 

// configuration 
ar =20; 
b=20; 
x0=-b; 
y0=0.0*b; 
z0=0; 

// mesh 
n=250; 
theta = linspace(0,%pi,n+1)+%pi/(n+1)*0.5;
theta = theta(1:$-1); 
s     = -(b/2)*cos(theta); 
ne    = 10; 
A = zeros(ne,1); 




//A = rand(ne,1)-0.5;
A(1) = 1; 

//A(4:10)=zeros(7,1);

nh = [1:ne]';


theta_sintable = sin(nh*theta); 
theta_costable = cos(nh*theta);
theta_nsintable = diag([1:ne])*theta_sintable; 
theta_ncostable = diag([1:ne])*theta_costable; 

A1 = A(1); A2=A(2); A3=A(3); 

 

l = 2*rho*u*u*b *(A1 * 1 * sin(theta) +  A2 * 1 * sin(2*theta) + A3 * 1 * sin(3*theta) ); 

ui = u * (A1 * 1 * sin(theta) +  A2 * 2 * sin(2*theta) + A3 * 3 * sin(3*theta) ) ./ sin(theta);
G = 2*u*b *(A1 * sin(theta) +  A2 * sin(2*theta) + A3 * sin(3*theta) ); 
dGds = 4*u * (A1 * 1 * cos(theta) +  A2 * 2 * cos(2*theta) + A3 * 3 * cos(3*theta) ) ./ sin(theta);

dgds = fdd(G,s); 


w_old = wakevelocity_old(n-1,b,s,G,u,x0,y0,z0);


[Gtest,dGdstest,wtest] = circulation(A,u,b); 
w=wakevelocity(G,dGds);






