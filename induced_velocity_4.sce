
clear; 

function [G,b] = circulation(n,rho,u,L,ar,cl_avg,cla,A1,A2,A3)

    q    = 0.5 * rho * u * u;
    area = L / (q * cl_avg); 
    b    = sqrt(ar*area); 
    chord0 = area / (%pi/4 * b); 

    s = ([0:1:n] + 0.5) * b / (n+1) - b/2; 
    ds = b / (n+1); 
    chord = chord0 * sqrt( 1 - (2*s/b).^2 );

    twist = A1 * cos(%pi * 2*s/b) + A2 * cos(%pi/2 * 2*s/b) + A3 * cos(%pi/4 * 2*s/b); 
    alpha0 = L/(q*area*cla) - sum(twist.*chord)*ds/area; 

    alpha = alpha0 + twist; 
    cl    = cla * alpha; 
    G     = cl .* chord * u/2; 

endfunction

function w=inducedvelocity(n,b,G)

    s = ([0:1:n] + 0.5) * b / (n+1) - b/2; 
    ds = b / (n+1); 

    dGds_aux  = ([G G($) G($)] - [G(1) G(1) G]) / ds; 
    dGds = dGds_aux(2:$-1); 
    dGds(1) = dGds(2); 
    dGds($) = dGds($-1); 

    eps=0.01;

    for j = 1:n+1; 
    
        y = s(j); 
        w_s = dGds .* (-y + s) ./((-y + s).^2 + eps) * 1/( 4*%pi );
        w(j)= sum(w_s) * ds; 

    end
    
endfunction;
    
function w=wakevelocity(n,b,G,u,x0,y0,z0)

// x, y, z = non dimensional distance between two A/C 

x1 = x0 * b; 
y1 = y0 * b; 
z1 = z0 * b; 

s = ([0:1:n] + 0.5) * b / (n+1) - b/2; 
ds = b / (n+1); 

dGds_aux  = ([G G($) G($)] - [G(1) G(1) G]) / ds; 
dGds = dGds_aux(2:$-1); 
dGds(1) = dGds(2); 
dGds($) = dGds($-1); 

tau = -x1 / u;
rc  = 0.52 * sqrt(tau);

for j = 1:n+1; 
    
    y = s(j) + y1; 

    w_s = dGds .* (- y + s) ./ ((y - s).^2 + z1^2 + rc^2) .* (1 - x1 *ones(s)./sqrt(x1^2 + (y - s).^2 + z1^2) ) * 1/( 4*%pi );
    w(j)= sum(w_s) * ds; 

end
    
endfunction;

n = 201; 
rho = 0.08; 
u = 12.0; 
L = 750.0; 
ar = 30; 
cl_avg = 0.75; 
cla = 2*%pi; 
    
[G,b] = circulation(n,rho,u,L,ar,cl_avg,cla,0.001,0.001,0.001);
ww    = inducedvelocity(n,b,G);
ww2   = wakevelocity(n,b,G,u,-11,1.00,0);


    
    
    
