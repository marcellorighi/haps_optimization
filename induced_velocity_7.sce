
// G defined w/o cl and aoa

clear; 

function y=h1(x)
    y = x; 
endfunction

function y=h2(x)

//    y = -1 + x^2; 
    y = 0.5*(3*x.^2 -1); 

endfunction

function y=h3(x)

//    y = -3*x + x^3; 
    y = 0.5*(5*x.^3 - 3*x); 

endfunction

function y=h4(x)

//    y = 3 - 6*x^2 + x^4; 
    y = 0.125*(35*x.^4 - 30*x.^2 + 3); 

endfunction

function y=h5(x)

//    y = 15*x - 10*x^3 + x^5; 
    y = 0.125*(63*x.^5 - 70*x.^3 + 15*x ); 

endfunction

function [G,b,s,ds] = circulation(n,rho,u,L,ar,cl_avg,cla,A1,A2,A3,A4)

    q    = 0.5 * rho * u * u;
    area = L / (q * cl_avg); 
    b    = sqrt(ar*area); 
    G0 = 4 * L / (%pi * rho * u * b); 

    s = ([0:1:n] + 0.5) * b / (n+1) - b/2; 
    ds = b / (n+1); 
//    G = G0 * sqrt( 1 - (2*s/b).^2 ) + A1 * cos(%pi * 2*s/b) + A2 * cos(2 * %pi * 2*s/b) + A3 * cos(4 * %pi * 2*s/b) + A4 * (cos(1 * %pi * 1 *s/b) - 0.5);

    G = G0 * sqrt( 1 - (2*s/b).^2 ) + A1 * h1(2*s/b) + A2 * h2(2*s/b) + A3 * h3(2*s/b) + A4 * h4(2*s/b);

endfunction

function w=inducedvelocity(n,b,G)

    tt=0:%pi/n:%pi;
    s = -b/2 * cos(tt);
//    s = ([0:1:n] + 0.5) * b / (n+1) - b/2; 
    ds_aux  = ([s s($) s($)] - [s(1) s(1) s]);
//    ds = b / (n+1); 
    ds = ds_aux(2:$-1); 

    dGds_aux  = ([G G($) G($)] - [G(1) G(1) G]) ./ ds_aux; 

    dGds = dGds_aux(2:$-1); 
    dGds(1) = dGds(2); 
    dGds($) = dGds($-1); 

    eps=0.0001;

    for j = 1:n+1; 
    
        y = s(j); 
        w_s = dGds .* (-y + s) ./((-y + s).^2 + eps) * 1/( 4*%pi );
//        w(j)= sum(w_s) * ds; 
        w(j)= sum(w_s .* ds) - 0.5 * (w_s(1)*ds(1) + w_s($)*ds($)); 

    end
    
endfunction;
    
function w=wakevelocity(n,b,bL,G,u,x0,y0,z0)

// bL = b lead A/C

// x, y, z = non dimensional distance between two A/C 

    x1 = x0 * b; 
    y1 = y0 * b; 
    z1 = z0 * b; 

    nL = length(G) - 1; 
// sL = ([0:1:nL] + 0.5) * bL / (nL+1) - bL/2; 
// dsL = bL / (nL+1); 

    tt=0:%pi/nL:%pi;
    sL = -bL/2 * cos(tt);
    dsL_aux  = ([sL sL($) sL($)] - [sL(1) sL(1) sL]);
    dsL = dsL_aux(2:$-1); 


    tt=0:%pi/n:%pi;
    s = -b/2 * cos(tt);
    ds_aux  = ([s s($) s($)] - [s(1) s(1) s]);
    ds = ds_aux(2:$-1); 

//s = ([0:1:n] + 0.5) * b / (n+1) - b/2; 
//ds = b / (n+1); 

    dGds_aux  = ([G G($) G($)] - [G(1) G(1) G]) ./ ds_aux; 
    dGds = dGds_aux(2:$-1); 

// dGds_aux  = ([G G($) G($)] - [G(1) G(1) G]) ./ ds; 
// dGds = dGds_aux(2:$-1); 
    dGds(1) = dGds(2); 
    dGds($) = dGds($-1); 

tau = -x1 / u;
rc  = 0.52 * sqrt(tau);

for j = 1:n+1; 
    
    y = s(j) + y1; 

    w_s = dGds .* (- y + sL) ./ ((y - sL).^2 + z1^2 + rc^2) .* (1 - x1 *ones(sL)./sqrt(x1^2 + (y - sL).^2 + z1^2) ) * 1/( 4*%pi );
    w(j)= sum(w_s .* dsL) - 0.5 * (w_s(1)*dsL(1) + w_s($)*dsL($)); 

end
    
endfunction;

function f = induceddrag(x)

    lateral_distance = x(1); 
    A1_LA            = x(2);
    A2_LA            = x(3);
    A3_LA            = x(4);
    A4_LA            = x(5);
    B1_LA            = 0.0; 
    B2_LA            = 0.0; 
    B3_LA            = 0.0; 
    A1               = x(6);
    A2               = x(7);
    A3               = x(8);
    A4               = x(9);
/*    B1               = x(10); 
    B2               = x(11); 
    B3               = x(12); 
*/
    n = 201; 
    rho = 0.08; 
    u = 20.0; 
//    L = 750.0; 
    L_LA = 850; 
    L    = 850; 
    ar = 10; 
    cl_avg = 0.75; 
    cla = 2*%pi; 
 
    [G,b,s,ds] = circulation(n,rho,u,L,ar,cl_avg,cla,A1,A2,A3,A4);
    ww    = inducedvelocity(n,b,G);
    [G_LA,b,s,ds] = circulation(n,rho,u,L_LA,ar,cl_avg,cla,A1_LA,A2_LA,A3_LA,A4_LA);
    ww2   = wakevelocity(n,b,b,G_LA,u,-11,lateral_distance,0);
    ww_LA = inducedvelocity(n,b,G_LA);

    f = -2*sum(G * rho * (ww + ww2))*ds - sum(G_LA * rho * ww_LA)*ds;
//    f = -1*sum(G * rho * ww ) * ds;
    
    di_solo =  -2*sum(G * rho * ww)*ds;

    disp([x(1) x(2) x(5) f di_solo]);

    scf(f1); clf; plot(s,ww,s,ww2,'r-'); 
    scf(f2); clf; plot(s,ww_LA);
    scf(f3); clf; plot(s,G); ylabel('$ G $'); 

endfunction; 


function [f, g, ind]=cost(x, ind)
    f = induceddrag (x);
    g = numderivative (induceddrag, x);
endfunction


/*
n = 201; 
rho = 0.08; 
u = 12.0; 
L = 750.0; 
ar = 30; 
cl_avg = 0.75; 
cla = 2*%pi; 
*/

    
/*
[G,b,s,ds] = circulation(n,rho,u,L,ar,cl_avg,cla,0.001,0.001,0.001);
ww    = inducedvelocity(n,b,G);
ww2   = wakevelocity(n,b,b,G,u,-11,0.95,0);

di = sum(G * rho * (ww - ww2));
*/
  
x0=[1.5; 0.0001; 0.0001; 0.0001; 0.0001; 0.0000; 0.0000; 0.0000; 0.0001]; 
 
x_min = [0.50; -1.5; -1.5; -1.5; -1.5; -1.50; -1.50; -1.50; -1.5];  
x_max=  [2.00; +1.5; +1.5; +1.5; +1.5; +1.5; +1.5; +1.5; +1.5];  

f1 = scf(1); clf; 
f2 = scf(2); clf; 
f3 = scf(3); clf; 

[fopt, xopt] = optim (cost, "b", x_min, x_max, x0)    
    
    
