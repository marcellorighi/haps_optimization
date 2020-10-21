/*
// 19.08.20 stretched grid
// G defined w/o cl and aoa
// GA
// 24.08 Re-aoa dependent drag added
*/

clear; 

function Hijkl=Psi(x1,x2,x3,x4)
    
// 0th order 
    Hijkl(1)=1.; 

// 1st order
    Hijkl(2)=h1(x1);
    Hijkl(3)=h1(x2);
    Hijkl(4)=h1(x3);
    Hijkl(5)=h1(x4);
    
// 2nd order 
    ii=5; 
    Hijkl(ii+1)=h2(x1);
    Hijkl(ii+2)=h2(x2);
    Hijkl(ii+3)=h2(x3);
    Hijkl(ii+4)=h2(x4);
    
    Hijkl(ii+5) =h1(x1)*h1(x2);
    Hijkl(ii+6) =h1(x1)*h1(x3);
    Hijkl(ii+7) =h1(x1)*h1(x4);
    Hijkl(ii+8) =h1(x2)*h1(x3);
    Hijkl(ii+9) =h1(x2)*h1(x4);
    Hijkl(ii+10)=h1(x3)*h1(x4);
    
// 3rd order 
    ii=15; 
    Hijkl(ii+1)=h3(x1);
    Hijkl(ii+2)=h3(x2);
    Hijkl(ii+3)=h3(x3);
    Hijkl(ii+4)=h3(x4);

    Hijkl(ii+5) =h1(x1)*h1(x2)*h1(x3);
    Hijkl(ii+6) =h1(x1)*h1(x3)*h1(x4);
    Hijkl(ii+7) =h1(x1)*h1(x2)*h1(x4);
    Hijkl(ii+8) =h1(x2)*h1(x3)*h1(x4);

    Hijkl(ii+9) =h2(x1)*h1(x2);
    Hijkl(ii+10)=h2(x1)*h1(x3);
    Hijkl(ii+11)=h2(x1)*h1(x4);
    Hijkl(ii+12)=h2(x2)*h1(x1);
    Hijkl(ii+13)=h2(x2)*h1(x3);
    Hijkl(ii+14)=h2(x2)*h1(x4);
    Hijkl(ii+15)=h2(x3)*h1(x1);
    Hijkl(ii+16)=h2(x3)*h1(x2);
    Hijkl(ii+17)=h2(x3)*h1(x4);
    Hijkl(ii+18)=h2(x4)*h1(x1);
    Hijkl(ii+19)=h2(x4)*h1(x2);
    Hijkl(ii+20)=h2(x4)*h1(x3);

// 4th order 
    ii=35; 
    Hijkl(ii+1)=h4(x1);
    Hijkl(ii+2)=h4(x2);
    Hijkl(ii+3)=h4(x3);
    Hijkl(ii+4)=h4(x4);

    Hijkl(ii+5) =h3(x1)*h1(x2);
    Hijkl(ii+6) =h3(x1)*h1(x3);
    Hijkl(ii+7) =h3(x1)*h1(x4);
    Hijkl(ii+8) =h3(x2)*h1(x1);
    Hijkl(ii+9) =h3(x2)*h1(x3);
    Hijkl(ii+10)=h3(x2)*h1(x4);
    Hijkl(ii+11) =h3(x3)*h1(x1);
    Hijkl(ii+12) =h3(x3)*h1(x2);
    Hijkl(ii+13)=h3(x3)*h1(x4);
    Hijkl(ii+14) =h3(x4)*h1(x1);
    Hijkl(ii+15) =h3(x4)*h1(x2);
    Hijkl(ii+16)=h3(x4)*h1(x3);
    
    Hijkl(ii+17)=h2(x1)*h2(x2);
    Hijkl(ii+18)=h2(x1)*h2(x3);
    Hijkl(ii+19)=h2(x1)*h2(x4);
    Hijkl(ii+20)=h2(x2)*h2(x3);
    Hijkl(ii+21)=h2(x2)*h2(x4);
    Hijkl(ii+22)=h2(x3)*h2(x4);
    
    Hijkl(ii+23)=h2(x1)*h1(x2)*h1(x3);
    Hijkl(ii+24)=h2(x1)*h1(x2)*h1(x4);
    Hijkl(ii+25)=h2(x1)*h1(x3)*h1(x4);
    Hijkl(ii+26)=h2(x2)*h1(x1)*h1(x3);
    Hijkl(ii+27)=h2(x2)*h1(x1)*h1(x4);
    Hijkl(ii+28)=h2(x2)*h1(x3)*h1(x4);
    Hijkl(ii+29)=h2(x3)*h1(x1)*h1(x2);
    Hijkl(ii+30)=h2(x3)*h1(x1)*h1(x4);
    Hijkl(ii+31)=h2(x3)*h1(x2)*h1(x4);
    Hijkl(ii+32)=h2(x4)*h1(x1)*h1(x2);
    Hijkl(ii+33)=h2(x4)*h1(x1)*h1(x3);
    Hijkl(ii+34)=h2(x4)*h1(x2)*h1(x3);

    Hijkl(ii+35)=h1(x1)*h1(x2)*h1(x3)*h1(x4);

// 5th order 
    ii=70; 
    Hijkl(ii+1)=h5(x1);
    Hijkl(ii+2)=h5(x2);
    Hijkl(ii+3)=h5(x3);
    Hijkl(ii+4)=h5(x4);

    Hijkl(ii+5) =h4(x1)*h1(x2);
    Hijkl(ii+6) =h4(x1)*h1(x3);
    Hijkl(ii+7) =h4(x1)*h1(x4);
    Hijkl(ii+8) =h4(x2)*h1(x1);
    Hijkl(ii+9) =h4(x2)*h1(x3);
    Hijkl(ii+10)=h4(x2)*h1(x4);
    Hijkl(ii+11)=h4(x3)*h1(x1);
    Hijkl(ii+12)=h4(x3)*h1(x2);
    Hijkl(ii+13)=h4(x3)*h1(x4);
    Hijkl(ii+14)=h4(x4)*h1(x1);
    Hijkl(ii+15)=h4(x4)*h1(x2);
    Hijkl(ii+16)=h4(x4)*h1(x3);
        
    Hijkl(ii+17)=h3(x1)*h2(x2);
    Hijkl(ii+18)=h3(x1)*h2(x3);
    Hijkl(ii+19)=h3(x1)*h2(x4);
    Hijkl(ii+20)=h3(x2)*h2(x1);
    Hijkl(ii+21)=h3(x2)*h2(x3);
    Hijkl(ii+22)=h3(x2)*h2(x4);
    Hijkl(ii+23)=h3(x3)*h2(x1);
    Hijkl(ii+24)=h3(x3)*h2(x2);
    Hijkl(ii+25)=h3(x3)*h2(x4);
    Hijkl(ii+26)=h3(x4)*h2(x1);
    Hijkl(ii+27)=h3(x4)*h2(x2);
    Hijkl(ii+28)=h3(x4)*h2(x3);

    Hijkl(ii+29)=h3(x1)*h1(x2)*h1(x3);
    Hijkl(ii+30)=h3(x1)*h1(x2)*h1(x4);
    Hijkl(ii+31)=h3(x1)*h1(x3)*h1(x4);
    Hijkl(ii+32)=h3(x2)*h1(x1)*h1(x3);
    Hijkl(ii+33)=h3(x2)*h1(x1)*h1(x4);
    Hijkl(ii+34)=h3(x2)*h1(x3)*h1(x4);
    Hijkl(ii+35)=h3(x3)*h1(x1)*h1(x2);
    Hijkl(ii+36)=h3(x3)*h1(x1)*h1(x4);
    Hijkl(ii+37)=h3(x3)*h1(x2)*h1(x4);
    Hijkl(ii+38)=h3(x4)*h1(x1)*h1(x2);
    Hijkl(ii+39)=h3(x4)*h1(x1)*h1(x3);
    Hijkl(ii+40)=h3(x4)*h1(x2)*h1(x3);

    Hijkl(ii+41)=h2(x1)*h1(x2)*h1(x3)*h1(x4);
    Hijkl(ii+42)=h2(x2)*h1(x1)*h1(x3)*h1(x4);
    Hijkl(ii+43)=h2(x3)*h1(x1)*h1(x2)*h1(x4);
    Hijkl(ii+44)=h2(x4)*h1(x1)*h1(x2)*h1(x3);

    Hijkl(ii+45)=h2(x1)*h2(x2)*h1(x3);
    Hijkl(ii+46)=h2(x1)*h2(x2)*h1(x4);
    Hijkl(ii+47)=h2(x1)*h2(x3)*h1(x2);
    Hijkl(ii+48)=h2(x1)*h2(x3)*h1(x4);
    Hijkl(ii+49)=h2(x1)*h2(x4)*h1(x2);
    Hijkl(ii+50)=h2(x1)*h2(x4)*h1(x3);
    Hijkl(ii+51)=h2(x2)*h2(x3)*h1(x1);
    Hijkl(ii+52)=h2(x2)*h2(x3)*h1(x4);
    Hijkl(ii+53)=h2(x2)*h2(x4)*h1(x1);
    Hijkl(ii+54)=h2(x2)*h2(x4)*h1(x3);
    Hijkl(ii+55)=h2(x3)*h2(x4)*h1(x1);
    Hijkl(ii+56)=h2(x3)*h2(x4)*h1(x2);

endfunction

function y=h1(x)
    y = x; 
endfunction

function y=h2(x)

//    y = -1 + x^2; 
    y = 0.5*(3*x.^2 -1); 
//    y = 2.0*x.^2 -1; 

endfunction

function y=h3(x)

//    y = -3*x + x^3; 
    y = 0.5*(5*x.^3 - 3*x); 
//    y = 4*x.^3 - 3*x; 

endfunction

function y=h4(x)

//    y = 3 - 6*x^2 + x^4; 
    y = 0.125*(35*x.^4 - 30*x.^2 + 3); 
//    y = 8*x.^4 - 8*x.^2 + 1; 

endfunction

function y=h5(x)

//    y = 15*x - 10*x^3 + x^5; 
    y = 0.125*(63*x.^5 - 70*x.^3 + 15*x ); 
//    y = 16*x.^5 - 20*x.^3 + 5*x; 

endfunction

function y=h6(x)

//    y = 15*x - 10*x^3 + x^5; 
    y = 0.0625*(231*x.^6 - 315*x.^4 + 105*x.^2 - 5 ); 
//    y = 32*x.^6 - 48*x.^4 + 18*x.^2 - 1; 

endfunction

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
    
    dfds = max(min(dfds,20),-20);

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

function [G,b,s] = circulation(n,rho,u,L,ar,cl_avg,A1,A2,A3,A4,A5,A6)

    q    = 0.5 * rho * u * u;
    area = L / (q * cl_avg); 
    b    = sqrt(ar*area); 
    G0 = 4 * L / (%pi * rho * u * b); 

//    s = ([0:1:n] + 0.5) * b / (n+1) - b/2; 
//    ds = b / (n+1); 

    tt=0:%pi/n:%pi;
    s = -b/2 * cos(tt);


//    G = G0 * sqrt( 1 - (2*s/b).^2 ) + A1 * cos(%pi * 2*s/b) + A2 * cos(2 * %pi * 2*s/b) + A3 * cos(4 * %pi * 2*s/b) + A4 * (cos(1 * %pi * 1 *s/b) - 0.5);

    G = G0 * sqrt( 1 - (2*s/b).^2 ) + A1 * h1(2*s/b) + A2 * h2(2*s/b) + A3 * h3(2*s/b) + A4 * h4(2*s/b) + A5 * h5(2*s/b) + A6 * h6(2*s/b);
    
//    G = G0 * ones(s) + A1 * h1(2*s/b) + A2 * h2(2*s/b) + A3 * h3(2*s/b) + A4 * h4(2*s/b) + A5 * h5(2*s/b) + A6 * h6(2*s/b);; 

    G(1) = 0; 
    G($) = 0; 

endfunction

function w=inducedvelocity(n,b,G)

//    s = ([0:1:n] + 0.5) * b / (n+1) - b/2; 
//    ds = b / (n+1); 

    tt=0:%pi/n:%pi;
    s = -b/2 * cos(tt);

//    dGds_aux  = ([G G($) G($)] - [G(1) G(1) G]) / (2*ds); 
//    dGds = dGds_aux(2:$-1); 
//    dGds(1) = dGds(2); 
//    dGds($) = dGds($-1); 

    dGds = fdd(G,s); 
    
//    disp([size(dGds) size(G)])

    eps=0.0001;

    for j = 1:n+1; 
    
        y = s(j); 
        w_s = dGds .* (-y + s) ./((-y + s).^2 + eps) * 1/( 4*%pi );
//        w(j)= (sum(w_s) - 0.5*(w_s(1) + w_s($)))* ds; 
        w(j) = ti(w_s,s);

    end

//    scf(f4); clf; plot(s,dGds,'b.-');

endfunction;
    
function w=wakevelocity(n,b,bL,G,u,x0,y0,z0)

// bL = b lead A/C

// x, y, z = non dimensional distance between two A/C 

x1 = x0 * b; 
y1 = y0 * b; 
z1 = z0 * b; 

nL = length(G) - 1; 
//sL = ([0:1:nL] + 0.5) * bL / (nL+1) - bL/2; 
//dsL = bL / (nL+1); 

tt=0:%pi/nL:%pi;
sL = -bL/2 * cos(tt);


// s = ([0:1:n] + 0.5) * b / (n+1) - b/2; 
// ds = b / (n+1); 

tt=0:%pi/n:%pi;
s = -b/2 * cos(tt);



/*dGds_aux  = ([G G($) G($)] - [G(1) G(1) G]) / (2*ds); 
dGds = dGds_aux(2:$-1); 
dGds(1) = dGds(2); 
dGds($) = dGds($-1); 
*/

dGds = fdd(G,s); 


tau = -x1 / u;
rc  = 0.52 * sqrt(tau);

for j = 1:n+1; 
    
    y = s(j) + y1; 

    w_s = dGds .* (- y + sL) ./ ((y - sL).^2 + z1^2 + rc^2) .* (1 - x1 *ones(sL)./sqrt(x1^2 + (y - sL).^2 + z1^2) ) * 1/( 4*%pi );

//    w(j)= (sum(w_s) - 0.5*(w_s(1) + w_s($))) * dsL; 

    w(j) = ti(w_s,sL);

end
    
endfunction;


function f = induced_drag(x)

    lateral_distance = x(1); 
    long_distance    = x(2);
    ar               = x(3); 
    cl_avg           = x(4); 
    u                = x(5); 

    A1_LA            = 0.0;
    A2_LA            = x(6);
    A3_LA            = 0.0;
    A4_LA            = x(7);
    A5_LA            = 0.0; 
    A6_LA            = x(8); 
    A1               = 0.0;
    A2               = x(9);
    A3               = 0.0;
    A4               = x(10);
    A5               = 0.0;
    A6               = x(11);

/*    ii=1;  
    A1_LA            = x(ii+2);
    A2_LA            = x(ii+3);
    A3_LA            = x(ii+4);
    A4_LA            = x(ii+5);
    A5_LA            = x(ii+6); 
    A6_LA            = x(ii+7); 
//    B3_LA            = 0.0; 
    A1               = x(ii+8);
    A2               = x(ii+9);
    A3               = x(ii+10);
    A4               = x(ii+11);
    A5               = x(ii+12);
    A6               = x(ii+13);
*/    
/*    B1               = x(10); 
    B2               = x(11); 
    B3               = x(12); 
*/


 
    [G,b,s] = circulation(n,rho,u,L,ar,cl_avg,A1,A2,A3,A4,A5,A6);
    ww    = inducedvelocity(n,b,G);
    ww_2  = ww($:-1:1);
    [G_LA,b,s] = circulation(n,rho,u,L_LA,ar,cl_avg,A1_LA,A2_LA,A3_LA,A4_LA,A5_LA,A6_LA);
    wwL2T   = wakevelocity(n,b,b,G_LA,u,-long_distance,lateral_distance,0);
    wwL2T_2 = wwL2T($:-1:1);
    wwT2L   = wakevelocity(n,b,b,G,u,+long_distance,-lateral_distance,0);
    wwT2L_2 = wwT2L($:-1:1); 
    ww_LA = inducedvelocity(n,b,G_LA);

    wwT2T = wakevelocity(n,b,b,G,u,-long_distance,lateral_distance,0);


// chord based on elliptic planform

    q    = 0.5 * rho * u * u;
    area = L / (q * cl_avg); 
    chord0 = ( area * 4) / (%pi * b);
    chord = chord0 * sqrt( 1 - (2*s/b).^2 ); 
    re    = rho * chord * u /nu; 

   for ispan=1:n+1; 
       Cd0(ispan) = linear_interpn(re(ispan),cl_avg,Re,Cl,v); 
   end; 

   fi = rho * u * u * u * Cd0' .* chord; 

//   disp([size(fi) size(s)])
   pwr_req_visc = ti(fi,s); 

   
// solo for ref
//    fi = -2*(rho * G' .* ww);
//    di_solo =  4* ti(fi,s); // -2*sum(G * rho * ww)*ds;


//    disp([size(G) size(ww) size(G_LA) size(ww_LA)])

//    f = -2*sum(rho * G' .* (ww + ww2))*ds - sum(rho * G_LA' .* ww_LA)*ds;
//    f = -1*sum(rho * G' .* ww ) * ds;

// chevron 4 drones     
    fi = -(rho * G' .* (ww + wwL2T)) - (rho * G' .* (ww_2 + wwL2T_2)) - (rho * G' .* (ww + wwT2T)) - (rho * G_LA' .* (ww_LA + wwT2L + wwT2L_2));    
//    fi = -2*(rho * G' .* (ww + wwL2T)) - (rho * G_LA' .* (ww_LA + wwT2L + wwT2L_2));    
//    fi = -2*(rho * G' .* (ww + wwL2T)) - (rho * G_LA' .* ww_LA)    

//    fi = -rho * G' .* ww; 
//   disp(size(fi))

    pwr_req = ti(fi',s) * u + pwr_req_visc;     
    f(1,1) = pwr_req * 12 / power_density; 
    
//    disp([f(1,1) pwr_req_visc])

//    f(1,1) = di_solo; 

// echelon 4 drones 
//    fi = -(rho * G' .* (ww + wwL2T)) - 2*(rho * G' .* (ww + wwT2T)) - (rho * G_LA' .* (ww_LA + wwT2L));    


// wig root bending moment 

    fi = rho * u * G .* abs(s) * 0.5;
//    disp([s(1) s($)])
    f(1,2) =  ti(fi,s) / mb0_ref; // -2*sum(G * rho * ww)*ds;






//    disp([x(1) x(2) x(5) f wwT2L(1) wwL2T(1)]);

//    scf(f1); clf; plot(s,ww,s,wwL2T,'r-',s,wwT2L,'g-'); 
//    scf(f2); clf; plot(s,ww_LA);
//    scf(f3); clf; plot(s,G); ylabel('$ G $'); 

endfunction; 

// cd data 
    
// cl = 0.90; 
aa = [5e5 0.00975; 
     3e5 0.015; 
     2e5 0.0150; //
     1.5e5 0.01525; //
     1e5 0.01575]; 
v(:,6) = aa($:-1:1,2);      
     
// cl = 0.80; 
aa = [5e5 0.0095; 
     3e5 0.01175; 
     2e5 0.0150; //
     1.5e5 0.01525; //
     1e5 0.01575]; 
v(:,5) = aa($:-1:1,2);      

// cl = 0.70; 
aa = [5e5 0.0075; 
     3e5 0.0085;
     2e5 0.0095; 
     1.5e5 0.0130
     1e5 0.0155]; 
v(:,4) = aa($:-1:1,2);      
     
// cl = 0.60; 
aa = [5e5 0.0060 
     3e5 0.0075
     2e5 0.0087
     1.5e5 0.0120
     1e5 0.0150]; 
v(:,3) = aa($:-1:1,2);      


// cl = 0.50; // redo qith ext'd AoA range
aa = [5e5 0.0062 
     3e5 0.0081
     2e5 0.0091
     1.5e5 0.0111
     1e5 0.0132]; 
v(:,2) = aa($:-1:1,2);      

v(:,1) = v(:,2); 
v(:,7) = v(:,6); 

v = [v(1,:); 
     v; 
     v($,:)]; 

Re=[1e4 1e5 1.5e5 2e5 3e5 5e5 1e6];
Cl=[0.0 0.5 0.6 0.7 0.8 0.9 1.2]; 

// 

n = 151; 
rho = 0.08; 
// u = 20.0; 
//    L = 750.0; 
L_LA = 850; 
L    = 850; 
// ar = 10; 
// cl_avg = 0.50; 
// cla = 2*%pi; 
di_ref = 4 * 0.5 * rho * (15 * 15) * (1.15^2) / (%pi * 30) * 25;
mb0_ref = L/2 * (15) / 2; 
nu = 1.5e-5; 
power_density = 435.; // Wh/kg


PopSize     = 500;
Proba_cross = 0.5;
Proba_mut   = 0.3;
NbGen       = 4;
NbCouples   = 110;
Log         = %T;
nb_disp     = 10; // Nb point to display from the optimal population
pressure    = 0.1;

ga_params = init_param();
ga_params = add_param(ga_params,'dimension',2);
ga_params = add_param(ga_params,'minbound',[0.5; 0.25; 10; 0.25; 15; -2.0*ones(6,1)]); //-0.5*ones(12,1)]);
ga_params = add_param(ga_params,'maxbound',[2.5;  25;  18; 1.25; 35; +2.0*ones(6,1)]); // +0.5*ones(12,1)]);

[pop_opt, fobj_pop_opt, pop_init, fobj_pop_init] = optim_nsga2(induced_drag, PopSize,NbGen, Proba_mut, Proba_cross, Log, ga_params)

[F_out,X_out,Ind_out] = pareto_filter(fobj_pop_opt,pop_opt);


for ii=1:PopSize; 
    a=pop_opt(ii); 
    lateral(ii)=a(1); 
    long(ii) =a(2);
    ar(ii)   =a(3);
    cl(ii)   =a(4);
    u(ii)    =a(5); 
//    A1LA(ii)=a(3);
    A2LA(ii)=a(6);
//    A3LA(ii)=a(5);
    A4LA(ii)=a(7);
//    A5LA(ii)=a(7);
    A6LA(ii)=a(8);
//    A1TA(ii)=a(9);
    A2TA(ii)=a(9);
//    A3TA(ii)=a(11);
    A4TA(ii)=a(10);
//    A5TA(ii)=a(13);
    A6TA(ii)=a(11);
    
end


f1=scf(1); clf; 
plot(fobj_pop_opt(:,1),fobj_pop_opt(:,2),'o'); 
plot(F_out(:,1),F_out(:,2),'r.');
xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE \mu_2$');

f2=scf(2); clf; 
f=get("current_figure") //get the handle of the current figure :
f.figure_size=[600,900]

subplot(511);
plot(fobj_pop_opt(:,1),lateral,'o'); 
//xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE Lat. $');

// f3=scf(3); 
subplot(512);
plot(fobj_pop_opt(:,1),long,'o'); 
//xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE Long. $');

//f4=scf(4); 
subplot(513);
plot(fobj_pop_opt(:,1),ar,'o'); 
//xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE \lambda $');

//f5=scf(5); 
subplot(514);
plot(fobj_pop_opt(:,1),cl,'o'); 
//xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE C_l $');

subplot(515);
plot(fobj_pop_opt(:,1),u,'o'); 
xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE u $');


f3=scf(3); clf; 
f=get("current_figure") //get the handle of the current figure :
f.figure_size=[600,900]

subplot(511);
plot(fobj_pop_opt(:,2),lateral,'o'); 
//xlabel('$\LARGE \mu_2$');
ylabel('$\LARGE Lat. $');

subplot(512);
plot(fobj_pop_opt(:,2),long,'o'); 
//xlabel('$\LARGE \mu_2$');
ylabel('$\LARGE Long. $');

subplot(513);
plot(fobj_pop_opt(:,2),ar,'o'); 
//xlabel('$\LARGE \mu_2$');
ylabel('$\LARGE \lambda $');

subplot(514);
plot(fobj_pop_opt(:,2),cl,'o'); 
//xlabel('$\LARGE \mu_2$');
ylabel('$\LARGE C_l $');


subplot(515);
plot(fobj_pop_opt(:,2),u,'o'); 
xlabel('$\LARGE \mu_2$');
ylabel('$\LARGE u $');


f4=scf(4); clf; 
//plot(fobj_pop_opt(:,1),A1LA,'ro'); 
plot(fobj_pop_opt(:,1),A2LA,'bo'); 
//plot(fobj_pop_opt(:,1),A3LA,'go'); 
plot(fobj_pop_opt(:,1),A4LA,'r*'); 
//plot(fobj_pop_opt(:,1),A5LA,'b*'); 
plot(fobj_pop_opt(:,1),A6LA,'g*'); 
xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE Leader $');

f5=scf(5); clf; 
//plot(fobj_pop_opt(:,1),A1TA,'ro'); 
plot(fobj_pop_opt(:,1),A2TA,'bo'); 
//plot(fobj_pop_opt(:,1),A3TA,'go'); 
plot(fobj_pop_opt(:,1),A4TA,'r*'); 
//plot(fobj_pop_opt(:,1),A5TA,'b*'); 
plot(fobj_pop_opt(:,1),A6TA,'g*'); 
xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE TRAILER $');





