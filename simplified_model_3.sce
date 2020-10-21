

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
    A1               = cl_avg / (%pi * ar); 
    Am               = x(6:(ne+6-2));
    A                =[A1; Am]; 
    
 // chord based on elliptic planform

    q    = 0.5 * rho * u * u;
    area = L / (q * cl_avg); 
    b    = sqrt(ar * area); 
    chord0 = ( area * 4) / (%pi * b);
    chord = chord0 * sqrt( 1 - (2*s/b).^2 ); 
    re    = rho * chord * u /nu; 
    
    s     = -(b/2)*cos(theta);

    [G,dGds,w] = circulation(A,u,b); 

//    disp([rho*u*ti(G,s) ])
    
// induced power 
// disp([size(G) size(w) size(s)])
    fi = (rho * G .* (w )); 
    pwr_req = ti(fi,s) * u;     
    f(1,1) = pwr_req; // * 12 / power_density; 
    
// wing root bending moment 

    fi = rho * u * G .* abs(s) * 0.5;
    f(1,2) =  ti(fi,s);

endfunction; 



// numerics
eps=1e-13; 
dfds_max = 10000; 

// flow 
rho = 0.08; 
u = 20; 
nu = 1.5e-5; 
q    = 0.5 * rho * u * u;

// configuration 
L = 750; 
ar =20; 
b=30; 
x0=-b; 
y0=0.0*b; 
z0=0; 
G0 = 4 * L / (%pi * rho * u * b); 
area = b^2 / ar; 
cl_avg = L / (q * area); 

// mesh 
n=250; 
ne    = 6; 
theta = linspace(0,%pi,n+1)+%pi/(n+1)*0.5;
theta = theta(1:$-1); 
nh = [1:ne]';
theta_sintable = sin(nh*theta); 
theta_costable = cos(nh*theta);
theta_nsintable = diag([1:ne])*theta_sintable; 
theta_ncostable = diag([1:ne])*theta_costable; 
s     = -(b/2)*cos(theta); 
A = zeros(ne,1); 
A(1) = cl_avg / (%pi * ar);

// check initial data 

[G,dGds,w] = circulation(A,u,b); 
lift_check = rho * u * ti(G,s);






// optimization 

PopSize     = 400;
Proba_cross = 0.5;
Proba_mut   = 0.3;
NbGen       = 48;
NbCouples   = 110;
Log         = %T;
nb_disp     = 10; // Nb point to display from the optimal population
pressure    = 0.1;


// 

ga_params = init_param();
ga_params = add_param(ga_params,'dimension',2);
ga_params = add_param(ga_params,'minbound',[2.5;  24.25; 15; 0.45; 15; -0.0010*ones(9,1)]); //-0.5*ones(12,1)]);
ga_params = add_param(ga_params,'maxbound',[3.5;  25;    32; 1.25; 25; +0.0010*ones(9,1)]); // +0.5*ones(12,1)]);

[pop_opt, fobj_pop_opt, pop_init, fobj_pop_init] = optim_nsga2(induced_drag, PopSize,NbGen, Proba_mut, Proba_cross, Log, ga_params)

[F_out,X_out,Ind_out] = pareto_filter(fobj_pop_opt,pop_opt);

PF_size = length(X_out);

AA = zeros(ne-1,PF_size);

for ii=1:PF_size; 
    a=X_out(ii); 
    lateral(ii)=a(1); 
    long(ii) =a(2);
    ar(ii)   =a(3);
    cl(ii)   =a(4);
    u(ii)    =a(5);
    area(ii) = L / (rho * u(ii) * u(ii) * cl(ii) ); 
    bb(ii)   = sqrt(ar(ii)*area(ii) ); 
    for jj=2:ne;
        AA(jj,ii) = a(5+jj-1);
    end;  
    AA(1,ii) = cl(ii)/(%pi * ar(ii)); 
    
    [Gpf(ii,:),dGdspf,wpf(ii,:)] = circulation(AA(:,ii),u(ii),bb(ii)); 
    
end


f1=scf(1); clf; 
plot(fobj_pop_opt(:,1),fobj_pop_opt(:,2),'o'); 
plot(F_out(:,1),F_out(:,2),'r.');
xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE \mu_2$');

f2=scf(2); clf; 
subplot(321)
plot(F_out(:,1),AA(2,:)','o'); 

subplot(322)
plot(F_out(:,1),AA(3,:)','o'); 

subplot(323)
plot(F_out(:,1),AA(4,:)','o'); 

subplot(324)
plot(F_out(:,1),AA(5,:)','o'); 

subplot(325)
plot(F_out(:,1),AA(6,:)','o'); 

//subplot(326)
//plot(F_out(:,1),AA(7,:)','o'); 


f7 = scf(7); clf; 
f8 = scf(8); clf; 

for ii= 1:PF_size; 

    scf(f7); 
    plot(s/b,Gpf(ii,:)');
    
    scf(f8); 
    plot(s/b,wpf(ii,:)'); 
    
end
