

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

function w=wakevelocity(G,dGds,s,x0,y0,z0,b);

x1 = x0 * b; 
y1 = y0 * b; 
z1 = z0 * b; 

// speed hard wired
tau = -x1 / 20;
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
    Amle             = x((ne+6-2)+1:(ne+6-2)+ne-1);
    Ale              =[A1; Amle]; 
    Ambm             = x((ne+6-2)+ne:(ne+6-2)+ne+ne-2);
    Abm              =[A1; Ambm]; 
    
 // chord based on elliptic planform

    q    = 0.5 * rho * u * u;
    area = L / (q * cl_avg); 
    b    = sqrt(ar * area); 
    s     = -(b/2)*cos(theta);
    chord0 = ( area * 4) / (%pi * b);
    chord = chord0 * sqrt( 1 - (2*s/b).^2 ); 
    
    re    = rho * chord * u /nu; 
    
// circulation around trailer drones wings 
    [G,dGds,w] = circulation(A,u,b); 
    cl_local = 2 * G ./ (chord * u); 

// circulation around leader drone wings 
    [Gle,dGdsle,wle] = circulation(Ale,u,b); 
    cl_local_le = 2 * Gle ./ (chord * u); 

// induced velocity on trailer 
    w_l2t=wakevelocity(Gle,dGdsle,s,long_distance,lateral_distance,0.0,b);

//    disp([rho*u*ti(G,s) ])

// induced power leader 

    fi = (rho * Gle .* (wle )); 
    pwr_req_ind_le = ti(fi,s) * u;     

// viscous drag leader

    for ispan=1:n; 
        Cd0(ispan) = linear_interpn(re(ispan),cl_local_le(ispan),Re,Cl,v); 
    end; 

    fi = q * u * Cd0' .* chord; 
    pwr_req_visc_le = ti(fi,s); 

// overall power required 
    pwr_req_le = pwr_req_visc_le + pwr_req_ind_le; 
    
// induced power trailer 

//    disp([size(w) size(w_l2t)])
    fi = (rho * G .* (w - w_l2t')); 
    pwr_req_ind = ti(fi,s) * u;     

// viscous drag trailer

    for ispan=1:n; 
        Cd0(ispan) = linear_interpn(re(ispan),cl_local(ispan),Re,Cl,v); 
    end; 

    fi = q * u * Cd0' .* chord; 

    pwr_req_visc = ti(fi,s); 
    
// overall power required trailer 

    pwr_req =  pwr_req_ind + pwr_req_visc;   

// objectives 
    f(1,1) = pwr_req; // * 12 / power_density; 
    
// wing root bending moment 

    [G,dGds,w] = circulation(Abm,u,b); 
    fi = rho * u * G .* abs(s) * 0.5;
    bending_m = ti(fi,s); 
    f(1,2) =  1.00/rho*1.35/cl_avg*bending_m; // pwr_req_le; //ti(fi,s);

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
ww = wakevelocity(G,dGds,s,-1,0.0,0.0,b);

/*
f11 = scf(11); clf; 
plot(s,w,s,ww,'b-'); 
*/






// optimization 

PopSize     = 400;
Proba_cross = 0.5;
Proba_mut   = 0.3;
NbGen       = 2; //48;
NbCouples   = 110;
Log         = %T;
nb_disp     = 10; // Nb point to display from the optimal population
pressure    = 0.1;


// 

ga_params = init_param();
ga_params = add_param(ga_params,'dimension',2);
ga_params = add_param(ga_params,'minbound',[0.5;   -10.0;   10; 0.30; 15; -0.0025*ones(ne-1,1); -0.0025*ones(ne-1,1); -0.0025*ones(ne-1,1)]); //-0.5*ones(12,1)]);
ga_params = add_param(ga_params,'maxbound',[1.5;    -0.15;  12; 1.25; 46; +0.0025*ones(ne-1,1); +0.0025*ones(ne-1,1); +0.0025*ones(ne-1,1)]); // +0.5*ones(12,1)]);

[pop_opt, fobj_pop_opt, pop_init, fobj_pop_init] = optim_nsga2(induced_drag, PopSize,NbGen, Proba_mut, Proba_cross, Log, ga_params)

[F_out,X_out,Ind_out] = pareto_filter(fobj_pop_opt,pop_opt);

PF_size = length(X_out);

AA = zeros(ne-1,PF_size);

disp('UQ: NIPC'); 

for ii=1:PF_size; 
    disp(ii); 
    a=X_out(ii); 
    lateral(ii)=a(1); 
    long(ii) =a(2);
    ar(ii)   =a(3);
    cl(ii)   =a(4);
    u(ii)    =a(5);
    area(ii) = L / (0.5 * rho * u(ii) * u(ii) * cl(ii) ); 
    bb(ii)   = sqrt(ar(ii)*area(ii) ); 
    for jj=2:ne;
        AA(jj,ii) = a(5+jj-1);
    end;  

    AA(1,ii) = cl(ii)/(%pi * ar(ii)); 

    for jj=2:ne;
        AAle(jj,ii) = a(ne+5+jj-2);
    end;  
    
    for jj=2:ne;
        AAbm(jj,ii) = a(ne+ne+5+jj-3);
    end;  
    
    AAle(1,ii) = cl(ii)/(%pi * ar(ii)); 
    AAbm(1,ii) = AAle(1,ii); 
    
    [Gpf(ii,:),dGdspf(ii,:),wpf(ii,:)] = circulation(AA(:,ii),u(ii),bb(ii)); 
    [Glepf(ii,:),dGdslepf(ii,:),wlepf(ii,:)] = circulation(AAle(:,ii),u(ii),bb(ii)); 

    s     = -(bb(ii)/2)*cos(theta);
    w_l2t=wakevelocity(Glepf(ii,:),dGdslepf(ii,:),s,long(ii),lateral(ii),0.0,bb(ii));

    chord0 = ( area(ii) * 4) / (%pi * bb(ii));
    chord = chord0 * sqrt( 1 - (2*s/bb(ii)).^2 ); 
    re    = rho * chord * u(ii) /nu; 
    cl_local = 2 * Gpf(ii,:) ./ (chord * u(ii)); 

    for ispan=1:n; 
        Cd0(ispan) = linear_interpn(re(ispan),cl_local(ispan),Re,Cl,v); 
    end; 

    fi = 0.5 * rho * u(ii) * u(ii) * u(ii) * Cd0' .* chord; 
    pwr_req_visc(ii) = ti(fi,s); 

    fi = (rho * Gpf(ii,:) .* (wpf(ii,:) - w_l2t')); 
    pwr_req_ind(ii) = ti(fi,s) * u(ii);    

    fi = rho * u(ii) * Gpf(ii,:) .* abs(s) * 0.5;
    bending(ii) = ti(fi,s); 
    
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


f3 = scf(3); clf; 
f=get("current_figure") //get the handle of the current figure :
f.figure_size=[600,900]

subplot(511);
plot(F_out(:,1),AA(2,:),'o'); 
//xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE AA(2) $');

// f3=scf(3); 
subplot(512);
plot(F_out(:,1),AA(3,:),'o'); 
//xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE AA(3) $');

//f4=scf(4); 
subplot(513);
plot(F_out(:,1),AA(4,:),'o'); 
//xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE AA(4) $');

//f5=scf(5); 
subplot(514);
plot(F_out(:,1),AA(5,:),'o'); 
//xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE AA(5) $');

subplot(515);
plot(F_out(:,1),AA(6,:),'o'); 
xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE AA(6) $');

f4 = scf(4); clf; 
f=get("current_figure") //get the handle of the current figure :
f.figure_size=[600,900]

subplot(511);
plot(F_out(:,2),AA(2,:),'o'); 
//xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE AA(2) $');

// f3=scf(3); 
subplot(512);
plot(F_out(:,2),AAbm(3,:),'o'); 
//xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE AA(3) $');

//f4=scf(4); 
subplot(513);
plot(F_out(:,2),AAbm(4,:),'o'); 
//xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE AA(4) $');

//f5=scf(5); 
subplot(514);
plot(F_out(:,2),AAbm(5,:),'o'); 
//xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE AA(5) $');

subplot(515);
plot(F_out(:,2),AAbm(6,:),'o'); 
xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE AA(6) $');

f5=scf(5); clf; 
f=get("current_figure") //get the handle of the current figure :
f.figure_size=[600,900]

subplot(611);
plot(F_out(:,1),lateral,'o'); 
//xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE Lat. $');

// f3=scf(3); 
subplot(612);
//plot(F_out(:,1),lift_trailers,'o'); 
plot(F_out(:,1),long,'o'); 
//xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE L_T. $');

//f4=scf(4); 
subplot(613);
plot(F_out(:,1),bb,'o'); 
//xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE b $');

//f5=scf(5); 
subplot(614);
plot(F_out(:,1),cl,'o'); 
//xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE C_l $');

subplot(615);
plot(F_out(:,1),u,'o'); 
xlabel('$ \mu_1$');
ylabel('$\LARGE u $');

subplot(616);
plot(F_out(:,2),area,'o'); 
xlabel('$ \mu_1$');
ylabel('$\LARGE area $');

f6=scf(6); clf; 
f=get("current_figure") //get the handle of the current figure :
f.figure_size=[600,900]

subplot(611);
plot(F_out(:,2),lateral,'o'); 
//xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE Lat. $');

// f3=scf(3); 
subplot(612);
//plot(F_out(:,1),lift_trailers,'o'); 
plot(F_out(:,2),long,'o'); 
//xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE L_T. $');

//f4=scf(4); 
subplot(613);
plot(F_out(:,2),bb,'o'); 
//xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE b $');

//f5=scf(5); 
subplot(614);
plot(F_out(:,2),cl,'o'); 
//xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE C_l $');

subplot(615);
plot(F_out(:,2),u,'o'); 
xlabel('$ \mu_1$');
ylabel('$\LARGE u $');

subplot(616);
plot(F_out(:,2),area,'o'); 
xlabel('$ \mu_1$');
ylabel('$\LARGE area $');




f7 = scf(7); clf; 
f8 = scf(8); clf; 

for ii= 1:PF_size; 

    scf(f7); 
    plot(s/b,Gpf(ii,:)');
    
    scf(f8); 
    plot(s/b,wpf(ii,:)'); 
    
end


f9=scf(9); clf; 
f=get("current_figure") //get the handle of the current figure :
f.figure_size=[600,900];
subplot(311); 
plot(F_out(:,1),pwr_req_ind,'ro',F_out(:,1),pwr_req_visc,'bo',F_out(:,1),(pwr_req_ind + pwr_req_visc),'g.');
legend('$\Large IND$','$\Large VISC$','$\Large TOT$'); 
xlabel('$\Large \mu_1$'); 
subplot(312); 
plot(ar,pwr_req_ind,'ro',ar,pwr_req_visc,'bo',ar,(pwr_req_ind + pwr_req_visc),'g.');
legend('$\Large IND$','$\Large VISC$','$\Large TOT$'); 
xlabel('$\Large \lambda$'); 
subplot(313); 
plot(cl,pwr_req_ind,'ro',cl,pwr_req_visc,'bo',cl,(pwr_req_ind + pwr_req_visc),'g.');
legend('$\Large IND$','$\Large VISC$','$\Large TOT$'); 
xlabel('$\Large C_l$'); 

// UQ 

exec('/home/rigm/SU2/runs/Scilab/airfoil/opt/multi_variate_polynomials.sce'); 

nactive = 4; 
nsamples = 102; 
eps = 0.050; 

for ii=1:PF_size; 
  disp([ii PF_size]);  
  RV = grand(nactive,"prm",(0:(nsamples-1))')/(nsamples) + grand(nsamples,nactive,"unf",0,1/nsamples)-0.5;
  for isample = 1:nsamples;       
      x(1) = lateral(ii) * (1 + eps * RV(isample,1)); 
      x(2) = long(ii); // * (1 + eps * RV(isample,2)); 
      x(3) = ar(ii) * (1 + eps * RV(isample,2)); 
      x(4) = cl(ii) * (1 + eps * RV(isample,3)); 
      x(5) = u(ii) * (1 + eps * RV(isample,4)); 
      
      for jj=2:ne;
          x(5+jj-1) = AA(jj,ii);
      end;  

      for jj=2:ne;
          x(ne+5+jj-2) = AAle(jj,ii);
      end;  
    
      for jj=2:ne;
          x(ne+ne+5+jj-3) = AAbm(jj,ii);
      end;  
      
      fsample(:,isample) = induced_drag(x);

      A5(isample,:) = Psi(RV(isample,1),RV(isample,2),RV(isample,3),RV(isample,4))'; 

  end; 

  A4 = A5(:,1:70); 
  A3 = A4(:,1:35); 
  A2 = A4(:,1:15); 
  A1 = A4(:,1:5); 

  aa1=inv(A1'*A1)*(A1'*fsample(1,:)'); 
  aa2=inv(A2'*A2)*(A2'*fsample(1,:)'); 
  aa3=inv(A3'*A3)*(A3'*fsample(1,:)'); 
  aa4=inv(A4'*A4)*(A4'*fsample(1,:)'); 

  bb1=inv(A1'*A1)*(A1'*fsample(2,:)'); 
  bb2=inv(A2'*A2)*(A2'*fsample(2,:)'); 
  bb3=inv(A3'*A3)*(A3'*fsample(2,:)'); 
  bb4=inv(A4'*A4)*(A4'*fsample(2,:)'); 

  mu_a1 = aa1(1); 
  mu_a2 = aa2(1); 
  mu_a3 = aa3(1); 
  mu_a4 = aa4(1); 
   
  mu_b1 = bb1(1); 
  mu_b2 = bb2(1); 
  mu_b3 = bb3(1); 
  mu_b4 = bb4(1); 

  mu_f1(ii) = mu_a2; 
  mu_f2(ii) = mu_b2; 

  sigma_a1 = sqrt(((aa1'.^2).*(ones(1,nsamples)*(A1.*A1)))*[0; ones(4,1)]); 
  sigma_a2 = sqrt(((aa2'.^2).*(ones(1,nsamples)*(A2.*A2)))*[0; ones(14,1)]); 
  sigma_a3 = sqrt(((aa3'.^2).*(ones(1,nsamples)*(A3.*A3)))*[0; ones(34,1)]); 
  sigma_a4 = sqrt(((aa4'.^2).*(ones(1,nsamples)*(A4.*A4)))*[0; ones(69,1)]); 

  sigma_b1 = sqrt(((bb1'.^2).*(ones(1,nsamples)*(A1.*A1)))*[0; ones(4,1)]); 
  sigma_b2 = sqrt(((bb2'.^2).*(ones(1,nsamples)*(A2.*A2)))*[0; ones(14,1)]); 
  sigma_b3 = sqrt(((bb3'.^2).*(ones(1,nsamples)*(A3.*A3)))*[0; ones(34,1)]); 
  sigma_b4 = sqrt(((bb4'.^2).*(ones(1,nsamples)*(A4.*A4)))*[0; ones(69,1)]); 

  sigma_f1(ii) = sigma_a2; 
  sigma_f2(ii) = sigma_b2; 

  sigma_f1_check(ii) = sigma_a4; 
  sigma_f2_check(ii) = sigma_b4; 

end; 

f20 = scf(20); clf; 
subplot(211); 
plot(mu_f1,sigma_f1,'ro');
plot(mu_f1,sigma_f1_check,'b.');
xlabel('$\Large \mu_1$');
ylabel('$\Large \sigma_1$');
legend('$\Large 2^{nd}\,order$','$\Large 4^{th}\,order$');

subplot(212); 
plot(mu_f2,sigma_f2,'ro');
plot(mu_f2,sigma_f2_check,'b.');
xlabel('$\Large \mu_2$');
ylabel('$\Large \sigma_2$');
legend('$\Large 2^{nd}\,order$','$\Large 4^{th}\,order$');

f21 = scf(21); clf; 
subplot(211); 
plot(ar,sigma_f1,'ro');
plot(ar,sigma_f1_check,'b.');
xlabel('$\Large \lambda$');
ylabel('$\Large \sigma_1$');
legend('$\Large 2^{nd}\,order$','$\Large 4^{th}\,order$');

subplot(212); 
plot(cl,sigma_f2,'ro');
plot(cl,sigma_f2_check,'b.');
xlabel('$\Large C_l$');
ylabel('$\Large \sigma_2$');
legend('$\Large 2^{nd}\,order$','$\Large 4^{th}\,order$');









