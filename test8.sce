
// rectangular wing

clear; 

function f=powerrequired(x)

//   f = 1+x(1)^2 + x(2)^2; //)(x(2)-x(1)^2)^2 + (1-x(1))^2;
   cl   = x(1); 
   u    = x(2);
   ar   = x(3);  
   
   Lift = 1500.; 

   rho = 0.08; // kg/m^(-3)
   nu = 1.5e-5; 
   powerdensity = 435.; // Wh/kg
   
// dyn pressure 

   q = 0.5 * rho * u^2;    

   Area = Lift / (cl * q); 
     
// span 

   span = sqrt( ar * Area); 

// root chord 

   chord = Area / span; 


// loop over n wing sections 

   n = 20; 
   yy = linspace(0,span/2,n); 
//   yy1 = linspace(0,span,n+1); 
//   yy2 = ([00 yy] + [yy yy($)])/2;
//   yy = yy2(2:$-1);
   c = ones(yy) * chord; //*sqrt( 1 - (2*yy/span).^2); 
   Re = rho * c * u / nu; 
//   Cd0 = 0.013 - 0.010*tanh((Re-2.0e5)/3e5); 

   for ispan=1:n; 
       Cd0(ispan) = linear_interpn(Re(ispan),cl,re,Cl,v); 
   end; 
//   disp(Cd0);
   
   
//   Cd1 = Cd0; 
//   Cd = Cd0 + Cd1 * cl; 
   e = 1 / (1.05 + 0.007*%pi*ar)

   Cdi = cl^2 /(%pi * e * ar); // disp(size(Cd0)); disp(size(c))
   
   Cdwing = sum(c'.*(Cd0 + Cdi))/n; 
//   Cdwing = 0.08 + Cdi; 

//   disp([Area cl sum(Cdi) sum(Cd)])

//   disp([cl u ar Area])   
// assess: AoA, Re for each section 


// calc Cl and Cd for each section incl induced Cd?


// calc L/D??? 

//   y = q * Area * Cdwing;

//   f =  q * Area * Cdwing; 

  pwrreq =  u * q * span * Cdwing; 
  energy_1night = pwrreq * 12; // Wh
  battery_mass = energy_1night / powerdensity; 
  
  f = battery_mass; 
//   disp([f cl u ar Area span sum(Cd) sum(Cdi)])
  disp([f cl u ar Area span chord(1) Re(1)])
   
//   f = -cl/Cdwing; 

endfunction

function [f, g, ind]=cost(x, ind)
    f = powerrequired (x);
    g = numderivative (powerrequired, x);
endfunction

// cl = 0.90; 
a = [5e5 0.00975; 
     3e5 0.015; 
     2e5 0.0150; //
     1.5e5 0.01525; //
     1e5 0.01575]; 
v(:,6) = a($:-1:1,2);      
     
// cl = 0.80; 
a = [5e5 0.0095; 
     3e5 0.01175; 
     2e5 0.0150; //
     1.5e5 0.01525; //
     1e5 0.01575]; 
v(:,5) = a($:-1:1,2);      

// cl = 0.70; 
a = [5e5 0.0075; 
     3e5 0.0085;
     2e5 0.0095; 
     1.5e5 0.0130
     1e5 0.0155]; 
v(:,4) = a($:-1:1,2);      
     
// cl = 0.60; 
a = [5e5 0.0060 
     3e5 0.0075
     2e5 0.0087
     1.5e5 0.0120
     1e5 0.0150]; 
v(:,3) = a($:-1:1,2);      


// cl = 0.50; // redo qith ext'd AoA range
a = [5e5 0.0062 
     3e5 0.0081
     2e5 0.0091
     1.5e5 0.0111
     1e5 0.0132]; 
v(:,2) = a($:-1:1,2);      

v(:,1) = v(:,2); 
v(:,7) = v(:,6); 

v = [v(1,:); 
     v; 
     v($,:)]; 

re=[1e4 1e5 1.5e5 2e5 3e5 5e5 1e6];
Cl=[0.0 0.5 0.6 0.7 0.8 0.9 1.2]; 


x0 = [0.80; 50; 20];

cl_min=0.15; 
cl_max=1.10; 

u_min=19; 
u_max=50; 

ar_min = 5; 
ar_max = 12; 

x_min=[cl_min; u_min; ar_min]; 
x_max=[cl_max; u_max; ar_max];

// Upper and lower bounds on x
// [fopt, xopt, gopt] = optim(cost, "b", [-1;0;2], [0.5;1;4], x0)

// )[fopt, xopt] = optim (rosenbrockCost, "b", [-2 -2], [-0.5 2], x0)
[fopt, xopt] = optim (cost, "b", x_min, x_max, x0)



