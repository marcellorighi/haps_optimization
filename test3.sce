
function f=powerrequired(x)

   f = 1+x(1)^2 + x(2)^2; //)(x(2)-x(1)^2)^2 + (1-x(1))^2;

   span = x(1); 
   u = x(2); 
   
   Lift = 500.; 
   Area = 22; // m 
   rho = 0.08; // kg/m^(-3)
   nu = 1.5e-5; 
   
// dyn pressure 

   q = 0.5 * rho * x(2)^2;    

   cl = Lift / (q * Area); 
      
// root chord 

   chord = 4 * Area / (%pi * x(1) ); 

// calc aspect ratio 

   ar = x(1)^2 / Area; 

// loop over n wing sections 

   n = 10; 
   yy = linspace(0,span,n); 
   c = chord*sqrt( 1 - (yy/span).^2); 
   Re = rho * c * u / nu; 
   Cd = 0.010 - 0.006*tanh((Re-4e5)/3e5); 
   Cdi = cl^2 /(%pi * ar); 
   
   Cdwing = sum(c.*(Cd + Cdi)); 
//   Cdwing = 0.08 + Cdi; 

   disp(cl)
   
// assess: AoA, Re for each section 


// calc Cl and Cd for each section incl induced Cd?


// calc L/D??? 

//   y = q * Area * Cdwing;

//  f = sum(Cdi); //Cdwing; 

endfunction

function [f, g, ind]=cost(x, ind)
    f = powerrequired (x);
    g = numderivative (cost, x);
endfunction

x0 = [20 20];
span_min=10; 
span_max=40; 

u_min=12; 
u_max=50; 

x_min=[span_min; u_min]; 
x_max=[span_max; u_max];

// Upper and lower bounds on x
// [fopt, xopt, gopt] = optim(cost, "b", [-1;0;2], [0.5;1;4], x0)

// )[fopt, xopt] = optim (rosenbrockCost, "b", [-2 -2], [-0.5 2], x0)
[fopt, xopt] = optim (cost, "b", x_min, x_max, x0)



