function y=banana(x)

// DV: x(1) = span, x(2) = speed

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

  y = cl / Cdwing; 
  
  
endfunction

// opt = optimset ( "PlotFcns" , optimplotfval );
// [x fval] = fminsearch ( banana , [15 40] , opt );

// Upper and lower bounds on x
x0 = [15; 30];
[fopt, xopt, gopt] = optim(banana, "b", [10; 20], [20; 40], x0)


