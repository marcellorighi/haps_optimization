
function f=induced_power(x)
   cl   = x(1); 
   u    = x(2);
   ar   = x(3);  
   
   Lift = 1500.; 
   rho = 0.08; // kg/m^(-3)
   nu = 1.5e-5; 
   powerdensity = 435.; // Wh/kg
   q = 0.5 * rho * u^2;    
   Area = Lift / (cl * q); 
   span = sqrt( ar * Area); 
   chord = Area / span; 
   e = 1 / (1.05 + 0.007*%pi*ar)
   cdi = cl^2 /(%pi * e * ar); // disp(size(Cd0)); disp(size(c))

  pwrreq =  u * q * span * cdi; 
  energy_1night = pwrreq * 12; // Wh
  battery_mass = energy_1night / powerdensity; 
  
  f(1,1) = battery_mass; 

  clmax  = 1.50; 
  f(1,2) = clmax/cl * (1.225/rho) * (10/u)^2;  

endfunction

PopSize     = 200;
Proba_cross = 0.5;
Proba_mut   = 0.3;
NbGen       = 8;
NbCouples   = 110;
Log         = %T;
nb_disp     = 10; // Nb point to display from the optimal population
pressure    = 0.1;

ga_params = init_param();
ga_params = add_param(ga_params,'dimension',2);
ga_params = add_param(ga_params,'minbound',[0.40; 10; 10]);
ga_params = add_param(ga_params,'maxbound',[1.20; 25; 35]);

[pop_opt, fobj_pop_opt, pop_init, fobj_pop_init] = optim_nsga2(induced_power, PopSize,NbGen, Proba_mut, Proba_cross, Log, ga_params)

f1=scf(1); 
plot(fobj_pop_opt(:,1),fobj_pop_opt(:,2),'o'); 
xlabel('$\LARGE \mu_1$');
ylabel('$\LARGE \mu_2$');

// ???? 

for ii=1:PopSize; 
    a=pop_opt(ii); 
    cl(ii)=a(1); 
    u(ii) =a(2);
    ar(ii)=a(3);
end



