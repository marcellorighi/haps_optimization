
function f=deb_1(x)
   f1_x1 = x(1);
   g_x2  = 1 + 9 * sum((x(2:$)-x(1)).^2) / (length(x) - 1);
   h     = 1 - sqrt(f1_x1 / g_x2);

   f(1,1) = f1_x1;
   f(1,2) = g_x2 * h;
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
ga_params = add_param(ga_params,'minbound',zeros(2,1));
ga_params = add_param(ga_params,'maxbound',ones(2,1));

[pop_opt, fobj_pop_opt, pop_init, fobj_pop_init] = optim_nsga2(deb_1, PopSize,NbGen, Proba_mut, Proba_cross, Log, ga_params)
