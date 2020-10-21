function y=banana(x)
  y = x(1)^2 + x(2)^2; //1+100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
endfunction
opt = optimset ( "PlotFcns" , optimplotfval );
[x fval] = fminsearch ( banana , [-1.2 1] , opt );
