
clear; 
 
cl = 0.90; 
a = [5e5 0.00975; 
     3e5 0.015; 
     2e5 0.0150; //
     1.5e5 0.01525; //
     1e5 0.01575]; 
v(:,6) = a($:-1:1,2);      
     
cl = 0.80; 
a = [5e5 0.0095; 
     3e5 0.01175; 
     2e5 0.0150; //
     1.5e5 0.01525; //
     1e5 0.01575]; 
v(:,5) = a($:-1:1,2);      

cl = 0.70; 
a = [5e5 0.0075; 
     3e5 0.0085;
     2e5 0.0095; 
     1.5e5 0.0130
     1e5 0.0155]; 
v(:,4) = a($:-1:1,2);      
     
cl = 0.60; 
a = [5e5 0.0060 
     3e5 0.0075
     2e5 0.0087
     1.5e5 0.0120
     1e5 0.0150]; 
v(:,3) = a($:-1:1,2);      


cl = 0.50; // redo qith ext'd AoA range
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
cl=[0.0 0.5 0.6 0.7 0.8 0.9 1.2]; 







