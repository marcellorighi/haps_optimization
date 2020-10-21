clear; 

nsamples = 20; 
nactive=4; 

RV=grand(nactive,"prm",(0:(nsamples-1))')/(nsamples) + grand(nsamples,nactive,"unf",0,1/nsamples)-0.5;


t=0.06; 

dt  = 0.025 * RV(:,1) * t; 
P1 = 1.1019*(t + dt).^2;
P2 = 0.30 * (1 + 0.025 * RV(:,2)); 
P3 = 0.5 * (t + dt); 
P4 = -0.2125 * (1 + RV(:,3)*0.025); 
P5 = P1; 
P6 = P2; 
P7 = P4; 
P8 = zeros(P1); 
P9 = 6.3e-4*2 * ones(P1); 
P10= zeros(P1); 
P11= 8*%pi/(180) * (1 + RV(:,4)*0.025); 



/*
p1=1.1019*t^2; //0.005; 
p2=0.30; 
p3=t/2; 
p4=-0.2125; 
p5=p1; 
p6=p2; 
p7=p4; 
p8=0; 
p9=6.3e-4*2; 
p10=0.; 
p11= 8*%pi/(180); 
*/

//x=linspace(0,1,100); 
tt=0:%pi/200:%pi/2;
x=1-cos(tt);
upper(:,1) = x; 
lower(:,1) = x; 


f1=scf(1); clf; 

for isample = 1:nsamples; 

    p1 = P1(isample);
    p2 = P2(isample);
    p3 = P3(isample);
    p4 = P4(isample);
    p5 = P5(isample);
    p6 = P6(isample);
    p7 = P7(isample);
    p8 = P8(isample);
    p9 = P9(isample);
    p10 = P10(isample);
    p11 = P11(isample);
    
    Cup = [ones(1,6); 
           p2^(1/2) p2^(3/2) p2^(5/2) p2^(7/2) p2^(9/2) p2^(11/2);
           1/2      3/2      5/2      7/2      9/2      11/2;
           1/2*p2^(-1/2) 3/2*p2^(1/2) 5/2*p2^(3/2) 7/2*p2^(5/2) 9/2*p2^(7/2) 11/2*p2^(9/2);
           -1/4*p2^(-3/2) 3/4*p2^(-1/2) 15/4*p2^(1/2) 35/4*p2^(3/2) 63/4*p2^(5/2) 99/4*p2^(7/2);
           1 zeros(1,5)]; 

    bup = [p8 + p9/2; 
           p3;
           tan(p10 -p11/2);
           0;
           p4; 
           sqrt(2*p1)]; 

    aup = inv(Cup)*bup; 

    yup=aup(1)*x.^(1/2)+aup(2)*x.^(3/2)+aup(3)*x.^(5/2)+aup(4)*x.^(7/2)+aup(5)*x.^(9/2)+aup(6)*x.^(11/2);
    yup_naca = t*5*(0.2969*x.^(1/2)-0.1260*x-0.3516*x.^2+0.2843*x.^3-0.1015*x.^4);

    plot(x,yup,'b-',x,yup_naca,'r-');


end; // loop over samples 






