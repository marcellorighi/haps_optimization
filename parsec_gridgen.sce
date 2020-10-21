clear; 

nsamples = 24; 
nactive=4; 

RV=grand(nactive,"prm",(0:(nsamples-1))')/(nsamples) + grand(nsamples,nactive,"unf",0,1/nsamples)-0.5;


t=0.06; 
eps = 0.075; 

dt  = eps * RV(:,1) * t; 
P1 = 1.1019*(t + dt).^2;
P2 = 0.30 * (1 + eps * RV(:,2)); 
P3 = 0.5 * (t + dt); 
P4 = -0.2125 * (1 + RV(:,3) * eps); 
P5 = P1; 
P6 = P2; 
P7 = P4; 
P8 = zeros(P1); 
P9 = 6.3e-4*2 * ones(P1); 
P10= zeros(P1); 
P11= 8*%pi/(180) * (1 + RV(:,4) * eps); 



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
xx=1-cos(tt);
upper(:,1) = xx; 
lower(:,1) = xx; 

exec('/home/rigm/SU2/runs/Scilab/airfoil/param/ogrid_func.sce');

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

    yup=aup(1)*xx.^(1/2)+aup(2)*xx.^(3/2)+aup(3)*xx.^(5/2)+aup(4)*xx.^(7/2)+aup(5)*xx.^(9/2)+aup(6)*xx.^(11/2);
    yup_naca = t*5*(0.2969*xx.^(1/2)-0.1260*xx-0.3516*xx.^2+0.2843*xx.^3-0.1015*xx.^4);

    plot(xx,yup,'b-',xx,yup_naca,'r-');

    upper(:,2)=yup;
    lower(:,2)=-yup;


    Dir='/home/rigm/SU2/runs/If/NACA0006/UQ/M015/pert_'+string(isample);

    unix('mkdir -p ' + Dir); // /home/staff/rigm/unix/SU2/runs/awp-testing/airfoil10/pert1/pert'+string(isample)); 

    fileNameG=Dir+'/mesh.su2'; 


    nPts=256;
    nPts2=256;
    nJ=128;
    hBL=5.0e-6;


    Radius=500.;

    [x,y,nX]=ogrid_func(upper,lower,nPts,nPts2,nJ,1.125,hBL,Radius);



f4=scf(4);
clf();
a=get("current_axes");
a.rotation_angles=[90,0];
title('$\Large \text{Airfoil}$');
mesh(x,y,zeros(x)); 
xgrid;

// numbering points 

for j=1:nJ; 
    for i=1:nX; 
       npoin(i,j)=i-1+(nX)*(j-1); 
    end; 
end;

// numbering elements

for j=1:nJ-1; 
    for i=1:nX; 
       nelem(i,j)=i-1+(nX-1)*(j-1); 
    end; 
end;

//

// saving Design Param perturbations 
fileNameParam=Dir+'/para.csv';

fd=mopen(fileNameParam,"wt");
mfprintf(fd,'DP   Pert\n');
mfprintf(fd,'%3.0f,%13.8f\n',[1:4]',RV(isample,:)');
mclose(fd);


fdg=mopen(fileNameG,'wt');
mfprintf(fdg,'NDIME=2\n');
mfprintf(fdg,'NELEM=%i\n',nX*(nJ-1));

for j=1:nJ-1
    for i=1:nX-1
        mfprintf(fdg,'%i %i %i %i %i %i\n',9,npoin(i,j),npoin(i,j+1),npoin(i+1,j+1),npoin(i+1,j),nelem(i,j)); 
    end;    
    mfprintf(fdg,'%i %i %i %i %i %i\n',9,npoin(nX,j),npoin(nX,j+1),npoin(1,j+1),npoin(1,j),nelem(nX,j)); 
end;

mfprintf(fdg,'NPOIN=%i\n',(nX)*(nJ));
for j=1:nJ
    for i=1:nX
        mfprintf(fdg,'%e %e %i\n',x(i,j),y(i,j),npoin(i,j)); 
    end;    
end;

nmark=2; 
mfprintf(fdg,'NMARK=%d\n',nmark);

mfprintf(fdg,'MARKER_TAG=%s\n','Airfoil'); 
mfprintf(fdg,'MARKER_ELEMS=%i\n',nX); 

for i=2:nX; 
    mfprintf(fdg,'%i %i %i\n',3,npoin(i,1),npoin(i-1,1));
end;
mfprintf(fdg,'%i %i %i\n',3,npoin(1,1),npoin(nX,1));

mfprintf(fdg,'MARKER_TAG=%s\n','Farfield'); 
mfprintf(fdg,'MARKER_ELEMS=%i\n',nX);


for i=1:nX-1; 
    mfprintf(fdg,'%i %i %i\n',3,npoin(i,nJ),npoin(i+1,nJ));
end;
mfprintf(fdg,'%i %i %i\n',3,npoin(nX,nJ),npoin(1,nJ));
 
mclose(fdg);
//

end; // loop over samples 






