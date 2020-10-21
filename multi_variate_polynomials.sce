function y=h1(x)
    y = x; 
endfunction

function y=h2(x)

    y = -1 + x^2; 
//    y = 0.5*(3*x^2 -1); 

endfunction

function y=h3(x)

    y = -3*x + x^3; 
//    y = 0.5*(5*x^3 - 3*x); 

endfunction

function y=h4(x)

    y = 3 - 6*x^2 + x^4; 
//    y = 0.125*(35*x^4 - 30*x^2 + 3); 

endfunction


function y=h5(x)

    y = 15*x - 10*x^3 + x^5; 
//    y = 0.125*(63*x^5 - 70*x^3 + 15*x ); 

endfunction



function Hijkl=Psi(x1,x2,x3,x4)
    
// 0th order 
    Hijkl(1)=1.; 

// 1st order
    Hijkl(2)=h1(x1);
    Hijkl(3)=h1(x2);
    Hijkl(4)=h1(x3);
    Hijkl(5)=h1(x4);
    
// 2nd order 
    ii=5; 
    Hijkl(ii+1)=h2(x1);
    Hijkl(ii+2)=h2(x2);
    Hijkl(ii+3)=h2(x3);
    Hijkl(ii+4)=h2(x4);
    
    Hijkl(ii+5) =h1(x1)*h1(x2);
    Hijkl(ii+6) =h1(x1)*h1(x3);
    Hijkl(ii+7) =h1(x1)*h1(x4);
    Hijkl(ii+8) =h1(x2)*h1(x3);
    Hijkl(ii+9) =h1(x2)*h1(x4);
    Hijkl(ii+10)=h1(x3)*h1(x4);
    
// 3rd order 
    ii=15; 
    Hijkl(ii+1)=h3(x1);
    Hijkl(ii+2)=h3(x2);
    Hijkl(ii+3)=h3(x3);
    Hijkl(ii+4)=h3(x4);

    Hijkl(ii+5) =h1(x1)*h1(x2)*h1(x3);
    Hijkl(ii+6) =h1(x1)*h1(x3)*h1(x4);
    Hijkl(ii+7) =h1(x1)*h1(x2)*h1(x4);
    Hijkl(ii+8) =h1(x2)*h1(x3)*h1(x4);

    Hijkl(ii+9) =h2(x1)*h1(x2);
    Hijkl(ii+10)=h2(x1)*h1(x3);
    Hijkl(ii+11)=h2(x1)*h1(x4);
    Hijkl(ii+12)=h2(x2)*h1(x1);
    Hijkl(ii+13)=h2(x2)*h1(x3);
    Hijkl(ii+14)=h2(x2)*h1(x4);
    Hijkl(ii+15)=h2(x3)*h1(x1);
    Hijkl(ii+16)=h2(x3)*h1(x2);
    Hijkl(ii+17)=h2(x3)*h1(x4);
    Hijkl(ii+18)=h2(x4)*h1(x1);
    Hijkl(ii+19)=h2(x4)*h1(x2);
    Hijkl(ii+20)=h2(x4)*h1(x3);

// 4th order 
    ii=35; 
    Hijkl(ii+1)=h4(x1);
    Hijkl(ii+2)=h4(x2);
    Hijkl(ii+3)=h4(x3);
    Hijkl(ii+4)=h4(x4);

    Hijkl(ii+5) =h3(x1)*h1(x2);
    Hijkl(ii+6) =h3(x1)*h1(x3);
    Hijkl(ii+7) =h3(x1)*h1(x4);
    Hijkl(ii+8) =h3(x2)*h1(x1);
    Hijkl(ii+9) =h3(x2)*h1(x3);
    Hijkl(ii+10)=h3(x2)*h1(x4);
    Hijkl(ii+11)=h3(x3)*h1(x1);
    Hijkl(ii+12)=h3(x3)*h1(x2);
    Hijkl(ii+13)=h3(x3)*h1(x4);
    Hijkl(ii+14)=h3(x4)*h1(x1);
    Hijkl(ii+15)=h3(x4)*h1(x2);
    Hijkl(ii+16)=h3(x4)*h1(x3);
    
    Hijkl(ii+17)=h2(x1)*h2(x2);
    Hijkl(ii+18)=h2(x1)*h2(x3);
    Hijkl(ii+19)=h2(x1)*h2(x4);
    Hijkl(ii+20)=h2(x2)*h2(x3);
    Hijkl(ii+21)=h2(x2)*h2(x4);
    Hijkl(ii+22)=h2(x3)*h2(x4);
    
    Hijkl(ii+23)=h2(x1)*h1(x2)*h1(x3);
    Hijkl(ii+24)=h2(x1)*h1(x2)*h1(x4);
    Hijkl(ii+25)=h2(x1)*h1(x3)*h1(x4);
    Hijkl(ii+26)=h2(x2)*h1(x1)*h1(x3);
    Hijkl(ii+27)=h2(x2)*h1(x1)*h1(x4);
    Hijkl(ii+28)=h2(x2)*h1(x3)*h1(x4);
    Hijkl(ii+29)=h2(x3)*h1(x1)*h1(x2);
    Hijkl(ii+30)=h2(x3)*h1(x1)*h1(x4);
    Hijkl(ii+31)=h2(x3)*h1(x2)*h1(x4);
    Hijkl(ii+32)=h2(x4)*h1(x1)*h1(x2);
    Hijkl(ii+33)=h2(x4)*h1(x1)*h1(x3);
    Hijkl(ii+34)=h2(x4)*h1(x2)*h1(x3);

    Hijkl(ii+35)=h1(x1)*h1(x2)*h1(x3)*h1(x4);

// 5th order 
    ii=70; 
    Hijkl(ii+1)=h5(x1);
    Hijkl(ii+2)=h5(x2);
    Hijkl(ii+3)=h5(x3);
    Hijkl(ii+4)=h5(x4);

    Hijkl(ii+5) =h4(x1)*h1(x2);
    Hijkl(ii+6) =h4(x1)*h1(x3);
    Hijkl(ii+7) =h4(x1)*h1(x4);
    Hijkl(ii+8) =h4(x2)*h1(x1);
    Hijkl(ii+9) =h4(x2)*h1(x3);
    Hijkl(ii+10)=h4(x2)*h1(x4);
    Hijkl(ii+11)=h4(x3)*h1(x1);
    Hijkl(ii+12)=h4(x3)*h1(x2);
    Hijkl(ii+13)=h4(x3)*h1(x4);
    Hijkl(ii+14)=h4(x4)*h1(x1);
    Hijkl(ii+15)=h4(x4)*h1(x2);
    Hijkl(ii+16)=h4(x4)*h1(x3);
        
    Hijkl(ii+17)=h3(x1)*h2(x2);
    Hijkl(ii+18)=h3(x1)*h2(x3);
    Hijkl(ii+19)=h3(x1)*h2(x4);
    Hijkl(ii+20)=h3(x2)*h2(x1);
    Hijkl(ii+21)=h3(x2)*h2(x3);
    Hijkl(ii+22)=h3(x2)*h2(x4);
    Hijkl(ii+23)=h3(x3)*h2(x1);
    Hijkl(ii+24)=h3(x3)*h2(x2);
    Hijkl(ii+25)=h3(x3)*h2(x4);
    Hijkl(ii+26)=h3(x4)*h2(x1);
    Hijkl(ii+27)=h3(x4)*h2(x2);
    Hijkl(ii+28)=h3(x4)*h2(x3);

    Hijkl(ii+29)=h3(x1)*h1(x2)*h1(x3);
    Hijkl(ii+30)=h3(x1)*h1(x2)*h1(x4);
    Hijkl(ii+31)=h3(x1)*h1(x3)*h1(x4);
    Hijkl(ii+32)=h3(x2)*h1(x1)*h1(x3);
    Hijkl(ii+33)=h3(x2)*h1(x1)*h1(x4);
    Hijkl(ii+34)=h3(x2)*h1(x3)*h1(x4);
    Hijkl(ii+35)=h3(x3)*h1(x1)*h1(x2);
    Hijkl(ii+36)=h3(x3)*h1(x1)*h1(x4);
    Hijkl(ii+37)=h3(x3)*h1(x2)*h1(x4);
    Hijkl(ii+38)=h3(x4)*h1(x1)*h1(x2);
    Hijkl(ii+39)=h3(x4)*h1(x1)*h1(x3);
    Hijkl(ii+40)=h3(x4)*h1(x2)*h1(x3);

    Hijkl(ii+41)=h2(x1)*h1(x2)*h1(x3)*h1(x4);
    Hijkl(ii+42)=h2(x2)*h1(x1)*h1(x3)*h1(x4);
    Hijkl(ii+43)=h2(x3)*h1(x1)*h1(x2)*h1(x4);
    Hijkl(ii+44)=h2(x4)*h1(x1)*h1(x2)*h1(x3);

    Hijkl(ii+45)=h2(x1)*h2(x2)*h1(x3);
    Hijkl(ii+46)=h2(x1)*h2(x2)*h1(x4);
    Hijkl(ii+47)=h2(x1)*h2(x3)*h1(x2);
    Hijkl(ii+48)=h2(x1)*h2(x3)*h1(x4);
    Hijkl(ii+49)=h2(x1)*h2(x4)*h1(x2);
    Hijkl(ii+50)=h2(x1)*h2(x4)*h1(x3);
    Hijkl(ii+51)=h2(x2)*h2(x3)*h1(x1);
    Hijkl(ii+52)=h2(x2)*h2(x3)*h1(x4);
    Hijkl(ii+53)=h2(x2)*h2(x4)*h1(x1);
    Hijkl(ii+54)=h2(x2)*h2(x4)*h1(x3);
    Hijkl(ii+55)=h2(x3)*h2(x4)*h1(x1);
    Hijkl(ii+56)=h2(x3)*h2(x4)*h1(x2);

endfunction


function S=sobol_indices(next,nint)
    
xitest1=grand(next,1,"nor",0.0,1.00); 

for j=1:next; 
    xitest24=grand(nint,3,"nor",0.0,1.00);
    for k=1:nint; 
        alpha = Psi(xitest1(j),xitest24(k,1),xitest24(k,2),xitest24(k,3)); 
        cltestSI(k)=aa3'*alpha(1:35); //aa3'*alpha;
    end;     
    E(j) = sum(cltestSI)/nint; 
end

EE1 = sum(E)/next; 
V1 = sum((E -EE1).^2)/next; 
S(1) = V1/Vtest; disp(S(1));

xitest1=grand(next,1,"nor",0.0,1.00); 

for j=1:next; 
    xitest24=grand(nint,3,"nor",0.0,1.00);
    for k=1:nint; 
        alpha = Psi(xitest24(k,1),xitest1(j),xitest24(k,2),xitest24(k,3)); 
        cltestSI(k)=aa3'*alpha(1:35); //aa3'*alpha;
    end;     
    E(j) = sum(cltestSI)/nint; 
end

EE1 = sum(E)/next; 
V1 = sum((E -EE1).^2)/next; 
S(2) = V1/Vtest; disp(S(2));

xitest1=grand(next,1,"nor",0.0,1.00); 

for j=1:next; 
    xitest24=grand(nint,3,"nor",0.0,1.00);
    for k=1:nint; 
        alpha = Psi(xitest24(k,1),xitest24(k,2),xitest1(j),xitest24(k,3)); 
        cltestSI(k)=aa3'*alpha(1:35); //aa3'*alpha;
    end;     
    E(j) = sum(cltestSI)/nint; 
end

EE1 = sum(E)/next; 
V1 = sum((E -EE1).^2)/next; 
S(3) = V1/Vtest; disp(S(3));

xitest1=grand(next,1,"nor",0.0,1.00); 

for j=1:next; 
    xitest24=grand(nint,3,"nor",0.0,1.00);
    for k=1:nint; 
        alpha = Psi(xitest24(k,1),xitest24(k,2),xitest24(k,3),xitest1(j)); 
        cltestSI(k)=aa3'*alpha(1:35); //aa3'*alpha;
    end;     
    E(j) = sum(cltestSI)/nint; 
end

EE1 = sum(E)/next; 
V1 = sum((E -EE1).^2)/next; 
S(4) = V1/Vtest; disp(S(4));

endfunction
