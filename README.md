# peco
clc
clear

//Ingreso de datos para el Problema
Eff= input('Ef (En GPa)=') ;
Ef=Eff*10^9 ;
Emm=input('Em (En GPa)=') ;
Em=Emm*10^9 ;
Gmm=input('Gm (En GPa)=') ;
Gm=Gmm*10^9 ;
Gff=input('Gf (En GPa)=') ;
Gf=Gff*10^9 ;
vf=input('Porcentaje Volumen Fibra (Entre 0 y 1)=') ;
vm=1-vf ;
vpf=input('Modulo Poisson Fibra=') ;
vpm=input('Modulo Poisson Matriz=') ;
cc=input('¿Ingrese la cantidad de capas?:    ') ;
z= input('¿Cual es el espesor de las capas? (Ingresar como matriz lineal) z=[]:    ') ;
L=length(z);
if L~=cc
   disp('Error, intentar de nuevo')
   return
end
ac=input('¿Que distribución tienen? (Indicar ordenadamente los ángulos a los que se encuentra cada capa) z=[]:    ') ;
L1=length(ac);
if L1~=cc
   disp('Error, intentar de nuevo')
   return
end

//Calculo de datos principales
v12= (vpf*vf) + (vpm*vm);
E1= (Ef*vf) + (Em*vm);
E2= 1/((vm/Em) + (vf/Ef));
v21=(E2/E1)*v12;
G12= 1/((vm/Gm) + (vf/Gf));
d=1-v12*v21;
Q11=E1/d;
Q12=(E2*v12)/d;
Q22=E2/d;
Q66=G12;
Q16=0;
Q26=0;
A11t=0;
A12t=0;
A16t=0;
A22t=0;
A26t=0;
A66t=0;
B11t=0;
B12t=0;
B16t=0;
B22t=0;
B26t=0;
B66t=0;
D11t=0;
D12t=0;
D16t=0;
D22t=0;
D26t=0;
D66t=0;
zi=sum(z)/2;

//Matriz a 0 y 90
      Q0=[Q11, Q12, Q16;
          Q12, Q22, Q26;
          Q16, Q26, Q66];
          
      Q90=[Q22, Q12, Q16;
          Q12, Q11, Q26;
          Q16, Q26, Q66];
      
//Matriz a 45 y -45
        x=45;
        Q11t= (Q11*cosd(x)^4)+((2*Q12+2*Q66)*(sind(x)^2)*(cosd(x)^2))+ (Q22*sind(x)^4);
        Q12t= (Q11+Q22+4*Q66)*(sind(x)^2)*(cosd(x)^2) + Q12*((sind(x)^4) + (cosd(x)^4));
        Q22t= (Q11*sind(x)^4)+((2*Q12+2*Q66)*(sind(x)^2)*(cosd(x)^2))+ (Q22*cosd(x)^4);
        Q16t= (Q11-Q12-2*Q66)*(sind(x))*(cosd(x)^3) + (Q12-Q22+2*Q66)*(sind(x)^3)*(cosd(x));
        Q26t= (Q11-Q12-2*Q66)*(sind(x)^3)*(cosd(x)) + (Q12-Q22+2*Q66)*(sind(x))*(cosd(x)^3);
        Q66t= (Q11+Q22-2*Q12-2*Q66)*(sind(x)^2)*(cosd(x)^2)+ Q66*((sind(x)^4)+(cosd(x)^4));
        
        Q45=[Q11t, Q12t, Q16t;
             Q12t, Q22t, Q26t;
             Q16t, Q26t, Q66t];
        
        Qm45=[Q11t, Q12t, -Q16t;
             Q12t, Q22t, -Q26t;
             -Q16t, -Q26t, Q66t];
         
//Matriz a 30
        x=30;
        Q11t= (Q11*cosd(x)^4)+((2*Q12+2*Q66)*(sind(x)^2)*(cosd(x)^2))+ (Q22*sind(x)^4);
        Q12t= (Q11+Q22+4*Q66)*(sind(x)^2)*(cosd(x)^2) + Q12*((sind(x)^4) + (cosd(x)^4));
        Q22t= (Q11*sind(x)^4)+((2*Q12+2*Q66)*(sind(x)^2)*(cosd(x)^2))+ (Q22*cosd(x)^4);
        Q16t= (Q11-Q12-2*Q66)*(sind(x))*(cosd(x)^3) + (Q12-Q22+2*Q66)*(sind(x)^3)*(cosd(x));
        Q26t= (Q11-Q12-2*Q66)*(sind(x)^3)*(cosd(x)) + (Q12-Q22+2*Q66)*(sind(x))*(cosd(x)^3);
        Q66t= (Q11+Q22-2*Q12-2*Q66)*(sind(x)^2)*(cosd(x)^2)+ Q66*((sind(x)^4)+(cosd(x)^4));
        
        Q30=[Q11t, Q12t, Q16t;
             Q12t, Q22t, Q26t;
             Q16t, Q26t, Q66t];
         
//Matriz a 60
        x=60;
        Q11t= (Q11*cosd(x)^4)+((2*Q12+2*Q66)*(sind(x)^2)*(cosd(x)^2))+ (Q22*sind(x)^4);
        Q12t= (Q11+Q22+4*Q66)*(sind(x)^2)*(cosd(x)^2) + Q12*((sind(x)^4) + (cosd(x)^4));
        Q22t= (Q11*sind(x)^4)+((2*Q12+2*Q66)*(sind(x)^2)*(cosd(x)^2))+ (Q22*cosd(x)^4);
        Q16t= (Q11-Q12-2*Q66)*(sind(x))*(cosd(x)^3) + (Q12-Q22+2*Q66)*(sind(x)^3)*(cosd(x));
        Q26t= (Q11-Q12-2*Q66)*(sind(x)^3)*(cosd(x)) + (Q12-Q22+2*Q66)*(sind(x))*(cosd(x)^3);
        Q66t= (Q11+Q22-2*Q12-2*Q66)*(sind(x)^2)*(cosd(x)^2)+ Q66*((sind(x)^4)+(cosd(x)^4));
        
        Q60=[Q11t, Q12t, Q16t;
             Q12t, Q22t, Q26t;
             Q16t, Q26t, Q66t];
         
QT= [Q0 Q30 Q45 Q60 Q90 Qm45];


//Calculo valores de matriz de fuerzas y momentos
for i = 1:cc
    if ac(i) == 0
        
        A11=QT(1,1)*(z(i));
        A12=QT(1,2)*(z(i));
        A16=QT(1,3)*(z(i));
        A22=QT(2,2)*(z(i));
        A26=QT(2,3)*(z(i));
        A66=QT(3,3)*(z(i));
        
       
        B11=QT(1,1)*(z(i))*(zi-z(i)/2);
        B12=QT(1,2)*(z(i))*(zi-z(i)/2);
        B16=QT(1,3)*(z(i))*(zi-z(i)/2);
        B22=QT(2,2)*(z(i))*(zi-z(i)/2);
        B26=QT(2,3)*(z(i))*(zi-z(i)/2);
        B66=QT(3,3)*(z(i))*(zi-z(i)/2);
        
        
        D11=QT(1,1)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D12=QT(1,2)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D16=QT(1,3)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D22=QT(2,2)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D26=QT(2,3)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D66=QT(3,3)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        
    elseif ac(i) == 30
        A11=QT(1,4)*(z(i));
        A12=QT(1,5)*(z(i));
        A16=QT(1,6)*(z(i));
        A22=QT(2,5)*(z(i));
        A26=QT(2,6)*(z(i));
        A66=QT(3,6)*(z(i));
        
        B11=QT(1,4)*(z(i))*(zi-z(i)/2);
        B12=QT(1,5)*(z(i))*(zi-z(i)/2);
        B16=QT(1,6)*(z(i))*(zi-z(i)/2);
        B22=QT(2,5)*(z(i))*(zi-z(i)/2);
        B26=QT(2,6)*(z(i))*(zi-z(i)/2);
        B66=QT(3,6)*(z(i))*(zi-z(i)/2);
        
        D11=QT(1,4)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D12=QT(1,5)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D16=QT(1,6)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D22=QT(2,5)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D26=QT(2,6)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D66=QT(3,6)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        
    elseif ac(i) == 45
        A11=QT(1,7)*(z(i));
        A12=QT(1,8)*(z(i));
        A16=QT(1,9)*(z(i));
        A22=QT(2,8)*(z(i));
        A26=QT(2,9)*(z(i));
        A66=QT(3,9)*(z(i));
        
        B11=QT(1,7)*(z(i))*(zi-z(i)/2);
        B12=QT(1,8)*(z(i))*(zi-z(i)/2);
        B16=QT(1,9)*(z(i))*(zi-z(i)/2);
        B22=QT(2,8)*(z(i))*(zi-z(i)/2);
        B26=QT(2,9)*(z(i))*(zi-z(i)/2);
        B66=QT(3,9)*(z(i))*(zi-z(i)/2);
        
        D11=QT(1,7)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D12=QT(1,8)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D16=QT(1,9)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D22=QT(2,8)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D26=QT(2,9)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D66=QT(3,9)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        
    elseif ac(i) == 60
        A11=QT(1,10)*(z(i));
        A12=QT(1,11)*(z(i));
        A16=QT(1,12)*(z(i));
        A22=QT(2,11)*(z(i));
        A26=QT(2,12)*(z(i));
        A66=QT(3,12)*(z(i));
        
        B11=QT(1,10)*(z(i))*(zi-z(i)/2);
        B12=QT(1,11)*(z(i))*(zi-z(i)/2);
        B16=QT(1,12)*(z(i))*(zi-z(i)/2);
        B22=QT(2,11)*(z(i))*(zi-z(i)/2);
        B26=QT(2,12)*(z(i))*(zi-z(i)/2);
        B66=QT(3,12)*(z(i))*(zi-z(i)/2);
        
        D11=QT(1,10)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D12=QT(1,11)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D16=QT(1,12)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D22=QT(2,11)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D26=QT(2,12)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D66=QT(3,12)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        
            elseif ac(i) == 90
        A11=QT(1,13)*(z(i));
        A12=QT(1,14)*(z(i));
        A16=QT(1,15)*(z(i));
        A22=QT(2,14)*(z(i));
        A26=QT(2,15)*(z(i));
        A66=QT(3,15)*(z(i));
        
        B11=QT(1,13)*(z(i))*(zi-z(i)/2);
        B12=QT(1,14)*(z(i))*(zi-z(i)/2);
        B16=QT(1,15)*(z(i))*(zi-z(i)/2);
        B22=QT(2,14)*(z(i))*(zi-z(i)/2);
        B26=QT(2,15)*(z(i))*(zi-z(i)/2);
        B66=QT(3,15)*(z(i))*(zi-z(i)/2);
        
        D11=QT(1,13)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D12=QT(1,14)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D16=QT(1,15)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D22=QT(2,14)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D26=QT(2,15)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D66=QT(3,15)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));

            elseif ac(i) == -45
        A11=QT(1,16)*(z(i));
        A12=QT(1,17)*(z(i));
        A16=QT(1,18)*(z(i));
        A22=QT(2,17)*(z(i));
        A26=QT(2,18)*(z(i));
        A66=QT(3,18)*(z(i));
        
        B11=QT(1,16)*(z(i))*(zi-z(i)/2);
        B12=QT(1,17)*(z(i))*(zi-z(i)/2);
        B16=QT(1,18)*(z(i))*(zi-z(i)/2);
        B22=QT(2,17)*(z(i))*(zi-z(i)/2);
        B26=QT(2,18)*(z(i))*(zi-z(i)/2);
        B66=QT(3,18)*(z(i))*(zi-z(i)/2);
        
        D11=QT(1,16)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D12=QT(1,17)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D16=QT(1,18)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D22=QT(2,17)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D26=QT(2,18)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
        D66=QT(3,18)*((z(i)*((zi-z(i)/2)^2))+((z(i)^3)/12));
    end
    
    zi=zi-z(i);
    
    A11t=A11t+A11;
    A12t=A12t+A12;
    A16t=A16t+A16;
    A22t=A22t+A22;
    A26t=A26t+A26;
    A66t=A66t+A66;
    
    B11t=B11t+B11;
    B12t=B12t+B12;
    B16t=B16t+B16;
    B22t=B22t+B22;
    B26t=B26t+B26;
    B66t=B66t+B66;
    
    D11t=D11t+D11;
    D12t=D12t+D12;
    D16t=D16t+D16;
    D22t=D22t+D22;
    D26t=D26t+D26;
    D66t=D66t+D66;
      
end

//Llenado de datos en matriz
Mdef= [ A11t A12t A16t B11t B12t B16t;
        A12t A22t A26t B12t B22t B26t;
        A16t A26t A66t B16t B26t B66t;
        B11t B12t B16t D11t D12t D16t;
        B12t B22t B26t D12t D22t D26t;
        B16t B26t B66t D16t D26t D66t];
    
Mfin=inv(Mdef);
        
Nxx=input('Indique las fuerzas que actuan sobre el sistema [Nx,Ny,Nxy]:    ');
Nx=Nxx';
Mxx=input('Indique los momentos que actuan sobre el sistema [Mx,My,Mxy]:    ');
Mx=Mxx';

Fu=[Nx;Mx];

Res=Mfin*Fu;

Disp=(Res)
