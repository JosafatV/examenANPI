#!/usr/bin/octave-cli

## Instituto Tecnológico de Costa Rica
## Área Académica de Ingeniería en Computadores
## CE-3102 Análisis Numérico para Ingeniería
## Prof. Pablo Alvarado
## II Semestre 2018
## Examen Final

## PROBLEMA 1

## NOMBRE: Josafat Vargas Gamboa
## CARNE:  2013030892

1;

global m = 0.1;   ## Masa de la partícula
global b = 0.05;  ## Coeficiente de atenuación
global k = 1;     ## Constante de Hook

## ################## 
## ## Problema 1.1 ##
## ################## 
## Fuerza aplicada en la partícula
global F = @(x,v,t) -b*v-k*x;

## Resuelva el sistema atenuado masa resorte usando Euler
## tn el último instante de tiempo
## Dt tamaño de paso temporal
function [t,x] = eulersys(tn,Dt)

  global m b k F;
  
  t=0:Dt:tn; ## Intervalo de simulación
  n=length(t)-1;

  ## Pre-reserve la memoria utilizada.
  x=zeros(size(t));
  v=zeros(size(t));

  ## Condiciones iniciales
  x(1)=-1;
  v(1)= 0;

  ## ################## 
  ## ## Problema 1.2 ##
  ## ################## 

  ## Resuelve el sistema de ecuaciones con Euler
  for i = 1:n
    ti = t(i);
    xi = x(i);
    vi = v(i);
    
    #x(i+1) = xi+F(xi, vi, ti)*Dt;
    x(i+1) = xi+F(xi, vi, ti)*Dt+(vi*m)*Dt*Dt*0.5; #Euler para orden superior
    v(i+1) = vi+F(xi, vi, ti)*Dt;
    
   endfor
   
endfunction
  
figure(1,"name","Euler");
hold off;
[t,x]=eulersys(10,0.05);
plot(t,x,"r;\\Delta t=0.05;");

hold on;
[t,x]=eulersys(10,0.01);
plot(t,x,"m;\\Delta t=0.01;");

[t,x]=eulersys(10,0.001);
plot(t,x,"b;\\Delta t=0.001;");

xlabel("t");
ylabel("x(t)");
axis([0,10,-2,2]);
grid on;


## Resuelva el sistema de ecuaciones con Runge-Kutta 4to orden
## tn Último instante de tiempo
## Dt Paso temporal (delta t)
function [t,x] = rksys(tn,Dt)
  global m b k F;

  t=0:Dt:tn; ## Intervalo de simulación
  n=length(t)-1;

  ## Pre-reserve la memoria utilizada.
  x=zeros(size(t));
  v=zeros(size(t));

  ## Condiciones iniciales
  x(1)=-1;
  v(1)= 0;

  ## ################## 
  ## ## Problema 1.4 ##
  ## ################## 

  for i = 1:n
      ti = t(i); #i=1: 0
      vi = v(i); #i=1: 0
      xi = x(i); #i=1: -1
      
          #calculate l and k values
    k1 = F(xi, vi, ti)*Dt;
        l1 = F(xi, vi, ti)*Dt;
    k2 = F(xi+k1*0.5, vi+l1*0.5, ti+0.5*Dt)*Dt;
        l2 = F(xi+k1*0.5, vi+l1*0.5, ti+0.5*Dt)*Dt;
    k3 = F(xi+k2*0.5, vi+l2*0.5, ti+0.5*Dt)*Dt;
        l3 = F(xi+k2*0.5, vi+l1*0.5, ti+0.5*Dt)*Dt;
    k4 = F(xi+k3, vi+l3, ti+Dt)*Dt;
        l4 = F(xi+k3, vi+l3, ti+Dt)*Dt;

    #calculate slope
    Xphi = (1/6)*(k1 + 2*k2 + 2*k3 + k4);
    Vphi = (1/6)*(l1 + 2*l2 + 2*l3 + l4);
    
    # store and step
    x(i+1) = xi+Xphi;
    v(i+1) = vi+Vphi;
    
  endfor

endfunction

figure(2,"name","RK");
hold off;
[t,x]=rksys(10,0.05);
plot(t,x,"r;\\Delta t=0.05;");

hold on;
[t,x]=rksys(10,0.01);
plot(t,x,"m;\\Delta t=0.01;");

[t,x]=rksys(10,0.001);
plot(t,x,"b;\\Delta t=0.001;");

xlabel("t");
ylabel("x(t)");
axis([0,10,-2,2]);
grid on;
