#!/usr/bin/octave-cli

## Instituto Tecnológico de Costa Rica
## Área Académica de Ingeniería en Computadores
## CE-3102 Análisis Numérico para Ingeniería
## Prof. Pablo Alvarado
## II Semestre 2018
## Examen Final

## PROBLEMA 1

## NOMBRE:
## CARNE:

1;

global m = 0.1;   ## Masa de la partícula
global b = 0.05;  ## Coeficiente de atenuación
global k = 1;     ## Constante de Hook

## ################## 
## ## Problema 1.1 ##
## ################## 
## Fuerza aplicada en la partícula
global F=@(x,v,t) 0;  ## <<< Ponga aquí su solución


## Resuelva el sistema atenuado masa resorte usando Euler
## tn el último instante de tiempo
## Dt paso temporal
function [t,x]=eulersys(tn,Dt)

  global m b k F;
  
  t=0:Dt:tn; ## Intervalo de simulación

  ## Pre-reserve la memoria utilizada.
  x=zeros(size(t));
  v=zeros(size(t));

  ## Condiciones iniciales
  x(1)=-1;

  ## ################## 
  ## ## Problema 1.2 ##
  ## ################## 

  ## Resuelva el sistema de ecuaciones con Euler

  ## >>> Ponga aquí su solución <<<
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

  ## Pre-reserve la memoria utilizada.
  x=zeros(size(t));
  v=zeros(size(t));

  ## Condiciones iniciales
  x(1)=-1;

  ## ################## 
  ## ## Problema 1.4 ##
  ## ################## 

  ## >>> Ponga aquí su solución <<<

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


