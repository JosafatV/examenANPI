#!/usr/bin/octave-cli

## Instituto Tecnológico de Costa Rica
## Área Académica de Ingeniería en Computadores
## CE-3102 Análisis Numérico para Ingeniería
## Prof. Pablo Alvarado
## II Semestre 2018
## Examen Final

## PROBLEMA 2

## NOMBRE: Josafat Vargas Gamboa
## CARNE:  2013030892

2;

## Construya algunos datos 2D para el problema
## N: Número de datos.  Cada fila de la matriz tendrá un punto.
##    La primera columna tendrá la coordenada x y la segunda columna
##    la coordenada y.
function points = createData(N)
  astep = 360/N;
  
  angles = (0:astep:360-astep)';
  anoise = rand(size(angles))*astep/4;
  rnoise = rand(size(angles))*1.5;
  
  radii  = 2+rnoise;
  angles = deg2rad(angles + anoise);
  
  points = [radii.*cos(angles) radii.*sin(angles)];
endfunction

## Calcule las segundas derivadas
## x: posiciones x de las muestras
## f: valores de la función en cada x
## retorne fpp con los valores de la segunda derivada en cada posición t
function fpp=findDerivs(x, f)
  #assert(length(f)==length(x));
  
  N=length(x); # Número de subintervalos
  
  ## Arme el sistema de ecuaciones
  M  =eye(N,N);
  fpp=zeros(N,1);
  b  =zeros(N,1);

  ## ################## 
  ## ## Problema 2.4 ##
  ## ################## 

  for i = 1:N
    
    # i plus one 
    ipo = i+1;
    
    # wrap-around
    if (ipo>N)
      ipo = 1;
    endif
    
    # i plus two
    ipt = ipo+1;
    
    # wrap-around
    if (ipt>N)
      ipt = 1;
    endif
    
    # insert values
    M(i,i)=x(ipo)-x(i);
    M(i,ipo)=2*(x(ipt)-x(i));
    M(i,ipt)=x(ipt)-x(ipo);

    b(i) = 6*( (f(ipt)-f(ipo))/(x(ipt)-x(ipo)) ) - 6*( (f(ipo)-f(i))/(x(ipo)-x(i)) );
    
  endfor

  ## Resuelva el sistema
  #fpp = M\b;
  fpp = linsolve( M, b ); # equivalent
  
endfunction

## Interpole los valores fi(xi) usando los puntos x y sus valores f(x)
## t: valores de soporte conocidos
## f: valores de la función en los x conocidos
## ts: valores en donde debe encontrarse la función interpolada
## retorna fs: valores de la función en los xs dados
function fs = interpole(t,f,ts)
  #assert(length(t)==length(f));

  ts=ts(:); ## Asegúrese de que es un vector columna
  fs=zeros(size(ts));
  N = length(ts)-1;
  
  ## Encuentre las segundas derivadas
  fpp=findDerivs(t,f);

  ## ##################
  ## ## Problema 2.5 ##
  ## ##################

  for h = 1:length(ts)
    
    # turns a value of ts into a valid index for t
    # allows the insertion of ts to get a vector with all the corresponding index
    i = lookup(t, ts(h), "l")
    l = i-1;
    
    # wrap-around
    if (l == 0)
      l = length(t);
    endif

    # common expression
    titl = t(i)-t(l);
    tlti = t(l)-t(i);

    # cubic parts
    tsi3 = (ts(h)-t(i)) * (ts(h)-t(i)) * (ts(h)-t(i));
    tsl3 = (ts(h)-t(l)) * (ts(h)-t(l)) * (ts(h)-t(l));
    
    # calculate in components
    A = (fpp(l) * tsi3) / (6*tlti);
    B = (fpp(i) * tsl3) / (6*titl);
    C = (f(l) / tlti) - (fpp(l) * (tlti) / 6);
    D = (f(i) / titl) - (fpp(i) * (titl) / 6);

    # calculate value
    fs(h) = (A + C*(ts(h)-t(i))) + (B + D*(ts(h)-t(l)));
  endfor
  
endfunction

## Depuración
figure(2,"name","Interpolación simple cerrada (depuración)");
x=[0,1,2,3];
f=[1,2,1,0.5];

hold off;
plot(x,f,'rx-;original;',"linewidth",2);

step=0.1;
xs=0:step:4-step;

fs=interpole(x,f,xs);

hold on;
plot(xs,fs,'bo-;interpolado;');
grid on;
xlabel("t");
ylabel("f(t)");

## El caso completo
N = 10;
D = createData(N);

## ##################
## ## Problema 2.6 ##
## ##################

D = [D; D(1,1), D(1,2)]; #Add first points at the end

figure(1,"name","Interpolación 2D cerrada");
hold off;
plot(D(:,1),D(:,2),'rx-',"linewidth",2);

step=0.1;
t=0:step:N-step;
xs=interpole([0:N-1]',D(:,1),t);
ys=interpole([0:N-1]',D(:,2),t);

hold on;
plot(xs,ys,'bo-');
xlabel("x");
ylabel("y");
grid;
