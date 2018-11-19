#!/usr/bin/octave-cli

## Instituto Tecnológico de Costa Rica
## Área Académica de Ingeniería en Computadores
## CE-3102 Análisis Numérico para Ingeniería
## Prof. Pablo Alvarado
## II Semestre 2018
## Examen Final

## PROBLEMA 3

## NOMBRE: Josafat Vargas Gamboa
## CARNE:  2013030892

3;

## Cargue los datos de distancias
##
## Cada columna de D tiene las distancias a un sensor fijo en el
## espacio. 
D=load("-ascii","dists.dat");

## Las posiciones de los emisores están dadas por las _columnas_
## de la matriz S
E=load("-ascii","emisors.dat");

## Función calcula posición a partir de distancias a emisores.
##
## dists:     Distancias del punto a cada emisor como vector fila.
## emisorPos: Posiciones de los emisores, cada emisor en una columna.
## option:    Si solo hay 3 emisores, cuál de las soluciones retorna:
##            1: use + en la cuadrática
##            2: use - en la cuadrática
##            3: decida cuál signo usar dependiendo de la velocidad y
##               aceleración
##
## La posición i de dists contiene la distancia al emisor en la
## columna i.
##
function p=calcPosition(dists,emisorPos,option=1)

  ## Número de emisores usados
  dim = min(columns(emisorPos),length(dists));

  ## ##################
  ## ## Problema 3.1 ##
  ## ##################
  
  M = zeros(dim,4); 
  
    ## Fill in the positions
  for i = 1:dim
    ei = [1; emisorPos(:,i)]; # make a 4x1 column with the position of one emitter
    M(i,:) = M(i,:) + transpose(ei); # fill the row by adding the position
  endfor

  ## ##################
  ## ## Problema 3.2 ##
  ## ##################
    
  # obtain magnitud of the position vector of emitters 
  b = zeros(dim,1);
  
  # obtain column with the square of distances
  dis2 = transpose (dists);
  dis2 = dis2 .* dis2;
  b = b .+ dis2;
  
  for i = 1:dim
      ei = norm( emisorPos(:,i) );
      ei = ei * ei;
      b(i) = b(i) - ei;
  endfor
  
  iM=pinv(M);

  ## ##################
  ## ## Problema 3.3 ##
  ## ##################

  ## Calcule la matriz seudo-inversa utilizando SVD
  [v, w ,u] = svd(M, "econ");
     
  # Get the inverse values of the diagonal. Rewrites w
  for i=1:length(w);
    if (w(i, i) == 0)
      wi = 0;
    else
      wi = 1/w(i, i);
    endif
    w(i,i) = wi;
  endfor
  
  #calculate inverse
  piM = transpose(v * w * transpose(u));
  #M = piM; # rewrite M as to not modify the code given

  ## Verifique que iM y pinv(M) son lo mismo
  if (norm(iM-piM,"fro") > 1e-6)
    error("Matriz inversa calculada con SVD incorrecta");
  endif

  ## ##################
  ## ## Problema 3.4 ##
  ## ##################
  
  ## Calcule la solución particular
  hatp = zeros(4,1);
  hatp = piM * b;  
    
  ## El caso de 3 dimensiones tiene dos posibles soluciones:
  if (dim==3)
    ## Con 3 emisores, calcule las dos posibles posiciones
    
    ## ##################
    ## ## Problema 3.5 ##
    ## ##################
    p=zeros(3,1);
    n=ones(4,1);
    n=null(M);
    
    a = n(2)*n(2) + n(3)*n(3) + n(4)*n(4);
    b = 2 * ( hatp(2)*n(2) + hatp(3)*n(3) + hatp(4)*n(4) - 0.5*hatp(1) );
    c = hatp(2)*hatp(2) + hatp(3)*hatp(3) + hatp(4)*hatp(4) - hatp(1);
    
    #Calcula una función cuadrática con la fórmula alternativa y precisión doble
  	disc = sqrt((b*b)-(4*a*c)); # calculo del discriminante 
    
	  if (option)
	  	lambda = -b+disc/(2*a);
	  else
	  	lambda = -b-disc/(2*a);
	  endif
    
    # lambda must be real
    if (imag(lambda) != 0)
      #According to Norrdine if the system is not solvable aproximate with the real part:
      real(lambda);
    endif
      
    p = hatp + lambda * n;
    p = p(2: 4); # remove distance value
    
  else 
    ## Caso general de más de tres dimensiones

    ## ##################
    ## ## Problema 3.6 ##
    ## ##################
    
    p = hatp(2: 4); # remove distance value
    
  endif
  
endfunction

## Calcule la trayectoria de posiciones para una matriz de
## distancias que contiene en cada fila el vector de distancias
## a cada emisor.
## 
## Sea N el número de datos y S el número de emisores
## dists:     matriz de distancias de tamaño N x S
## emisorPos: matrix de tamaño 3xS
## option:    si S==3, entonces cuál de las dos soluciones devolver
##            option=1 implica usar + en cuadrática, y otra cosa usa -
function p=calcPositions(dists,emisorPos,option=1)
  n=min(columns(emisorPos),columns(dists));
  k=rows(dists);

  ## Reserve memoria para todos los puntos en la trayectoria
  p=zeros(k,3);
  dp=zeros(k,3);
  d2p=zeros(k,3);

  if (option==3)
    
    predPos=[2;0;1]; ## Predicción inicial de posición
    
    ## Para cada punto en la secuencia
    for i=1:k

      ## Calcule las dos posibilidades
      sol1=calcPosition(dists(i,1:n),emisorPos(:,1:n),1)';
      sol2=calcPosition(dists(i,1:n),emisorPos(:,1:n),2)';
      
      ## Verifique cuál solución está más cerca de la predicción
      if (norm(predPos-sol1) < norm(predPos-sol2))
	      p(i,:) = sol1;
      else
	      p(i,:) = sol2;
      endif

      ## ##################
      ## ## Problema 3.7 ##
      ## ##################

      # asumo una distancia de paso constante y un lambda subrelajado
      dt = 0.1;
      lamda = 0.8;
      p(1,:) = transpose(predPos);
      
      #aproximación de la primer y segunda derivada
      if ( i>=2 )
        dp(i,:) = (p(i,:)-p(i-1,:)) / dt;
        dp(i,:) = dp(i,:)*lamda - lamda * dp(i-1,:);
      endif
      
      if ( i>=3 )
          d2p(i,:) = (p(i,:)-2*p(i-1,:)+p(i-2,:)) / dt*dt;
          d2p(i,:) = lamda.*d2p(i,:) - (1-lamda).*d2p(i-1,:);
      endif

        ## Actualice la predicción (CRASHES OCTAVE -> data values greater than float capacity)
      #p(i,:) = p(i,:) + dp(i,:)*dt + 0.5*(d2p(i,:)*dt*dt);
      
  endfor
  
  else
    ## Para todos los puntos en la trayectoria
    for i=1:k
      ## Calcule la posición del punto, dadas las distancias a los emisores
      p(i,:)=calcPosition(dists(i,1:n),emisorPos(:,1:n),option)';
    endfor
  endif

endfunction

## Función calcula distancia de un punto a todos los emisores
## pos:       Posición del objeto (en las filas)
## emisorPos: Posiciones de los emisores, en las columnas.
function d=calcDistances(pos,emisorPos)
  if (columns(pos)!=3)
    error("Position must be a 3D vector");
  endif

  ## Calcule la distancia a cada emisor, dada la posición "pos" del
  ## objeto y las posiciones de los emisores.
  
  d=zeros(rows(pos),5);
  
  ## ##################
  ## ## Problema 3.8 ##
  ## ##################
  
  for i = 1:rows(pos)
      #calculo para el primer emisor
    dx = pos(i, 1) - emisorPos(1, 1);
    dy = pos(i, 2) - emisorPos(2, 1);
    dz = pos(i, 3) - emisorPos(3, 1);
    di = sqrt(dx*dx+dy*dy+dz*dz);
    d(i, 1) = di;
    
        #calculo para el segundo emisor
    dx = pos(i, 1) - emisorPos(1, 2);
    dy = pos(i, 2) - emisorPos(2, 2);
    dz = pos(i, 3) - emisorPos(3, 2);
    di = sqrt(dx*dx+dy*dy+dz*dz);
    d(i, 2) = di;
    
        #calculo para el tercer emisor
    dx = pos(i, 1) - emisorPos(1, 3);
    dy = pos(i, 2) - emisorPos(2, 3);
    dz = pos(i, 3) - emisorPos(3, 3);
    di = sqrt(dx*dx+dy*dy+dz*dz);
    d(i, 3) = di;
    
        #calculo para el cuarto emisor
    dx = pos(i, 1) - emisorPos(1, 4);
    dy = pos(i, 2) - emisorPos(2, 4);
    dz = pos(i, 3) - emisorPos(3, 4);
    di = sqrt(dx*dx+dy*dy+dz*dz);
    d(i, 4) = di;
    
        #calculo para el quinto emisor
    dx = pos(i, 1) - emisorPos(1, 5);
    dy = pos(i, 2) - emisorPos(2, 5);
    dz = pos(i, 3) - emisorPos(3, 5);
    di = sqrt(dx*dx+dy*dy+dz*dz);
    d(i, 5) = di;
    
  endfor
  
endfunction


## #################################
## ## Pruebas con varios emisores ##
## #################################

## Pruebe todo el esquema usando el número de emisores dado
## n:    número de emisores que debe ser 3, 4 o 5
## D:    datos de distancia
## E:    posición de emisores
## fig1: número de primera figura a utilizar
function testCase(n,D,E,fig1)

  assert(n>2 && n<=columns(D));

  printf("Probando el caso de %i emisores\n",n);

  ## Prueba unitaria: calcule la posición del primer dato
  p=calcPosition(D(1,1:n),E);
  
  ## Calcule las posiciones a partir de las distancias
  p=calcPositions(D(:,1:n),E);

  ## Distancias a posiciones estimadas
  
  ed=calcDistances(p,E);
  
  ## Errores
  err=(ed-D);

  ## Promedio de errores
  error_average = sum(err)/rows(err)

  ## Muestre los errores de las distancias a los primeros tres emisores
  figure(fig1,"name",cstrcat("Errores con ",num2str(n)," emisores"));
  hold off;
  plot(err(:,1),"r");
  hold on;
  plot(err(:,2),"b");
  plot(err(:,3),"k");
  
  ## Grafique las posiciones encontradas
  figure(fig1+1,"name",...
	 cstrcat("Trayectoria estimada con ",num2str(n)," emisores"));
  hold off;
  
  ## Grafique la trayectoria en azul
  plot3(p(:,1),p(:,2),p(:,3),'b');
  xlabel("x");
  ylabel("y");
  zlabel("z");
  
  hold on;
  
  if (n==3)  
    ## Si n==3 hay otros posibles valores para la trayectoria

    ## Use - en la cuadrática
    p=calcPositions(D(:,1:n),E,2);
    plot3(p(:,1),p(:,2),p(:,3),'r');

    ## Use la predicción para auto-seleccionar cuál solución usar y
    ## píntela con cuadrados magenta
    p=calcPositions(D(:,1:n),E,3);
    plot3(p(:,1),p(:,2),p(:,3),'ms');
    
    ## A la fuerza use n=5 para ver la trayectoria verdadera
    ## y píntela usando circulos verdes
    p=calcPositions(D(:,1:5),E);
    plot3(p(:,1),p(:,2),p(:,3),'go');
  endif

  ## Grafique los emisores en postes
  plot3([E(1,1:n);E(1,1:n)],...
	[E(2,1:n);E(2,1:n)],...
	[E(3,1:n);zeros(1,n)],...
	"linewidth",3);
  
  plot3(E(1,1:n),E(2,1:n),E(3,1:n),'rx');
endfunction;

## Fuerce probar con 3, 4 y 5 emisores
testCase(3,D,E,1); ## 3 emisores, en las figuras 1 y 2
testCase(4,D,E,3); ## 4 emisores, en las figuras 3 y 4
testCase(5,D,E,5); ## 5 emisores, en las figuras 5 y 6
