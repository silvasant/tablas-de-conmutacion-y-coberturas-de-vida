getwd()
setwd(dir=choose.dir())
getwd()


#Defino la funcion para generar la tabla de valores de conmutacion
crear.tabla<-function(tabla,tasa=0.04,lx.init=10000000){
  
  X<-tabla[,1]
  qx<-tabla[,2]
  lx<-rep(0,length(X))
  dx<-rep(0,length(X))
  lx[1]<-lx.init
  Dx<-rep(0,length(X))
  Cx<-rep(0,length(X))
  Nx<-rep(0,length(X))
  Mx<-rep(0,length(X))
  Rx<-rep(0,length(X))
  Sx<-rep(0,length(X))
  px1<-1-qx
  
  
  
  v<-1/(1+tasa)
  
  q1<-function(x){
    q1<-match(x,tabla[,2])  
    return(q1)}
  
  d<-function(x){                
    q1(x)*l(x)
  }
  
  l<-function(x){
    li<-match((x-1),lx)
    px1i<-match(x,px1)                    
    l<-li*px1i
  }
  
  lx[1]<-lx.init                       
  dx[1]<-lx[1]*qx[1]                   
  for (i in 2:length(X)) {
    lx[i]<-lx[i-1]-dx[i-1]
    dx[i]<-round(qx[i]*lx[i],0)  
  }                                  
  
  
  for (i in 1:length(X)) {
    Dx[i]<-round(lx[i]*v^(i-1),0)
    Cx[i]<-round(dx[i]*v^(i),0)
  }  
  
  for(i in X){
    Nx[i]<-sum(Dx[i:length(X)])
    Mx[i]<-sum(Cx[i:length(X)])
  }
  
  for(i in X){
    Sx[i]<-sum(Nx[i:length(X)])
    Rx[i]<-sum(Mx[i:length(X)])
  }
  
  Nx[length(X)]<-Dx[length(X)];Sx[length(X)]<-Dx[length(X)] 
  Mx[length(X)]<-Cx[length(X)];Rx[length(X)]<-Cx[length(X)] 
  
  
  return(tabla.nueva<-data.frame(cbind(X,qx,lx,dx,Dx,Nx,Sx,Cx,Mx,Rx)))
}

#############################

### Funciones de búsqueda ###

D<-function(x,tabla){
  D<-tabla[(x+1),5]
  return(D)}
N<-function(x,tabla){
  N<-tabla[(x+1),6]
  return(N)}
S<-function(x,tabla){
  S<-tabla[(x+1),7]
  return(S)}
C<-function(x,tabla){
  C<-tabla[(x+1),8]
  return(C)}
M<-function(x,tabla){
  M<-tabla[(x+1),9]
  return(M)}
R<-function(x,tabla){
  R<-tabla[(x+1),10]
  return(R)}
### Coberturas de capitales constantes ###
a<-function(x,h,n,tabla){
  o<-x+h+n
  if(o>=length(tabla[,1])){
    a<-N(x+h,tabla)/D(x,tabla)
    }
  else {
    a<-(N(x+h,tabla)-N(o,tabla))/D(x,tabla)   
    }
  return(a)
  }                                        

A<-function(x,h,n,tabla){
  o<-x+h+n
  if(o>=length(tabla[,1])){
    A<-M(x+h,tabla)/D(x,tabla)
    }
  else{
    A<-(M(x+h,tabla)-M(o,tabla))/D(x,tabla)
    }
  return(A)
 }                               


E<-function(x,h,tabla){
  E<-D(x+h,tabla)/D(x,tabla)
  return(E)
}

### Coberturas de capitales variables ###

aI<-function(x,h,n,tabla){
  o<-x+h+n
  if(o>=length(tabla[,1])){
    aI<-S(x+h,tabla)/D(x,tabla)
  }else { aI<-(S(x+h,tabla)-S(o,tabla)-n*N(o,tabla))/D(x,tabla)}
  return(aI)
}


AI<-function(x,h,n,tabla){
  o<-x+h+n
  if(o>=length(tabla[,1])){
    AI<-R(x+h,tabla)/D(x,tabla)
  }else { AI<-(R(x+h,tabla)-R(o,tabla)-n*M(o,tabla))/D(x,tabla)}
  return(AI)
}

av<-function(x,h,n,r,tabla){
  av<-a(x,h,n,tabla)+r*aI(x,h+1,n-1,tabla)			
  return(av)
  }

Av<-function(x,h,n,r,tabla){
  Av<-A(x,h,n,tabla)+r*AI(x,h+1,n-1,tabla)			
  return(Av)
  }

avg<-function(x,h,n,r,tabla,tasa=0.04){				
  if(r<tasa){											
    v_1<-(1+r)/(1+tasa)
    tasa_2<-(1/v_1)-1
    tabla_2<-crear.tabla(tabla,tasa=tasa_2,lx.init=10000000)	
    avg<-a(x,h,n,tabla=tabla_2)
  }	else{
    if(r==tasa){
      avg<-sum(tabla[(x+h+1):(x+h+n),3])/tabla[(x+1),3]			
    } else {print("r>tasa")}
  }												
  return(avg)
}


Avg<-function(x,h,n,r,tabla,tasa=0.04){				
  if(r<tasa){											
    v_1<-(1+r)/(1+tasa)
    tasa_2<-(1/v_1)-1
    tabla_2<-crear.tabla(tabla,tasa=tasa_2,lx.init=10000000)	
    Avg<-A(x,h,n,tabla=tabla_2)
  }	else{
    if(r==tasa){
      v<-1/(1+tasa)
      Avg<-v*sum(tabla[(x+h+1):(x+h+n),2])				
    } else {print("r>tasa")}
  }												
  return(Avg)
}

##############################################################

#Obtengo tabla de conmutacion a usar como referencia
CSO1980<-read.csv2("Tabla CSO1980.csv",sep=";",skip=2,header=T);View(CSO1980)
data.class(CSO1980)
#Construyo init.data , los datos iniciales a partir de los cuales voy a construir la tabla de prueba, 
#a partir de los datos inicales de la tabla de control
qx<-CSO1980[,11]
X<-CSO1980[,1]
init.data<-cbind(X,qx)
#Genero una tabla:
CSO80_valconm<-crear.tabla(tabla = init.data,lx.init=10000000)
View(CSO80_valconm)

###Calculo de coberturas ###

E(35,10,CSO80_valconm)
a(35,5,5,CSO80_valconm)
A(35,5,5,CSO80_valconm)
aI(35,5,5,CSO80_valconm)
AI(35,5,5,CSO80_valconm)
av(35,5,5,5,CSO80_valconm)
Av(35,5,5,5,CSO80_valconm)
avg(35,0,10,0.03,CSO80_valconm)
Avg(35,0,10,0.03,CSO80_valconm)
avg(35,0,5,0.04,CSO80_valconm)
Avg(35,0,5,0.04,CSO80_valconm)


### Graficos ###

edades<-35:45

rtas.vida<-rep(NA,10)
  for (i in edades){
  rtas.vida[i-34]<-a(i,0,10,CSO80_valconm)
}

temp.muerte<-rep(NA,10)
for (i in edades){
  temp.muerte[i-34]<-A(i,0,10,CSO80_valconm)
}


par(mfrow=c(1,2))
plot(edades,rtas.vida,xlim = c(35,45),col="red",xlab = "x: Edad al momento de contratación",ylab = "a(x,0,10)",main = "Cobertura de vida")
plot(edades,temp.muerte,xlim = c(35,45),col="blue",xlab = "x: Edad al momento de contratación",ylab = "A(x,0,10)",main = "Cobertura de muerte")


variacion.plazo<-2:10
rtas.vida.plazo<-rep(NA,9)
for (i in variacion.plazo){
  rtas.vida.plazo[i-1]<-a(35,0,i,CSO80_valconm)
};rtas.vida.plazo

temp.muerte.plazo<-rep(NA,9)
for (i in variacion.plazo){
  temp.muerte.plazo[i-1]<-A(35,0,i,CSO80_valconm)
};temp.muerte.plazo


par(mfrow=c(1,2))
plot(variacion.plazo,rtas.vida.plazo,xlim = c(2,10),col="red",xlab = "n: plazo de cobertura",ylab = "a(35,0,n)",main = "Cobertura de vida")
plot(variacion.plazo,temp.muerte.plazo,xlim = c(2,10),col="blue",xlab = "n: plazo de cobertura",ylab = "A(35,0,n)",main = "Cobertura de muerte")


variacion.tasa<-0.01
rtas.vida.tasa<-rep(NA,10)
for (i in 1:10){
  rtas.vida.tasa[i]<-aI(35,5,10,crear.tabla(tabla = init.data,tasa = i*variacion.tasa))
};rtas.vida.tasa
temp.muerte.tasa<-rep(NA,10)
for (i in 1:10){
  temp.muerte.tasa[i]<-AI(35,5,10,crear.tabla(tabla = init.data,tasa = i*variacion.tasa))
};temp.muerte.tasa

par(mfrow=c(1,2))
plot((1:10)*variacion.tasa,rtas.vida.tasa,xlim = c(0.01,0.1),col="red",xlab = "i: tasa de interés",ylab = "aI(35,5,10)",main = "Cobertura de vida")
plot((1:10)*variacion.tasa,temp.muerte.tasa,xlim = c(0.01,0.1),col="blue",xlab = "i: tasa de interés",ylab = "AI(35,5,10)",main = "Cobertura de muerte")




