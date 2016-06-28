#' Un script para analizar el Centro de Presión desde Cero. Las funciones .ao son para Artoficio
#'
#' Se aplican funciones a listas dentro de listas
#' Función 'cargar.ao' permite cargar todos los archivos dentro de una carpeta (1º nivel)
#' Función 'CoP.ao' permite calcular el desplazamiento, velocidad y Aceleración (2º nivel), en ejes X e Y.
#' Función 'lllaply_p' permite graficar todos los niveles en una matriz m=c(3x2)
#' Función 'lllaply' permite aplicar FUN=fun al 3º nivel de la lista. Argumento 'plot'
#' permite graficar y guardar archivo '.EPS' donde el título y nombre de archivo comienza con fun="". ACF por defecto,
#' se puede utilizar PACF, fdGPH...
#' Función 'lllaply_d' permite diferenciar en base a parámetros de diferenciación 'd_par' obtenidos con función lllaply(x,fdGPH)
#' Función 'lllaply_frac' permite obtener parámetros AR y MA
#' Función 'cortar' elimina o reemplaza datos al comienzo y al final. Se debe ingresar el número de puntos o
#' el tiempo en segundos si se entrega la frecuencia de muestreo.

cargar.ao<-function(i=1){
  if(is.null(i))
    warning("Se cargará el primer archivo de la carpeta")
  archivos <-list.files()
  if(is.null(i)) {
    i<-length(archivos) # Número de sujetos evaluados}
    x <- lapply(archivos,read.table,header=F)
    names(x) <- paste('Sujeto', 1:i, sep = " ")
    output=x
  }else{
    x <- lapply(archivos[1:i],read.table,header=F)
    names(x) <- paste('Sujeto', 1:i, sep = " ")
    output=x
  }}


cargar_tot.ao<-function(i=1,freq=40){
  if(missing(i))
    warning("Se cargará el primer archivo de la carpeta")
  if(missing(freq))
    warning("Se utiliza una frecuencia de 40 hz para calcular derivadas")
  archivos <-list.files()
  if(is.null(i)) {
  n<-length(archivos) # Número de sujetos evaluados}
  CoP_orig <- lapply(archivos,read.table,header=F)
  names<-c("A","B","C","D")
  CoP<-lapply(CoP_orig, function(x){colnames(x)<-names
                                    return(x)  })
  CoP_xy<-lapply(CoP_orig,FUN=CoP) #Aplica función CoP
  names(CoP_xy) <- paste('Sujeto', 1:i, sep = " ")
  output=CoP_xy
}else{
  CoP_orig <- lapply(archivos[1:i],read.table,header=F)
  names<-c("A","B","C","D")
  CoP<-lapply(CoP_orig, function(x){colnames(x)<-names
                                    return(x)  })
  CoP_xy<-lapply(CoP_orig,function(x){
    dt=1/freq
    Peso=mean((x[,1]+x[,2]+x[,3]+x[,4])/dim(x[1])[2])
    a<-(x[,3]+x[,4]-x[,1]-x[,2])/Peso
    b<-(x[,1]+x[,3]-x[,2]-x[,4])/Peso
    Eje_X<-(a-mean(a))
    Eje_Y<-(b-mean(b))
    Velocidad_X<-diff(Eje_X)/dt
    Velocidad_Y<-diff(Eje_Y)/dt
    Aceleración_X<-diff(Velocidad_X)/dt
    Aceleración_Y<-diff(Velocidad_Y)/dt
    output=list(Desplazamiento=data.frame("AP"=Eje_Y,"ML"=Eje_X),
                Velocidad=data.frame("AP"=Velocidad_Y,"ML"=Velocidad_X),
                Aceleración=data.frame("AP"=Aceleración_Y,"ML"=Aceleración_X))
  }
  ) #Aplica función CoP
  names(CoP_xy) <- paste('Sujeto', 1:i, sep = " ")
  output=CoP_xy
}}
CoP.ao<-function(x,freq=40){
  if(missing(freq))
    warning("Se utiliza una frecuencia de 40 hz para calcular derivadas")
  dt=1/freq
  Peso=mean((x[,1]+x[,2]+x[,3]+x[,4])/dim(x[1])[2])
  a<-(x[,3]+x[,4]-x[,1]-x[,2])/Peso
  b<-(x[,1]+x[,3]-x[,2]-x[,4])/Peso
  Eje_X<-(a-mean(a))
  Eje_Y<-(b-mean(b))
  Velocidad_X<-diff(Eje_X)/dt
  Velocidad_Y<-diff(Eje_Y)/dt
  Aceleración_X<-diff(Velocidad_X)/dt
  Aceleración_Y<-diff(Velocidad_Y)/dt
  output=list(Desplazamiento=data.frame("AP"=Eje_Y,"ML"=Eje_X),
              Velocidad=data.frame("AP"=Velocidad_Y,"ML"=Velocidad_X),
              Aceleración=data.frame("AP"=Aceleración_Y,"ML"=Aceleración_X))
}
lllaply_p<-function(x,m=c(3,2),f=40,col.main="blue",type="l",...){
  if(!is.null(m))
    warning("Graficado en m=c(3,2)")
  if(!is.null(dt))
    warning("Se ha utilizado una frecuencia de sampleo de 40Hz")
  library(quantmod)
  dt=1/f # Intervalo de registro
  lapply(names(x), function(i) {
    setEPS()
    postscript(paste(i,".eps",sep=""))
    par(mfrow = m)
    lapply(names(x[[i]]), function(d)
      lapply(names(x[[i]][[d]]), function(e)
        (plot(x=seq(from=1,length.out=length(x[[i]][[d]][[e]]),by=dt),y=x[[i]][[d]][[e]],
              xlab="Tiempo (s)",col.main="blue", main=paste(i,", Eje ",e,sep=""),ylab=d,col.main=col.main,type=type))
      ));dev.off();})}


# ACF, donde acf.CoP_xy[[i]] tiene la información del i-ésimo sujeto [[i]][[d]][[,1=x,2=y]],[[i]][[,x,y]]$acf. Usar names y as.numeric para obtener los resultados
lllaply<-function(x,plot=NULL,FUN=acf,fun="acf",...){
  library(fracdiff)
  sapply(names(x), function(i)(
    if((!is.null(plot))){
      setEPS()
      postscript(paste(fun,i,".eps",sep=""))
      par(mfrow = c(3, 2))
      sapply(names(x[[i]]), function(d){
        sapply(names(x[[i]][[d]]),USE.NAMES = TRUE,simplify=FALSE ,function(e)
          FUN(x[[i]][[d]][[e]],plot=(!is.null(plot)), main=paste(fun,i,", Eje ",e,sep=""),
              xlab="Rezago",ylab=paste(fun,d,sep=" ")))})
      ;dev.off()
    } else {
      sapply(names(x[[i]]), function(d)
        sapply(names(x[[i]][[d]]), function(e)
          FUN(x[[i]][[d]][[e]],...),USE.NAMES = TRUE,simplify=FALSE)
        ,USE.NAMES = TRUE,simplify=FALSE)})
    ,USE.NAMES = TRUE,simplify=FALSE)
}

lllaply_d<-function(x,d_par){
  library(fracdiff)
  sapply(names(x),USE.NAMES = TRUE,simplify=FALSE,function(i)(
    sapply(names(x[[i]]),USE.NAMES = TRUE,simplify=FALSE, function(d)(
      sapply(names(x[[i]][[d]]),USE.NAMES = TRUE,simplify=FALSE ,function(e){
        diffseries(x[[i]][[d]][[e]],d_par[[i]][[d]][[e]][[1]])})))))}

# Para obtener los parámetros AR,MA y d
lllaply_frac<-function(x,nar,nma){
  sapply(names(x),USE.NAMES = TRUE,simplify=FALSE,function(i)(
    sapply(names(x[[i]]),USE.NAMES = TRUE,simplify=FALSE, function(d)(
      sapply(names(x[[i]][[d]]),USE.NAMES = TRUE,simplify=FALSE ,function(e){
        fracdiff(x[[i]][[d]][[e]],nar,nma)})))))}

lllaply_frac<-function (x, nar, nma,...) {
  sapply(names(x), USE.NAMES = TRUE, simplify = FALSE, function(i) (
    sapply(names(x[[i]]), USE.NAMES = TRUE, simplify = FALSE, function(d) (
      sapply(names(x[[i]][[d]]),USE.NAMES = TRUE, simplify = FALSE, function(e) {
        fracdiff(x[[i]][[d]][[e]], nar, nma)})))))}

lllaply_arima<-function (x, order,...){
  sapply(names(x), USE.NAMES = TRUE, simplify = FALSE, function(i) 
    (sapply(names(x[[i]]), USE.NAMES = TRUE, simplify = FALSE, function(d) 
      (sapply(names(x[[i]][[d]]), USE.NAMES = TRUE, simplify = FALSE, function(e) {
        arima(x[[i]][[d]][[e]],order,...)})))))}

cortar<-function(x,i=0,f=0,freq=1,r.na=F){
  if(missing(i)&&missing(f))
    stop("No hay datos a eliminar")
  if(missing(freq))
    warning("Frecuencia por defecto son 1hz")
  if(r.na==F){
    if(is.data.frame(x)){
      if(!missing(i)){
        if(!missing(f)){
          x[-c(1:(i*freq),(dim(x)[1]-(f*freq)+1):(dim(x)[1])),n]}
          else
          {x[-c(1:(i*freq)),]}
      }
      else{x[-c((dim(x)[1]-(f*freq)+1):(dim(x)[1])),]}}
    else{
      if(!missing(i)){
        if(!missing(f)){
        x[-c(1:(i*freq),(length(x)-(f*freq)+1):(length(x)))]}
        else
        {x[-c(1:(i*freq))]}}
  else{x[-c((length(x)-(f*freq)+1):(length(x)))]}}}

#ra.na1
else{
  if(is.data.frame(x)){
    if(!missing(i)){
      if(!missing(f)){
        inicio<-rep(NA,times=i*freq)
        fin<-rep(NA,times=f*freq)
        c(inicio,x[c((i*freq)+1):(dim(x)[1]-(f*freq)),],fin) # Con esta linea relleno con NA
      }else
      {
        inicio<-rep(NA,times=i*freq)
        c(inicio,x[c((i*freq)+1):(dim(x)[1]-(f*freq)),]) # Con esta linea relleno con NA    
      }#missing(f)
    }#missnig(i)
    else{
      fin<-rep(NA,times=f*freq)
      c(x[c((i*freq)+1):(dim(x)[1]-(f*freq)),],fin) # Con esta linea relleno con NA
    }}#data.frame1
  else{
    if(!missing(i)){
      if(!missing(f)){
        inicio<-rep(NA,times=i*freq)
        fin<-rep(NA,times=f*freq)
        c(inicio,x[c((i*freq)+1):(length(x)-(f*freq))],fin) # Con esta linea relleno con NA
      }else
      {
        inicio<-rep(NA,times=i*freq)
        c(inicio,x[c((i*freq)+1):(length(x)-(f*freq))]) # Con esta linea relleno con NA    
      }#missing(f)
    }#missnig(i)
    else{
      fin<-rep(NA,times=f*freq)
      c(x[c((i*freq)+1):(length(x)-(f*freq))],fin) # Con esta linea relleno con NA
    }}#data.frame2
}#ra.na2
}#function

sig.fracdiff<-function(x){
  (1-pnorm(abs(unlist(x[c('d','ar','ma')])/
                 unlist(sqrt(diag(x$covariance.dpq))))))*2}

sig.arima<-function(x){
  (1-pnorm(abs(unlist(x$coef)/
                 unlist(sqrt(diag(x$var.coef))))))*2}

count.p<-function(x){
  suma<-sum(x$acf>=qnorm((1 + 0.95)/2)/sqrt(x$n.used))
  }
# Extraer información desde modelo ajustado (residuos, coef...) luego aplicar Normalidad, ACF, PACF, etc...
# res_d<-lllaply(par_d,FUN=acf,fun="acf")
# coef = do.call("rbind", sapply(names(par_d), USE.NAMES = TRUE,simplify=T,
#                    function(i) sapply(names(par_d[[i]]), USE.NAMES = TRUE,simplify=T,
#                                       function(d) lapply(par_d[[i]][[d]], "[[", 'coef'))))
#resid = sapply(names(par_d), USE.NAMES = TRUE,simplify=F,
#                               function(i) sapply(names(par_d[[i]]), USE.NAMES = TRUE,simplify=F,
#                                                  function(d) lapply(par_d[[i]][[d]], "[[", 'residuals')))

# Coeficientes ordenados por derivada y eje para luego aplicar estadística
#pos<-seq(from=1, to=n*6,by=6)
#coef<-list(dx=coef[c(pos),],dy=coef[c(pos+1),],vx=coef[c(pos+2),],vy=coef[c(pos+3),],ax=coef[c(pos+4),],ay=coef[c(pos+5),])

#d_CoP_xy<-sapply(names(CoP_xy), function(i)
#  sapply(names(CoP_xy[[i]][[1]]), function(e){
#            diffseries(CoP_xy[[i]][[1]][e],d=d[[i]][[e]])},
#            USE.NAMES = TRUE,simplify=FALSE)
#  ,USE.NAMES = TRUE,simplify=FALSE)


#Box-test para ver efecto GARCH
#b_xy<-lapply(CoP_xy, function(x) lapply(x, function(x) lapply(x, function(x) Box.test(x)$p.value)))
#b_xy<-matrix(unlist(b_xy),ncol=n,byrow=T)
#colnames(b_xy) <- paste('Sujeto', 1:n, sep = " ")

#Cálculo del Área del Centro de Presión. Se puede resumir con function(x)...
#cov_mat<-lapply(CoP_xy,function(x) lapply(x,FUN=cov)$Desplazamiento) # cov_mat<-(cov(cbind(x,y)))
#autov<-lapply(cov_mat, function(x){eigen(x)$values}) # autov<-(eigen(cov_mat))$values
#area<-lapply(autov, function(x){(pi*prod(2.4478*sqrt(x)))})  # area<-(pi*prod(2.4478*sqrt(autov)))

# comparar con función por defecto
#me <- apply((cbind(x,y)), 2, mean)
#v <- var(cbind(x,y))
#rad <- sqrt(2*qf(0.5, 2, nrow(cbind(x,y))-1))
#z <- ellipse(me, v, rad, segments=1001)
#dist2center <- sqrt(rowSums((t(t(z)-me))^2))
#pi*min(dist2center)*max(dist2center)

# Función para cortar los primeros y últimos 5 segundos
#x<-(lllaply(x,FUN=cortar,fun="cortar",i=5,f=5,freq=40))


save.excel <-function(.list, default = 'var1', path = ''){
  require("XLConnect")
  .name <- as.list(match.call())[2]
  if(is.language(.name[[1]])) wb_name <- paste0(paste0(path, default, collapse = '/'), '.xlsx')
  if(is.symbol(.name[[1]])) wb_name <- paste0(paste0(path, as.character(.name), collapse = '/'), '.xlsx')
  wb <- loadWorkbook(wb_name, create = TRUE)
  createSheet(wb, names(.list))
  writeWorksheet(wb,.list, names(.list),header=FALSE)
  saveWorkbook(wb)
}

cargar<-function (i = 1, freq = 100,col.names=NULL,header=F,...) 
{
  if (missing(i)) 
    warning("Se cargará el primer archivo de la carpeta")
  if (missing(freq)) 
    warning("Se utiliza una frecuencia de 100 hz para calcular derivadas")
  archivos <- list.files()
  if (is.null(i)) {
    dat <- sapply(archivos, function(x){read.table(x,header=header,...)},
                  USE.NAMES = TRUE,simplify=FALSE)
  nombre<-lapply(strsplit(names(dat),"\\."), function(x) paste(x[1:(length(x)-1)], collapse="."))#    lappl
  names(dat)<-nombre}
  else {
    dat <- sapply(archivos[1:i], function(x){read.table(x,header=header,...)},
                  USE.NAMES = TRUE,simplify=FALSE)  
    nombre<-lapply(strsplit(names(dat),"\\."), function(x) paste(x[1:(length(x)-1)], collapse="."))#    lappl
    names(dat)<-nombre}
  if (!is.null(col.names)&&header==F) {
    dat<-lapply(dat,function(x){names(x)<-col.names
                                return(x)})

    }else{return(dat)}
  dat<-append(dat,list("Sampling"=(data.frame(freq,dt = 1/freq))))
}

d.t<-function(x,d=1,cols=NULL){
  if((d>=3)||(d<=0))
    stop("Solo se pueden realizarán las 2 primeras derivadas")
  n<-length(x)
  if(missing(cols))
{
  dx.dt<-sapply(x[1:(n-1)],USE.NAMES = TRUE,simplify=FALSE, function(i) {apply(i,function(j){diff(j)/(x[["Sampling"]]$dt)},MARGIN=2)})
  if(d==1){
    return(list(Vel=dx.dt))}
  else
  {vx.dt<-sapply(x[1:(n-1)],USE.NAMES = TRUE,simplify=FALSE, function(i){apply(i,function(j){diff(j,differences=2)/(x[["Sampling"]]$dt)},MARGIN=2)})}
  return(list(Vel=dx.dt,Acc=vx.dt))}
else
{
  dx.dt<-sapply(x[1:(n-1)],USE.NAMES = TRUE,simplify=FALSE, function(i) {apply(i[cols],function(j){diff(j)/(x[["Sampling"]]$dt)},MARGIN=2)})
  if(d==1){
    return(list(Vel=dx.dt))}
  else
  {vx.dt<-sapply(x[1:(n-1)],USE.NAMES = TRUE,simplify=FALSE, function(i) {apply(i[cols],function(j){diff(j,differences=2)/(x[["Sampling"]]$dt)},MARGIN=2)})}
  return(list(Vel=dx.dt,Acc=vx.dt))}
}

# Cálculo del Área del Centro de Presión.
area<-function(x,y){
#cov_mat<-lapply(CoP_xy,function(x) lapply(x,FUN=cov)$Desplazamiento)
cov_mat<-(cov(cbind(x,y)))
#autov<-lapply(cov_mat, function(x){eigen(x)$values})
autov<-(eigen(cov_mat))$values
#area<-lapply(autov, function(x){(pi*prod(2.4478*sqrt(x)))})
area<-(pi*prod(2.4478*sqrt(autov)))
area<-pi*r^2
r<-sqrt(area/pi)#?
return(area)
}

# Función para graficar área del CoP
plot.a<-function(x,y,...){
a<-area(x,y)
r<-sqrt(a/pi)#?
theta <- seq(0, 2 * pi, length = 200)
xlim=c(min(((r * cos(theta))+mean(x)),x),max(((r * cos(theta))+mean(x)),x))
ylim=c(min(((r * sin(theta))+mean(y)),y),max(((r * sin(theta))+mean(y)),y))
plot(x,y,type="l",
     xlim=xlim,ylim=ylim,...)
lines(x = (r * cos(theta))+mean(CoPx), y = (r * sin(theta))+mean(CoPy),col="red")
}

# Función para el FFT
fft.<-function(x,samplingRate){
  library(tuneR)
  n <- length(x)
  p <- fft(x)
  nUniquePts <- ceiling((n+1)/2)
  p <- p[1:nUniquePts] #select just the first half since the second half  is a mirror image of the first
  p <- abs(p)  #take the absolute value, or the magnitude 
  p <- p / n #scale by the number of points so that the magnitude does not depend on the length of the signal or on its sampling frequency  
  p <- p^2  # square it to get the power 
  # multiply by two (see technical document for details)
  # odd nfft excludes Nyquist point
  if (n %% 2 > 0){
    p[2:length(p)] <- p[2:length(p)]*2 # we've got odd number of points fft
  } else {
    p[2: (length(p) -1)] <- p[2: (length(p) -1)]*2 # we've got even number of points fft
  }
  freqArray <- (0:(nUniquePts-1)) * (as.numeric(samplingRate) / n) #  create the frequency array 
  plot(freqArray, p, type='l', col='black', xlab='Frequency (Hz)', ylab='Power',xlim=c(0,10),ylim=c(0,max(p[10:4000])))
}
