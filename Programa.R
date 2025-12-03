
# -----------------------------------------------------------------------------#
# INICIO  ----
# -----------------------------------------------------------------------------#

# -----------------------------------------------------------------------------#
# Evaluación pronósticos fuera de muestra y
# Creacion de pronosticos condicionados a partir de modelo ARIMAX
# Siguiendo la metodologia de Guerrero (1998, 2001)
# Adaptacion: Arturo Gonzalez ----
# -----------------------------------------------------------------------------#



# carga de librerias (funciones para el desarrollo practico)
require(pacman) # library(pacman)

p_load(readxl, sandwich, car, lmtest, TSstudio, lmtest, forecast
       , tseries, TSA, tsoutliers, GGally, xts, ggplot2, dplyr
       , MASS, nortest, tsoutliers, astsa, Metrics, matlib, ltsa
       , bestglm)


## Importacion datos ----
#------------------------------------------------------------------------------#
Base_ent <- read_xlsx("Base_Modelo_ARIMAX.xlsx"
                              , sheet = "Base"
                              , range = "a2:e312", col_names = T)
Base_ent_ts <- ts(Base_ent[,-1],start = c(2000,1),frequency = 12)
tail(Base_ent_ts)


## Elementos gráficos ----
#------------------------------------------------------------------------------#

windows()
Base_ent_ts %>% autoplot(facets=T) + xlab("Fecha") + ylab("Valores")
Base_ent_ts_lx <- log(Base_ent_ts)


# Transformaciones ----
#------------------------------------------------------------------------------#

# diferencia log mes -> es simil con la inflación mensual
dlx_IPC <- diff(log(Base_ent_ts[,1]), lag = 1, differences = 1)

# Aceleracion mensual de la inflacion mensual
ddlx_IPC <- diff(dlx_IPC, lag = 1, differences = 1)
diff(log(Base_ent_ts[,1]), lag = 1, differences = 2)
# diferencia log anual -> simil inflacion anual
slx_IPC <- diff(log(Base_ent_ts[,1]), lag = 12, differences = 1)

# Aceleracion mensual inflacion anual
dslx_IPC <- diff(slx_IPC, lag = 1, differences = 1)

# Aceleracion anual inflacion anual
sslx_IPC <- diff(slx_IPC, lag = 12, differences = 1)


# Test de estacionariedad ----
#------------------------------------------------------------------------------#
adf.test(Base_ent_ts[,1]) 
kpss.test(Base_ent_ts[,1])
pp.test(Base_ent_ts[,1])  
# se necesita una diferenciar para lograr la estacionariedad

adf.test(diff(Base_ent_ts[,1]))
kpss.test(diff(Base_ent_ts[,1]))
pp.test(diff(Base_ent_ts[,1]))  

ndiffs(Base_ent_ts[,1])      # -> d=2
ndiffs(log(Base_ent_ts[,1])) # -> d=1 en log
# diferenciado una vez se logra la estacionariedad


# Identificar estructura ----
#------------------------------------------------------------------------------#
windows()
tsdisplay(diff(log(Base_ent_ts[,1])))
ts_cor(diff(log(Base_ent_ts[,1])))
nsdiffs(diff(log(Base_ent_ts[,1]))) # indica el valor del segmento estacional "D"
ts_seasonal(diff(log(Base_ent_ts[,1])),type = "all")
eacf(diff(log(Base_ent_ts[,1])))



# Primer modelo modelo IPC sin Exogenas
p <- 4 # valor para la parte autoregresiva
d <- 1 # num. de diferencia que voy aplicar internamente para lograr estacionariedad
q <- 2 # valor para la parte de medias moviles (memoria)
P <- 1 # valor para la parte autoregresiva estacional
D <- 0 # num. de diferencias estacionales para lograr estacionariedad (parte estacional)
Q <- 0 # valor parte de medias moviles estacionales
mod1 <- Arima(log(Base_ent_ts[,1]),order = c(p,d,q),seasonal = c(P,D,Q), method = "ML")
summary(mod1)
lmtest::coeftest(mod1)
checkresiduals(mod1)
shapiro.test(mod1$residuals)


# segundo modelo incorporando exógenas
mod2 <- auto.arima(log(Base_ent_ts[,1]),d=1, D = 0, stationary = F
                   , stepwise = T, xreg = log(Base_ent_ts[,-1])
                   , trace = T)
summary(mod2)
lmtest::coeftest(mod2)
checkresiduals(mod2)
shapiro.test(mod2$residuals)



# Pronos exogenas ----
#------------------------------------------------------------------------------#
exogenas_pry <- ts(read_xlsx("Base_Modelo_ARIMAX.xlsx"
                             , sheet = "Exogenas",range = "b3:d17")
                   , start = c(2025,11), frequency = 12)
exogenas_pry_lx <- log(exogenas_pry)
dim(exogenas_pry_lx)

# Pronosticos libres -----
# -----------------------------------------------------------------------------#

# Pronostico libre modelo 1
horizonte <- 14
fore_mod1 <- forecast(mod1,h=horizonte)
windows()
autoplot(fore_mod1)
exp(fore_mod1$mean)

windows()
plot(fore_mod1, main = "Pronos libres log IPC mod1")
lines(mod1$fitted,col="green")

# Pronstico libre modelo 2
fore_mod2 <- forecast(mod2, xreg = exogenas_pry_lx)
exp(fore_mod2$mean)

windows()
plot(fore_mod2,main="Pronos libres log IPC mod2")
lines(mod2$fitted,col="red")



# Eval fuera de muestra ----
# -----------------------------------------------------------------------------#


# -----------------------------------------------------------------------------#

# Definicion de las matrices de errores RMSE para cada estructura
# 1. modelo arimax estimado criterio propio (mod1)
# 2. modelo arimax estimado algoritmo auto.arima() y con exogenas (mod2)


# datos de entrada
base_in <- na.omit(Base_ent_ts_lx)
exogenas <- base_in[,-1]

recorte_fm <- (nrow(base_in)-24)/nrow(base_in)   # 92% de los datos
horizonte  <- 12       # maximo horizonte de pronostico (pasos de pronostico)
pronostico <- 24       # simulaciones por cada horizonte de pronostico

rmse_err_mod1  <- rmse_err_mod2 <- matrix(NA, horizonte,pronostico)

rownames(rmse_err_mod1) <- c("err_1paso","err_2paso","err_3paso","err_4paso"
                             ,"err_5paso","err_6paso","err_7paso","err_8paso"
                             ,"err_9paso","err_10paso","err_11paso","err_12paso")
colnames(rmse_err_mod1) <- paste("Simul",seq(1:pronostico),sep = "_")
 

# Errores modelo propio mod1 ----
# -----------------------------------------------------------------------------#
i <- j <- 1
for(j in 1:horizonte){
  base_train <- ts(base_in[(1:round(nrow(base_in)*recorte_fm)),]
                   , start = c(2000,1),frequency = 12)
  for(i in 1:(pronostico-j+1)) {
    mod1_sim <- Arima(base_train[,1], order = c("Valor p","Valor d","Valor q"),seasonal = c("Valor P","Valor D","Valor Q")
                      , xreg = base_train[,-1]
                      , method="CSS")  
    
    if(j==1){pronos_ind <- t(exogenas[(nrow(as.matrix(base_train))+1):(nrow(as.matrix(base_train))+j),])}
    else {pronos_ind <- exogenas[(nrow(as.matrix(base_train))+1):(nrow(as.matrix(base_train))+j),]}
    
    fore_mod1_sim <- forecast(mod1_sim
                              , xreg=pronos_ind)
    #calculo del error
    rmse_err  <- rmse(base_in[(nrow(as.matrix(base_train))+1):(nrow(as.matrix(base_train))+j),1],fore_mod1_sim$mean)
    rmse_err_mod1[j,i] <- rmse_err
    base_train <- rbind(as.matrix(base_train),base_in[(nrow(as.matrix(base_train))+1):(nrow(as.matrix(base_train))+1),])
  }
  print(t(rmse_err_mod1))
}


View(t(rmse_err_mod1))

write.csv(t(rmse_err_mod1),"Datos_sal/estructura_errores_mod1.csv")
openxlsx::write.xlsx(t(rmse_err_mod1),"Datos_sal/estructura_errores_mod1.xlsx")

# Errores modelo propio mod2 ----
# -----------------------------------------------------------------------------#
rownames(rmse_err_mod2) <- c("err_1paso","err_2paso","err_3paso","err_4paso"
                             ,"err_5paso","err_6paso","err_7paso","err_8paso"
                             ,"err_9paso","err_10paso","err_11paso","err_12paso")
colnames(rmse_err_mod2) <- paste("Simul",seq(1:pronostico),sep = "_")


j <- i <- 1 
for(j in 1:horizonte){
  base_train <- base_in[(1:round(nrow(base_in)*recorte_fm)),1]
  base_train_exo <- exogenas[(1:round(nrow(exogenas)*recorte_fm)),]
  for(i in 1:(pronostico-j+1)) {
    mod2_sim <- auto.arima(base_train, xreg=base_train_exo,d=1, stationary = F
                           , stepwise = T , trace = F)
    
    if(j==1){pronos_ind <- t(exogenas[(nrow(as.matrix(base_train))+1):(nrow(as.matrix(base_train))+j),])}
    else {pronos_ind <- exogenas[(nrow(as.matrix(base_train))+1):(nrow(as.matrix(base_train))+j),]}
    
    fore_mod2_sim <- forecast(mod2_sim, xreg=pronos_ind)
    
    rmse_err  <- rmse(base_in[(nrow(as.matrix(base_train))+1):(nrow(as.matrix(base_train))+j),1]
                      ,fore_mod2_sim$mean)
    rmse_err_mod2[j,i] <- rmse_err
    base_train <- rbind(as.matrix(base_train),base_in[(nrow(as.matrix(base_train))+1):(nrow(as.matrix(base_train))+1),1])
    base_train_exo <- rbind(base_train_exo,exogenas[(nrow(base_train_exo)+1):(nrow(base_train_exo)+1),])
    #print(rmse_err)
  }
  print(t(rmse_err_mod2))
}


View(t(rmse_err_mod2))

write.csv(t(rmse_err_mod2),"Datos_sal/estructura_errores_mod2.csv")
openxlsx::write.xlsx(t(rmse_err_mod2),"Datos_sal/estructura_errores_mod2.xlsx")


# -----------------------------------------------------------------------------#
# Evaluacion de la estructura de errores fuera de muestra de cada modelo
#
# se crea un vector con el estadistico y el p-valor de la prueba para cada hori-
# zonte de tiempo, desde 1 paso adelante, hasta 12 pasos adelante
# -----------------------------------------------------------------------------#

# Pruebas de Dieblod & Mariano ----
salidaDM <- NULL


i<-1
for(i in 1:horizonte){
  salida <-  cbind("esta. prueba"=dm.test(rmse_err_mod1[i,],rmse_err_mod2[i,], h = i
                                          ,alternative=c("g"),varestimator = "bartlett")$statistic
                   ,"p-valor"=dm.test(rmse_err_mod1[i,],rmse_err_mod2[i,], h = i
                                      ,alternative=c("g"),varestimator = "bartlett")$p.value)
  salidaDM <- rbind(salidaDM,salida)
}
rownames(salidaDM) <- c("1_paso","2_pasos","3_pasos","4_pasos"
                        ,"5_pasos","6_pasos","7_pasos","8_pasos"
                        ,"9_pasos","10_pasos","11_pasos","12_pasos")

write.csv(salidaDM,"Datos_sal/Prueba_Diebold_Mariano_modelos_Arima.csv")
View(salidaDM)

# Para el tema de varianza negativa en la prueba de diebold / Mariano ver
# https://www.sciencedirect.com/science/article/abs/pii/S0169207017300559

# Fin eval fuera de muestra ----/
# -----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#




# -----------------------------------------------------------------------------#
# Metodologia pronos restringidos ----
# implementacion metodologica del paper descrito por Guerrero (1989)
# -----------------------------------------------------------------------------#


# Generacion impulso respuesta ----
irf_mod1a <- FitARMA::ImpulseCoefficientsARMA(phi = mod1$model$phi
                                              ,theta = mod1$model$theta
                                              , lag.max = (dim(exogenas_pry_lx)[1]-1))
irf_mod1a

irf_mod2a <- FitARMA::ImpulseCoefficientsARMA(phi = mod2$model$phi
                                              ,theta = mod2$model$theta
                                              , lag.max = (dim(exogenas_pry_lx)[1]-1))
irf_mod2a


# Matriz diagonal H ----

mat_H <- function(irf_input, tamanho) {
  H_mat <- diag(irf_input[1], nrow = tamanho,ncol = tamanho)
    for (j in 1:tamanho) {
    H_mat[,j] <- c(rep(0,j-1),head(irf_input,tamanho-(j-1)))
    }
  return(H_mat)
}

H_mat1 <- mat_H(irf_mod1a, dim(exogenas_pry_lx)[1])
H_mat2 <- mat_H(irf_mod2a, dim(exogenas_pry_lx)[1])
View(H_mat2)


# Definicion restricciones ----

# Definicion de las restricciones y matriz C

# Objetivo
# Restriccion1: finalizar 2024 la inflacion sea 5%
# Restriccion2: finalizar 2025 con inflacion sea 3%


Z1 <-  base::matrix(fore_mod1$mean, ncol = 1
                    ,nrow = base::length(fore_mod1$mean)
                    ,byrow = T)
num_restri <- 2

C_in <- base::rep(0,num_restri*base::length(fore_mod1$mean))

C_def <- base::matrix(C_in, ncol = num_restri
                      , nrow = base::length(fore_mod1$mean), byrow = TRUE)



base::colnames(C_def) <- base::paste("Rest"
                                     ,base::seq(1:num_restri),sep = "_")
C_def <- xts(C_def
             , order.by = seq.Date(as.Date("2025-06-1"), as.Date("2026-12-1")
                                   , "month") )
  

# indicador de fechas de afectacion ----
# -----------------------------------------------------------------------------#


C_def[7,1] <- 1  # Bandera donde se va a restringir el pronostico
C_def[19,2] <- 1 # Restriccion 2


# Definicion de la restriccion puntual. Esto se hace basado 
# en una fecha de referencia para el caso de la inflacion
# Es decir, se quiere que al cierre de 2025 la inflacion sea 5%.
# entonces valor de log(IPC 2025:12) = log(IPC 2024:12)*1.05
#                   log(IPC 2026:12) = log(IPC 2024:12)*1.05*1.03

valor_log_ref <- tail(Base_ent_ts,6)[1,1]   # valor log de IPC 2024:12

# Definición valores a los que debe estar condicionado el pronóstico
fore_obj <- base::matrix(c(valor_log_ref*1.05, valor_log_ref*1.05*1.03)
                        , ncol = 1, nrow = num_restri
                        , byrow = T)

fore_obj <- log(fore_obj)


# Calculo de pronos restringidos modelo 1 ----
# -----------------------------------------------------------------------------#

res_e <- fore_obj - base::t(C_def)%*%Z1
res_b <- solve(t(C_def) %*% t(H_mat1) %*% H_mat1 %*% C_def)
res_a <- base::t(H_mat1)%*%(H_mat1)%*%C_def%*%res_b

pronos_rest <- fore_mod1$mean + res_a%*%res_e



indice_pry <-  base::matrix(c("IPC"=base::as.matrix(Base_ent_ts[,1])
                              ,"IPC"=base::as.matrix(exp(pronos_rest))))

var_yoy <- base::round((indice_pry[13:nrow(indice_pry)]/indice_pry[1:(nrow(indice_pry)-12)]-1)*100,2)
var_mom <- base::round((indice_pry[2:nrow(indice_pry)]/indice_pry[1:(nrow(indice_pry)-1)]-1)*100,2)

salida1 <- base::cbind("IPC"=indice_pry,"var_mom"=c(NA,var_mom)
                      ,"var_yoy"=c(rep(NA,12),var_yoy))
salida1 <- xts(salida1
              , order.by = seq.Date(as.Date("2000-1-1"), as.Date("2026-12-1")
                                    , "month")) 
salida1
write.csv(salida1,"Datos_sal/pronostico_restringido_mod1.csv")


# Calculo de pronos restringidos modelo 2 ----
# -----------------------------------------------------------------------------#
Z2 <-  base::matrix(fore_mod2$mean, ncol = 1
                    ,nrow = base::length(fore_mod2$mean)
                    ,byrow = T)


res_e <- fore_obj - base::t(C_def)%*%Z2
res_b <- solve(base::t(C_def)%*%base::t(H_mat2)%*%H_mat2%*%C_def)
res_a <- base::t(H_mat2)%*%(H_mat2)%*%C_def%*%res_b

pronos_rest2 <- fore_mod2$mean + res_a%*%res_e


indice_pry2 <-  base::matrix(c("IPC"=base::as.matrix(Base_ent_ts[,1])
                              ,"IPC"=base::as.matrix(base::exp(pronos_rest2))))

var_yoy2 <- base::round((indice_pry2[13:nrow(indice_pry2)]/indice_pry2[1:(nrow(indice_pry2)-12)]-1)*100,2)
var_mom2 <- base::round((indice_pry2[2:nrow(indice_pry2)]/indice_pry2[1:(nrow(indice_pry2)-1)]-1)*100,2)

salida2 <- base::cbind("IPC"=indice_pry2,"var_mom"=c(NA,var_mom2)
                      ,"var_yoy"=c(rep(NA,12),var_yoy2))
salida2 <- xts(salida2
              , order.by = seq.Date(as.Date("2000-1-1"), as.Date("2026-12-1")
                                    , "month")) 
salida2
write.csv(salida2,"Datos_sal/pronostico_restringido_mod2.csv")

# Grafica de comparación
union <- cbind.xts(salida1, salida2)
View(union)
colnames(union) <- c("IPC_mod1", "var_mom_mod1", "var_yoy_mod1", "IPC_mod2", "var_mom_mod2", "var_yoy_mod2")
ts_plot(union)

  
# (por mejorar, no ejecutar) ----
#  Test de compatibilidad de las restricciones con la historia

num_restri

# test_rest <- function(fore_obj, fore_libre, sigma_modelo
#                       , num_restri, H_entra, C_entra, num_datos){
#   Z1 <-  base::matrix(fore_libre, ncol = 1,nrow = base::length(fore_libre), byrow = T)
#   res_e <- fore_obj - base::t(C_entra)%*%Z1
#   res_b <- matlib::inv(base::t(C_entra)%*%base::t(H_entra)%*%H_entra%*%C_entra)
#   res_a <- base::t(H_entra)%*%(H_entra)%*%C_entra%*%res_b
#   H1 <- (1/sigma_modelo)*(t(res_e)%*%res_b%*%res_e)
#   H2 <- (1/num_restri)*(t(res_e)%*%res_b%*%res_e)*H1
# 
#   print("Prueba de compatilidad de las restriciones I")
#   print(paste("Punto_critico: ", round(H2,digits = 4)," "
#               ,"p-valor:",stats::pchisq(H2,num_restri)))
#   
#   print("Prueba de compatilidad de las restriciones II")
#   print(paste("Punto_critico: ", round(H1,digits = 4)," "
#               ,"p-valor:",round(stats::pchisq(H1,num_restri),digits = 4)))
# 
# }

# test_rest(fore_obj = fore_obj, fore_libre = fore_mod1$mean
#                          , sigma_modelo = mod1$sigma2
#                          , num_restri = num_restri 
#                          , H_entra = H_mat1
#                          , C_entra = C_def, num_datos = dim(Base_ent_ts)[1])
# 
# 
# test_rest(fore_obj = fore_obj, fore_libre = fore_mod2$mean
#           , sigma_modelo = mod2$sigma2
#           , num_restri = num_restri 
#           , H_entra = H_mat2
#           , C_entra = C_def, num_datos = dim(Base_ent_ts)[1])



# -----------------------------------------------------------------------------#
# FIN DE PROGRAMA ----
# -----------------------------------------------------------------------------#
