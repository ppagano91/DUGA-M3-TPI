# ---- Estadística descriptiva de Precipitaciones ----
library(readxl)
library(ggplot2)
library(akima)
library(RColorBrewer)
library(geoR)
# install.packages("scatterplot3d")
library(scatterplot3d)
# install.packages("rgl")
library(rgl)
# install.packages("plotly")
library(plotly)

library(car)
library(fBasics)

library(moments)

# 1) Cargar datos
df <- read_excel("data/Precipitaciones.xlsx", sheet = 1)
View(df)
summary(df)
names(df) <- trimws(names(df))

# 2) Limpiar NAs en la variable de interés
df <- df[!is.na(df$Precipitaciones), ]

# 3) Estadísticos básicos
prec <- df$Precipitaciones

n <- length(prec)
min_val <- min(prec)
max_val <- max(prec)
mean_val <- mean(prec)
median_val <- median(prec)
sd_val <- sd(prec)
var_val <- var(prec)
range_val <- max_val - min_val

cat("Cantidad de registros:", n, "\n")
cat("Mínimo:", round(min_val,2), "\n")
cat("Máximo:", round(max_val,2), "\n")
cat("Promedio:", round(mean_val,2), "\n")
cat("Mediana:", round(median_val,2), "\n")
cat("Desvío estándar:", round(sd_val,2), "\n")
cat("Varianza:", round(var_val,2), "\n")
cat("Rango:", round(range_val,2), "\n")

# 4) Histogramas y boxplot

# p1 <- ggplot(df, aes(x = Precipitaciones)) +
#   geom_histogram(bins = 10, fill = "skyblue", color = "black") +
#   geom_vline(xintercept = mean_val, color = "red", linetype = "dashed") +
#   geom_vline(xintercept = median_val, color = "blue", linetype = "dotted") +
#   labs(title = "Histograma de Precipitaciones",
#        x = "Precipitación (mm/año)", y = "Frecuencia")
# 
# p2 <- ggplot(df, aes(y = Precipitaciones)) +
#   geom_boxplot(fill = "lightgreen") +
#   labs(title = "Diagrama de caja de Precipitaciones", y = "Precipitación (mm/año)")

par(mfrow = c(1,2))
hist(df$Precipitaciones,
     probability = TRUE,
     main = "Histograma",
     xlab = "Precipitación (mm/año)",
     ylab = "Frecuencia",
     col = "skyblue", border = "black")
lines(density(df$Precipitaciones, na.rm = TRUE),
      col = "red", lwd = 2)

boxplot(df$Precipitaciones,
        main = "Diagrama de Caja",
        ylab = "Precipitación (mm/año)",
        col = "lightgreen")
par(mfrow = c(1,1))

interp_ppt <- interp(x = df$Longitud,
                     y = df$Latitud,
                     z = df$Precipitaciones,
                     duplicate = "mean",
                     nx = 100, ny = 100)  # resolución

filled.contour(interp_ppt,
               color.palette = colorRampPalette(brewer.pal(10, "YlGnBu")),
               xlab = "Longitud",
               ylab = "Latitud",
               main = "Distribución espacial preliminar de Precipitaciones (mm/año)",
               plot.axes = {
                 axis(1); axis(2)
                 contour(interp_ppt, add = TRUE, drawlabels = FALSE, lwd = 0.5)
                 points(df$Longitud, df$Latitud, pch = 20, col = "black")
               },
               asp = 1)


df_geodata <- as.geodata(df, coords.col = c("Longitud", "Latitud"), data.col = "Precipitaciones")

duplicated(df_geodata)
dup.coords(df_geodata)

plot(df_geodata)

plot(df_geodata, trend = "1st")


# ---- Exploración de tendencias espaciales ----

# Scatter Longitud vs Precipitaciones
plot(df$Longitud, df$Precipitaciones,
     pch = 20, col = "blue",
     main = "Longitud vs Precipitaciones",
     xlab = "Longitud", ylab = "Precipitación (mm/año)")
abline(lm(Precipitaciones ~ Longitud, data = df), col = "red", lwd = 2)

# Scatter Latitud vs Precipitaciones
plot(df$Latitud, df$Precipitaciones,
     pch = 20, col = "darkgreen",
     main = "Latitud vs Precipitaciones",
     xlab = "Latitud", ylab = "Precipitación (mm/año)")
abline(lm(Precipitaciones ~ Latitud, data = df), col = "red", lwd = 2)

# Modelos lineales de tendencia
m_lon <- lm(Precipitaciones ~ Longitud, data = df)
m_lat <- lm(Precipitaciones ~ Latitud, data = df)

summary(m_lon)
summary(m_lat)

# ---- Superficie de tendencia (modelo espacial global) ----
# Ajuste polinómico 1er orden: Precip = a + bX + cY
m_trend <- lm(Precipitaciones ~ Longitud + Latitud, data = df)
summary(m_trend)

# Crear grilla para predecir
lon_seq <- seq(min(df$Longitud), max(df$Longitud), length.out = 100)
lat_seq <- seq(min(df$Latitud), max(df$Latitud), length.out = 100)
grid <- expand.grid(Longitud = lon_seq, Latitud = lat_seq)
grid$pred <- predict(m_trend, newdata = grid)

# Mapa de la superficie de tendencia
z <- matrix(grid$pred, nrow = 100, ncol = 100)
filled.contour(lon_seq, lat_seq, z,
               color.palette = colorRampPalette(brewer.pal(9, "YlGnBu")),
               xlab = "Longitud", ylab = "Latitud",
               main = "Superficie de tendencia (modelo lineal)",
               plot.axes = {
                 axis(1); axis(2)
                 points(df$Longitud, df$Latitud, pch = 20, col = "black")
               },
               asp = 1)

# La superficie de tendencia de primer orden evidencia un gradiente espacial en sentido oeste–este, con valores crecientes de precipitación hacia el este. Este patrón sugiere la presencia de una tendencia lineal global, de intensidad moderada, que deberá considerarse al momento de ajustar el modelo geoestadístico.


filled.contour(lon_seq, lat_seq, z,
               color.palette = colorRampPalette(brewer.pal(9, "YlGnBu")),
               xlab = "Longitud", ylab = "Latitud",
               main = "Superficie de tendencia (modelo lineal)",
               plot.axes = {
                 axis(1); axis(2)
                 points(df$Longitud, df$Latitud, pch = 20, col = "black")
                 text(df$Longitud, df$Latitud,
                      labels = round(df$Precipitaciones,0), cex=0.7, pos=3)
               },
               asp = 1)


scatterplot3d(df$Longitud, df$Latitud, df$Precipitaciones,
              pch = 20, color = "blue", main = "Plano de tendencia")
s3d <- scatterplot3d(df$Longitud, df$Latitud, df$Precipitaciones)
s3d$plane3d(m_trend, draw_polygon = TRUE, polygon_args = list(col=rgb(0.2,0.5,1,0.3)))



# Scatter de puntos 3D
plot3d(df$Longitud, df$Latitud, df$Precipitaciones,
       col = "blue", size = 5,
       xlab = "Longitud", ylab = "Latitud", zlab = "Precipitaciones")

# Plano de tendencia
grid$pred <- predict(m_trend, newdata = grid)
z <- matrix(grid$pred, nrow = length(lon_seq), ncol = length(lat_seq))
surface3d(lon_seq, lat_seq, z, color = "lightblue", alpha = 0.5)



plot_ly(df, x = ~Longitud, y = ~Latitud, z = ~Precipitaciones,
        type = "scatter3d", mode = "markers",
        marker = list(size = 4, color = "blue")) %>%
  add_surface(x = lon_seq, y = lat_seq, z = matrix(grid$pred, 100, 100),
              colorscale = "Blues", opacity = 0.5) %>%
  layout(scene = list(
    xaxis = list(title = "Longitud"),
    yaxis = list(title = "Latitud"),
    zaxis = list(title = "Precipitaciones")
  ))




# ---- Análisis de normalidad ----
# Tests estadísticos
shapiro <- shapiro.test(df$Precipitaciones)  # Shapiro-Wilk
cat("Shapiro-Wilk W =", round(shapiro$statistic,3),
    "p =", signif(shapiro$p.value,3), "\n")

skewness(df$Precipitaciones)   # asimetría
kurtosis(df$Precipitaciones)  

# Transformaciones
df$prec_log  <- ifelse(df$Precipitaciones > 0, log(df$Precipitaciones), NA)
df$prec_sqrt <- sqrt(df$Precipitaciones)

sh_log  <- shapiro.test(na.omit(df$prec_log))
sh_sqrt <- shapiro.test(df$prec_sqrt)

cat("Shapiro log  p =", signif(sh_log$p.value,3), "\n")
cat("Shapiro sqrt p =", signif(sh_sqrt$p.value,3), "\n")

# ---- QQ-plots ----
par(mfrow=c(1,3))
qqnorm(df$Precipitaciones, main="QQ-Plot Original")
qqline(df$Precipitaciones, col="red")

qqnorm(df$prec_log, main="QQ-Plot Log", ylab="log(precipitaciones)")
qqline(df$prec_log, col="red")

qqnorm(df$prec_sqrt, main="QQ-Plot Sqrt", ylab="sqrt(precipitaciones)")
qqline(df$prec_sqrt, col="red")
par(mfrow=c(1,1))


hist(df$Precipitaciones, prob = TRUE,
     main = "Histograma con curva normal",
     ylab = "Frecuencia",
     xlab = "Precipitación (mm/año)")
x <- seq(min(df$Precipitaciones), max(df$Precipitaciones), length = 40)
f <- dnorm(x, mean = mean(df$Precipitaciones), sd = sd(df$Precipitaciones))
lines(x, f, col = "red", lwd = 2)


shapiro <- shapiro.test(df$prec_log)  # Shapiro-Wilk
cat("Shapiro-Wilk W =", round(shapiro$statistic,3),
    "p =", signif(shapiro$p.value,3), "\n")
skewness(df$prec_log)   # asimetría
kurtosis(df$prec_log)   # curtosis

shapiro <- shapiro.test(df$prec_sqrt)  # Shapiro-Wilk
cat("Shapiro-Wilk W =", round(shapiro$statistic,3),
    "p =", signif(shapiro$p.value,3), "\n")
skewness(df$prec_sqrt)   # asimetría
kurtosis(df$prec_sqrt)   # curtosis


bc_trans <- powerTransform(df$Precipitaciones)
summary(bc_trans)

# aplicar la transformación sugerida
df$prec_bc <- bcPower(df$Precipitaciones, coef(bc_trans))
shapiro.test(df$prec_bc)

hist(df$prec_bc, prob = TRUE,
     main = "Histograma con curva normal",
     ylab = "Frecuencia",
     xlab = "Precipitación (mm/año)")
x <- seq(min(df$prec_bc), max(df$prec_bc), length = 40)
f <- dnorm(x, mean = mean(df$prec_bc), sd = sd(df$prec_bc))
lines(x, f, col = "red", lwd = 2)




# ---- ANÁLISIS GEOESTADÍSTICO ----
# Requisitos: geoR, gstat, sp, sf, raster
# install.packages(c("geoR","gstat","sp","sf","raster","spdep")) # si falta alguno

library(geoR)
library(gstat)
library(sp)
library(sf)
library(raster)

# 1) Preparar datos: si detectaste tendencia fuerte, modelala y usa residuos
# (si no, comentá la parte de la regresión y usá Precipitaciones directamente)
m_trend <- lm(Precipitaciones ~ Longitud + Latitud, data = df)
summary(m_trend)
# residuos
df$resid_trend <- residuals(m_trend)

# Opción: usar residuos para variograma (remueve tendencia global)
use_residuals_for_variogram <- FALSE

if(use_residuals_for_variogram){
  data_for_geo <- df
  data_for_geo$zvar <- df$resid_trend
  message("Usando residuos del modelo de tendencia para el variograma.")
} else {
  data_for_geo <- df
  data_for_geo$zvar <- df$Precipitaciones
  message("Usando Precipitaciones originales para el variograma.")
}

# 2) Crear objeto geoR
geod <- as.geodata(data_for_geo, coords.col = c("Longitud","Latitud"), data.col = "zvar")

# 3) Verificar duplicados
dup.coords(geod)   # imprime info si hay duplicados

# 4) Variograma experimental
maxd <- max(dist(geod$coords))/2  # max dist considered
uvec <- seq(0, maxd, length.out = 12)  # lag classes
vario_exp <- variog(geod, uvec = uvec, tol.hor = uvec[2]/2)
plot(vario_exp, main = "Variograma experimental (residuos si aplicable)")

# 5) Ajuste de modelos teóricos con likfit (ML y REML)
# Modelos a probar: spherical (sph), exponential (exp), gaussian (gau). 'wave' opcional.
models <- c("sph","exp","gau")
emp_var <- var(geod$data, na.rm = TRUE)
ini_cov_pars <- c(emp_var*0.6, 0.25 * max(dist(geod$coords)))  # sigma2, range: heurística

fits <- list()
for(mod in models){
  fits[[mod]] <- list()
  # tratar de ajustar REML y ML
  fits[[mod]]$ml <- tryCatch(likfit(geod, cov.model = mod, ini.cov.pars = ini_cov_pars,
                                    lik.method = "ML", nugget = emp_var*0.05),
                             error=function(e) NULL)
  fits[[mod]]$reml <- tryCatch(likfit(geod, cov.model = mod, ini.cov.pars = ini_cov_pars,
                                      lik.method = "REML", nugget = emp_var*0.05),
                               error=function(e) NULL)
}

# 6) Mostrar resultados resumidos y comparar AIC
aic_tab <- data.frame(model=character(), method=character(), AIC=numeric(), stringsAsFactors=FALSE)
for(mod in names(fits)){
  if(!is.null(fits[[mod]]$ml))  aic_tab <- rbind(aic_tab, data.frame(model=mod, method="ML", AIC=fits[[mod]]$ml$AIC))
  if(!is.null(fits[[mod]]$reml)) aic_tab <- rbind(aic_tab, data.frame(model=mod, method="REML", AIC=fits[[mod]]$reml$AIC))
}
print(aic_tab)
best_row <- aic_tab[which.min(aic_tab$AIC), ]
cat("Mejor por AIC:", best_row$model, best_row$method, "\n")
best_fit <- fits[[as.character(best_row$model)]][[tolower(as.character(best_row$method))]]

# 7) Graficar variograma experimental con la(s) curva(s) ajustadas
plot(vario_exp, main = "Variograma experimental con ajustes")
cols <- c("red","blue","green")
i <- 1
for(mod in names(fits)){
  if(!is.null(fits[[mod]]$reml)){
    lines(fits[[mod]]$reml, col = cols[i], lwd = 2)
    i <- i + 1
  }
}
legend("topright", legend = names(fits), col = cols[1:length(names(fits))], lty = 1, bty = "n")

# 8) Interpretación directa de parámetros (ejemplo)
# best_fit$cov.pars -> c(sigma2, range)
# best_fit$nugget -> nugget
if(!is.null(best_fit)){
  sigma2 <- best_fit$cov.pars[1]
  phi <- best_fit$cov.pars[2]
  nug <- best_fit$nugget
  cat(sprintf("Parámetros del mejor ajuste (%s): pepita (nugget)=%.4f, sigma2(partial sill)=%.4f, rango(phi)=%.4f\n",
              as.character(best_row$model), nug, sigma2, phi))
  cat("Interpretación: pepita = variación a distancia casi 0 (errores+microvariabilidad).\n")
  cat("Meseta (sill total) = nugget + sigma2. Rango = distancia a la que la autocorrelación se atenúa.\n")
}

# 9) Validación por leave-one-out con gstat
# Convertir a SpatialPointsDataFrame
library(sp)

coords_sp <- data.frame(Longitud = data_for_geo$Longitud, Latitud = data_for_geo$Latitud)
spdf <- SpatialPointsDataFrame(coords = coords_sp, data = data.frame(zvar = data_for_geo$zvar))
proj4string(spdf) <- CRS("+init=epsg:4326")

# Convertir best_fit a modelo gstat (si best_fit existe)
library(gstat)
library(sf)
if(!is.null(best_fit)){
  # geoR: cov.pars = c(sigma2, phi); nugget
  psill <- best_fit$cov.pars[1]   # partial sill
  range_g <- best_fit$cov.pars[2]
  nugget_g <- best_fit$nugget
  # vgm_gstat <- vgm(psill, model = as.character(best_row$model), range = range_g, nugget = nugget_g)
  vgm_gstat <- vgm(psill, model = "Sph", range = range_g, nugget = nugget_g)
  # formula (if used residuals we krige residuos; predictions then add tendencia de vuelta)
  g <- gstat(formula = zvar ~ 1, data = spdf, model = vgm_gstat)
  cv <- krige.cv(zvar ~ 1, spdf, model = vgm_gstat, nmax = 30)
  rmse <- sqrt(mean((cv$observed - cv$var1.pred)^2, na.rm=TRUE))
  bias <- mean(cv$observed - cv$var1.pred, na.rm=TRUE)
  cat(sprintf("LOOCV: RMSE = %.4f | Bias = %.4f\n", rmse, bias))
}

# 10) Kriging: generar grilla de predicción (en lon/lat; idealmente proyectar a UTM antes)
lon_min <- min(df$Longitud); lon_max <- max(df$Longitud)
lat_min <- min(df$Latitud); lat_max <- max(df$Latitud)
nx <- 150; ny <- 150
xg <- seq(lon_min, lon_max, length.out = nx)
yg <- seq(lat_min, lat_max, length.out = ny)
pred_grid <- expand.grid(Longitud = xg, Latitud = yg)
coordinates(pred_grid) <- ~ Longitud + Latitud
proj4string(pred_grid) <- CRS("+init=epsg:4326")

# Kriging ordinario con gstat (krige)
if(!is.null(best_fit)){
  kriged <- predict(g, pred_grid)
  # Si usamos residuos: sumamos la tendencia de vuelta
  if(use_residuals_for_variogram){
    # predecir tendencia en la grilla y sumarla
    trend_pred <- predict(m_trend, newdata = as.data.frame(pred_grid))
    kriged@data$var1.pred_trend <- kriged@data$var1.pred + trend_pred
    kriged@data$var1.var <- kriged@data$var1.var
    # Para exportar, usaremos var1.pred_trend como predicción final
    kriged@data$final_pred <- kriged@data$var1.pred_trend
  } else {
    kriged@data$final_pred <- kriged@data$var1.pred
  }
}

library(raster)
# 11) Convertir a raster y graficar predicción + varianza
if(exists("kriged")){
  pred_spdf <- SpatialPixelsDataFrame(points = coordinates(pred_grid), data = data.frame(kriged@data))
  r_pred <- raster(pred_spdf, layer = "final_pred")
  r_var  <- raster(pred_spdf, layer = "var1.var")
  # Plot sencillo
  plot(r_pred, main = "Kriging - Predicción de Precipitaciones (mm/año)")
  points(df$Longitud, df$Latitud, pch = 20)
  plot(r_var, main = "Kriging - Varianza (incertidumbre)")
  points(df$Longitud, df$Latitud, pch = 20)
  # Guardar GeoTIFF
  writeRaster(r_pred, "predicciones_precipitaciones_kriging.tif", format="GTiff", overwrite=TRUE)
  writeRaster(r_var, "varianzas_precipitaciones_kriging.tif", format="GTiff", overwrite=TRUE)
  cat("Rasters exportados: predicciones_precipitaciones_kriging.tif, varianzas_precipitaciones_kriging.tif\n")
}

# 12) Exportar contornos y puntos de muestreo si querés (similar a lo que tenías)
# ... (usar contourLines y st_write como en tus ejemplos anteriores)
