# instancia_integradora_kriging.R
# Requisitos: instalar paquetes si no los tenés: install.packages(c("readxl","dplyr","sf","sp","geoR","gstat","raster","spatialEco","ggplot2","gridExtra"))
library(readxl)
library(dplyr)
library(geoR)
library(sp)
library(sf)
library(gstat)
library(raster)
library(ggplot2)
library(gridExtra)

# ---- 1) Cargar datos ----
df <- read_excel("data/Precipitaciones.xlsx", sheet = 1)
# Normalizar nombres de columnas (por si acaso)
names(df) <- trimws(names(df))

# Comprueba que estén las columnas esperadas:
# "Longitud", "Latitud", "Precipitaciones"
if(!all(c("Longitud","Latitud","Precipitaciones") %in% names(df))){
  stop("El archivo debe contener las columnas: Longitud, Latitud, Precipitaciones")
}

# Elimina filas con NA en coordenadas o precipitaciones
df <- df %>% filter(!is.na(Longitud) & !is.na(Latitud) & !is.na(Precipitaciones))

# ---- 2) Exploración rápida ----
summary(df$Precipitaciones)
hist(df$Precipitaciones, main="Histograma de Precipitaciones", xlab="Precipitaciones", breaks=10)
qqnorm(df$Precipitaciones); qqline(df$Precipitaciones)

# Test de normalidad (Shapiro-Wilk) - atención: para n>500 puede dejar de ser apropiado
shapiro_test <- tryCatch(shapiro.test(df$Precipitaciones), error = function(e) NULL)
print(shapiro_test)

# Si no es normal, probar transformaciones (log, sqrt)
# Guardamos transformaciones para comparar
df <- df %>%
  mutate(prec_log = ifelse(Precipitaciones > 0, log(Precipitaciones), NA),
         prec_sqrt = sqrt(ifelse(Precipitaciones >= 0, Precipitaciones, NA)))

# Plot comparativo
p1 <- ggplot(df, aes(Precipitaciones)) + geom_histogram(bins=15) + ggtitle("Original")
p2 <- ggplot(df, aes(prec_log)) + geom_histogram(bins=15) + ggtitle("Log")
p3 <- ggplot(df, aes(prec_sqrt)) + geom_histogram(bins=15) + ggtitle("Sqrt")
grid.arrange(p1,p2,p3, ncol=3)

# Elegir variable a modelar
# Si log mejora normalidad, usar log; si no, usar original.
use_log <- FALSE
if(!is.null(shapiro_test)){
  if(shapiro_test$p.value < 0.05){
    # no normal -> probar log p-value (si no-NA)
    sh_log <- tryCatch(shapiro.test(na.omit(df$prec_log)), error=function(e) NULL)
    if(!is.null(sh_log) && sh_log$p.value > shapiro_test$p.value){
      use_log <- TRUE
      df$Y <- df$prec_log
      message("Se usará log(Precipitaciones) para el análisis.")
    } else {
      df$Y <- df$Precipitaciones
      message("Se usará Precipitaciones (sin transformar) para el análisis.")
    }
  } else {
    df$Y <- df$Precipitaciones
    message("Variable aproximadamente normal (shapiro p>=0.05). Usamos la original.")
  }
} else {
  df$Y <- df$Precipitaciones
  message("Shapiro no pudo ejecutarse. Usamos la original.")
}

# ---- 3) Convertir a geoR (para usar geoR como en clase) ----
# geoR espera coords.col = c(long, lat) y data.col la variable
geod <- as.geodata(df, coords.col = c("Longitud","Latitud"), data.col = "Y")

# Duplicados de coordenadas
dup.coords(geod)  # imprime si hay duplicados

# ---- 4) Variograma experimental ----
bins <- seq(0, max(dist(geod$coords))/2, length.out = 12) # 12 lag classes
variog_exp <- variog(geod, uvec = bins)
plot(variog_exp, main = "Variograma experimental (Precipitaciones)")

# Posible ajuste a ojo con eyefit (opcional)
# windows(); eyefit(variog_exp); dev.off()   # si usás RStudio en Windows, descomentar

# ---- 5) Ajuste de modelos: probar esférico, exponencial, gaussian y wave si querés ----
# Valores iniciales de cov.pars: sill_partial (sigma^2), range (phi)
# Tomar varianza empírica como punto de partida
emp_var <- var(geod$data, na.rm=TRUE)
init_cov_pars <- c(emp_var*0.6, 0.1 * max(dist(geod$coords))) # heurística

# Ajuste ML y REML con modelos comunes
mods <- list()
models_to_try <- c("sph","exp","gau","wave")  # 'sph' = spherical in geoR notation
for(mod in models_to_try){
  # intentar ML
  fit_ml <- tryCatch(
    likfit(geod, cov.model = mod, ini.cov.pars = init_cov_pars, lik.method = "ML", nugget = emp_var*0.1),
    error = function(e) NULL
  )
  # intentar REML
  fit_reml <- tryCatch(
    likfit(geod, cov.model = mod, ini.cov.pars = init_cov_pars, lik.method = "REML", nugget = emp_var*0.1),
    error = function(e) NULL
  )
  mods[[mod]] <- list(ml = fit_ml, reml = fit_reml)
}

# Mostrar resultados resumidos
for(mod in names(mods)){
  cat("----- Modelo:", mod, "-----\n")
  if(!is.null(mods[[mod]]$ml)){
    print(summary(mods[[mod]]$ml))
  } else cat("ML: fallo o no convergió\n")
  if(!is.null(mods[[mod]]$reml)){
    print(summary(mods[[mod]]$reml))
  } else cat("REML: fallo o no convergió\n")
}

# Graficar variograma experimental y líneas de los ajustes (si convergieron)
plot(variog_exp, main="Variograma experimental con ajustes")
cols <- c("red","blue","green","purple")
i <- 1
for(mod in names(mods)){
  if(!is.null(mods[[mod]]$ml)){
    lines(mods[[mod]]$ml, col = cols[i], lwd = 2)
    i <- i + 1
  }
}
legend("topright", legend = names(mods), col = cols[1:length(names(mods))], lty=1, bty="n")

# ---- 6) Seleccionar el mejor modelo según criterio: AIC / criterio de validación cruzada ----
# Usar AIC de los objetos likfit si disponibles
aic_tab <- data.frame(model = character(), method = character(), aic = numeric(), stringsAsFactors = FALSE)
for(mod in names(mods)){
  if(!is.null(mods[[mod]]$ml)){
    aic_tab <- rbind(aic_tab, data.frame(model=mod, method="ML", aic=mods[[mod]]$ml$AIC))
  }
  if(!is.null(mods[[mod]]$reml)){
    aic_tab <- rbind(aic_tab, data.frame(model=mod, method="REML", aic=mods[[mod]]$reml$AIC))
  }
}
print(aic_tab)
# Escoger mínimo AIC
best_row <- aic_tab[which.min(aic_tab$aic),]
cat("Mejor ajuste por AIC:", paste(best_row$model, best_row$method), "\n")

best_model_name <- as.character(best_row$model)
best_method <- as.character(best_row$method)
best_fit <- mods[[best_model_name]][[tolower(best_method)]]

# ---- 7) Validación cruzada (geoR o gstat) ----
# Usamos leave-one-out con gstat para comparar predicciones
# Primero convertir df a SpatialPointsDataFrame
coordinates(df) <- ~ Longitud + Latitud
proj4string(df) <- CRS("+init=epsg:4326")

# gstat object para kriging con modelo elegido
# Convertir parámetros de geoR a gstat: gstat needs sill, range, nugget in specific format
# GeoR's cov.pars = c(sigma2, phi). For gstat, model = vgm(psill, model, range, nugget)
sigma2 <- best_fit$cov.pars[1]
phi <- best_fit$cov.pars[2]
nugget_val <- best_fit$nugget

# En gstat psill = partial sill (sigma2), nugget separate, range = phi
# vgm_model <- vgm(psill = sigma2, model = best_model_name, range = phi, nugget = nugget_val)
vgm_model <- vgm(psill = sigma2, model = "Sph", range = phi, nugget = nugget_val)

# gstat object
g <- gstat(formula = Y ~ 1, data = df, model = vgm_model)

# LOOCV
cv <- gstat::krige.cv(Y ~ 1, df, model = vgm_model, nmax = 30)  # nmax neighbors
# Métricas
rmse <- sqrt(mean((cv$observed - cv$var1.pred)^2, na.rm=TRUE))
bias <- mean(cv$observed - cv$var1.pred, na.rm=TRUE)
cat(sprintf("LOOCV: RMSE = %.4f, Bias = %.4f\n", rmse, bias))

# ---- 8) Grilla de predicción ----
# Definir grilla con resolución (por ejemplo 100 x 100)
lon_min <- min(df$Longitud); lon_max <- max(df$Longitud)
lat_min <- min(df$Latitud); lat_max <- max(df$Latitud)
nx <- 150; ny <- 150
xg <- seq(lon_min, lon_max, length.out = nx)
yg <- seq(lat_min, lat_max, length.out = ny)
pred_grid <- expand.grid(Longitud = xg, Latitud = yg)
coordinates(pred_grid) <- ~ Longitud + Latitud
proj4string(pred_grid) <- CRS("+init=epsg:4326")

# Kriging ordinario con gstat
kriged <- predict(g, pred_grid)
# kriged tiene columns var1.pred (pred), var1.var (kriging var)

# ---- 9) Exportar raster de predicciones y varianzas ----
# Convertir a raster
# Primero convertir a SpatialPixelsDataFrame
pred_spdf <- SpatialPixelsDataFrame(points = coordinates(pred_grid), data = data.frame(kriged@data))
r_pred <- raster(pred_spdf, layer = "var1.pred")
r_var  <- raster(pred_spdf, layer = "var1.var")

# Guardar GeoTIFF
writeRaster(r_pred, filename = "predicciones_precipitaciones.tif", format="GTiff", overwrite=TRUE)
writeRaster(r_var, filename = "varianzas_precipitaciones.tif", format="GTiff", overwrite=TRUE)
cat("Rasteres escritos: predicciones_precipitaciones.tif, varianzas_precipitaciones.tif\n")

# ---- 10) Exportar contornos (opcional) y puntos de muestreo como shapefiles ----
library(reshape2)
# Contornos
zmat <- matrix(kriged@data$var1.pred, nrow = nx, ncol = ny)
conts <- contourLines(x = xg, y = yg, z = zmat)
# Convertir a sf líneas
lines_list <- lapply(conts, function(cl){
  st_linestring(cbind(cl$x, cl$y))
})
if(length(lines_list) > 0){
  sfc_lines <- st_sfc(lines_list, crs = 4326)
  levels <- sapply(conts, function(c) c$level)
  contours_sf <- st_sf(level = levels, geometry = sfc_lines)
  st_write(contours_sf, "contornos_prediccion.shp", delete_layer = TRUE)
  cat("Contornos escritos: contornos_prediccion.shp\n")
}

# Puntos de muestreo
puntos_sf <- st_as_sf(as.data.frame(df), coords = c("Longitud","Latitud"), crs = 4326)
st_write(puntos_sf, "puntos_muestreo.shp", delete_layer = TRUE)
cat("Puntos de muestreo escritos: puntos_muestreo.shp\n")

# ---- 11) Plots finales (mapa con puntos + raster) ----
plot(r_pred, main = "Predicciones - Precipitaciones")
plot(puntos_sf$geometry, add = TRUE, pch = 20)

# Mapa de varianza
plot(r_var, main = "Varianzas de Kriging")
plot(puntos_sf$geometry, add = TRUE, pch = 20)

# ---- 12) Guardar resultados del modelo ----
save(best_fit, file = "mejor_modelo_geoR.RData")
write.csv(aic_tab, "tabla_aic_modelos.csv", row.names = FALSE)
write.csv(as.data.frame(cv), "cv_loocv_predicciones.csv", row.names = FALSE)

cat("Proceso terminado. Archivos generados.\n")
