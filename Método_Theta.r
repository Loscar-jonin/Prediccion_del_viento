# ============================================================
#  MODELO THETA – PREDICCIÓN DE VELOCIDAD Y DIRECCIÓN DEL VIENTO
#  Frecuencia: Mensual  |  Horizonte de predicción: 12 meses
#
#  El Método Theta descompone la serie en dos "líneas theta":
#    θ=0 (tendencia lineal) y θ=2 (amplifica la curvatura local).
#  El pronóstico es el promedio de ambas + suavizamiento exponencial.
#  Ganó el concurso M3 de predicción de series temporales (2000).
#
#  Prerrequisito: ejecutar datos_originales.R,
#                 Interpolacion_vel_viento.R e
#                 Interpolacion_dir_viento.R
# ============================================================

# ============================================================
# ── 1. RECONSTRUIR SERIES INTERPOLADAS ──────────────────────
# ============================================================

inicio_vv <- c(year(min(vv_completo$Fecha)), month(min(vv_completo$Fecha)))
inicio_dv <- c(year(min(dv_completo$Fecha)), month(min(dv_completo$Fecha)))

ts_vel <- na.interp(ts(vv_completo$Valor, start = inicio_vv, frequency = 12))
ts_dir <- na.interp(ts(dv_completo$Valor, start = inicio_dv, frequency = 12))

cat("Serie velocidad – longitud:", length(ts_vel), "meses\n")
cat("Serie dirección – longitud:", length(ts_dir), "meses\n")

# ============================================================
# ── 2. PARÁMETROS GLOBALES ──────────────────────────────────
# ============================================================

H <- 12   # horizonte de predicción (meses)

# ============================================================
# ── 3. FUNCIÓN PRINCIPAL: AJUSTAR Y EVALUAR THETA ───────────
# ============================================================

ajustar_theta <- function(serie, etiqueta, unidad_y) {
  
  n        <- length(serie)
  n_train  <- n - H
  ts_train <- window(serie, end   = time(serie)[n_train])
  ts_test  <- window(serie, start = time(serie)[n_train + 1])
  
  cat("\n══════════════════════════════════════════════\n")
  cat(" THETA –", etiqueta, "\n")
  cat("══════════════════════════════════════════════\n")
  cat(" Observaciones totales :", n, "\n")
  cat(" Entrenamiento         :", n_train, "meses\n")
  cat(" Prueba                :", H, "meses\n")
  
  # ── 3a. AJUSTE: thetaf ──────────────────────────────────────
  # thetaf selecciona automáticamente si aplica ajuste estacional
  # antes de aplicar el método Theta estándar.
  # level: intervalos de confianza al 80 % y 95 %
  pron <- thetaf(
    y     = ts_train,
    h     = H,
    level = c(80, 95)
  )
  
  # thetaf() no expone $par ni $drift directamente en el objeto forecast.
  # Se extraen explícitamente:
  #   alpha → ajustando SES sobre la serie de entrenamiento
  #   drift → pendiente de la regresión lineal (línea theta=0)
  ses_interno <- ses(ts_train, h = 1)
  alpha_val   <- round(ses_interno$model$par["alpha"], 5)
  
  t_vec_temp  <- seq_len(length(ts_train))
  lm_temp     <- lm(as.numeric(ts_train) ~ t_vec_temp)
  drift_val   <- round(coef(lm_temp)[2], 5)
  
  cat("\nParámetro alpha (SES) :", alpha_val, "\n")
  cat("Drift (tendencia lineal):", drift_val, "\n")
  
  real     <- as.numeric(ts_test)
  predicho <- as.numeric(pron$mean)
  
  # ── 3b. MÉTRICAS ────────────────────────────────────────────
  metricas <- data.frame(
    Metrica = c("RMSE", "MSE", "MAE", "MAPE (%)", "SMAPE (%)"),
    Theta   = c(
      Metrics::rmse(real, predicho),
      Metrics::mse(real,  predicho),
      Metrics::mae(real,  predicho),
      Metrics::mape(real, predicho) * 100,
      smape_fn(real, predicho)
    )
  ) %>% mutate(across(where(is.numeric), ~round(.x, 5)))
  
  cat("\n── Métricas de predicción ──────────────────\n")
  print(metricas)
  
  # ── 3c. GRÁFICA 1: Serie histórica + pronóstico ─────────────
  df_hist <- data.frame(
    Fecha = as.Date(as.yearmon(time(serie))),
    Valor = as.numeric(serie)
  )
  
  df_pred <- data.frame(
    Fecha  = as.Date(as.yearmon(time(pron$mean))),
    Valor  = as.numeric(pron$mean),
    Lo80   = as.numeric(pron$lower[, 1]),
    Hi80   = as.numeric(pron$upper[, 1]),
    Lo95   = as.numeric(pron$lower[, 2]),
    Hi95   = as.numeric(pron$upper[, 2])
  )
  
  p1 <- ggplot() +
    geom_line(data = df_hist,
              aes(x = Fecha, y = Valor), colour = "#1d3557", linewidth = 0.7) +
    geom_ribbon(data = df_pred,
                aes(x = Fecha, ymin = Lo95, ymax = Hi95, fill = "IC 95 %"),
                alpha = 0.4) +
    geom_ribbon(data = df_pred,
                aes(x = Fecha, ymin = Lo80, ymax = Hi80, fill = "IC 80 %"),
                alpha = 0.4) +
    scale_fill_manual(
      name   = "Intervalo de confianza",
      values = c("IC 95 %" = "#a8dadc",
                 "IC 80 %" = "#457b9d")
    ) +
    geom_line(data = df_pred,
              aes(x = Fecha, y = Valor),
              colour = "#e63946", linewidth = 0.9, linetype = "dashed") +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(
      title    = paste("Modelo Theta –", etiqueta),
      subtitle = paste("alpha =", alpha_val,
                       "| IC: 80 % y 95 %"),
      x = "Fecha", y = unidad_y
    ) +guides(
      colour = guide_legend(order = 1, title = NULL),
      fill   = guide_legend(order = 2)
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))
  
  print(p1)
  
  # ── 3d. GRÁFICA 2: Real vs Predicho ─────────────────────────
  df_comp <- data.frame(
    Fecha    = df_pred$Fecha,
    Real     = real,
    Predicho = predicho
  ) %>%
    pivot_longer(cols = c(Real, Predicho),
                 names_to = "Serie", values_to = "Valor")
  
  p2 <- ggplot(df_comp, aes(x = Fecha, y = Valor, colour = Serie)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    scale_colour_manual(values = c("Real" = "#1d3557", "Predicho" = "#e63946")) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
    labs(
      title    = paste("Real vs. Predicho – Theta –", etiqueta),
      subtitle = paste("RMSE =", round(Metrics::rmse(real, predicho), 4),
                       "| MAE =", round(Metrics::mae(real, predicho), 4)),
      x = "Fecha", y = unidad_y, colour = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))
  
  print(p2)
  
  # ── 3e. GRÁFICA 3: Dispersión Real vs Predicho ──────────────
  df_scatter <- data.frame(Real = real, Predicho = predicho)
  
  p3 <- ggplot(df_scatter, aes(x = Real, y = Predicho)) +
    geom_point(colour = "#e63946", size = 3, alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "#1d3557", linewidth = 0.8) +
    labs(
      title    = paste("Diagrama de dispersión – Theta –", etiqueta),
      subtitle = "Línea punteada = predicción perfecta (y = x)",
      x = paste("Real –", unidad_y),
      y = paste("Predicho –", unidad_y)
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
  
  print(p3)
  
  # ── 3f. VISUALIZACIÓN INTERNA: líneas theta (θ=0 y θ=2) ─────
  # Líneas theta para la descomposición visual
  # (alpha_val y drift_val ya fueron calculados arriba)
  n_tr  <- length(ts_train)
  t_vec <- seq_len(n_tr)
  
  # Línea theta=0: regresión lineal sobre la serie
  theta0 <- lm(as.numeric(ts_train) ~ t_vec)$fitted.values
  
  # Línea theta=2: 2*serie - theta0 (amplifica curvatura local)
  theta2 <- 2 * as.numeric(ts_train) - theta0
  
  df_theta <- data.frame(
    t      = as.Date(as.yearmon(time(ts_train))),
    Serie  = as.numeric(ts_train),
    Theta0 = theta0,
    Theta2 = theta2
  ) %>%
    pivot_longer(cols = c(Serie, Theta0, Theta2),
                 names_to = "Linea", values_to = "Valor")
  
  p4 <- ggplot(df_theta, aes(x = t, y = Valor, colour = Linea)) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(
      values = c("Serie"  = "#1d3557",
                 "Theta0" = "#2a9d8f",
                 "Theta2" = "#e9c46a"),
      labels = c("Serie"  = "Serie original",
                 "Theta0" = "θ = 0 (tendencia lineal)",
                 "Theta2" = "θ = 2 (curvatura amplificada)")
    ) +
    labs(
      title    = paste("Descomposición Theta –", etiqueta),
      subtitle = "θ=0 captura la tendencia; θ=2 la dinámica local",
      x = "Fecha", y = unidad_y, colour = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))
  
  print(p4)
  
  invisible(list(forecast = pron, metricas = metricas))
}

# ============================================================
# ── 4. APLICAR THETA A AMBAS VARIABLES ──────────────────────
# ============================================================

res_theta_vel <- ajustar_theta(
  serie    = ts_vel,
  etiqueta = "Velocidad del Viento",
  unidad_y = "Vel. Viento (m/s)"
)

res_theta_dir <- ajustar_theta(
  serie    = ts_dir,
  etiqueta = "Dirección del Viento",
  unidad_y = "Dir. Viento (sector)"
)

# ============================================================
# ── 5. TABLA COMPARATIVA FINAL ──────────────────────────────
# ============================================================

cat("\n══════════════════════════════════════════════\n")
cat("  RESUMEN COMPARATIVO – THETA\n")
cat("══════════════════════════════════════════════\n")

tabla_final <- data.frame(
  Metrica   = res_theta_vel$metricas$Metrica,
  Velocidad = res_theta_vel$metricas$Theta,
  Direccion = res_theta_dir$metricas$Theta
)

print(tabla_final)
