# ============================================================
#  MÉTODO ARIMA – PREDICCIÓN DE VELOCIDAD Y DIRECCIÓN DEL VIENTO
#  Frecuencia: Mensual  |  Horizonte de predicción: 12 meses
#  Prerrequisito: ejecutar datos_originales.R,
#                 Interpolacion_vel_viento.R e
#                 Interpolacion_dir_viento.R
# ============================================================

# ============================================================
# ── 1. RECONSTRUIR SERIES INTERPOLADAS ──────────────────────
# ============================================================

# Se usa la interpolación estacional (STL + Kalman) por ser la
# más robusta para datos con estacionalidad anual.
# Asegúrate de que vv_completo y dv_completo existen en el entorno.

inicio_vv <- c(year(min(vv_completo$Fecha)), month(min(vv_completo$Fecha)))
inicio_dv <- c(year(min(dv_completo$Fecha)), month(min(dv_completo$Fecha)))

ts_vel <- na.interp(ts(vv_completo$Valor, start = inicio_vv, frequency = 12))
ts_dir <- na.interp(ts(dv_completo$Valor, start = inicio_dv, frequency = 12))

cat("Serie velocidad  – longitud:", length(ts_vel), "meses\n")
cat("Serie dirección  – longitud:", length(ts_dir), "meses\n")

# ============================================================
# ── 2. PARÁMETROS GLOBALES DEL EXPERIMENTO ──────────────────
# ============================================================

H <- 12   # horizonte de predicción (meses)

# ============================================================
# ── 3. FUNCIÓN PRINCIPAL: AJUSTAR Y EVALUAR ARIMA ───────────
# ============================================================

ajustar_arima <- function(serie, etiqueta, unidad_y) {

  n        <- length(serie)
  n_train  <- n - H
  ts_train <- window(serie, end   = time(serie)[n_train])
  ts_test  <- window(serie, start = time(serie)[n_train + 1])

  cat("\n══════════════════════════════════════════════\n")
  cat(" ARIMA –", etiqueta, "\n")
  cat("══════════════════════════════════════════════\n")
  cat(" Observaciones totales :", n, "\n")
  cat(" Entrenamiento         :", n_train, "meses\n")
  cat(" Prueba                :", H, "meses\n")

  # ── 3a. AJUSTE: auto.arima sin componente estacional ────────
  # (seasonal = FALSE distingue ARIMA puro de SARIMA)
  modelo <- auto.arima(
    ts_train,
    seasonal      = FALSE,   # ARIMA no estacional
    stepwise      = FALSE,   # búsqueda exhaustiva
    approximation = FALSE,
    ic            = "aicc",
    trace         = FALSE
  )

  cat("\nModelo seleccionado:", as.character(modelo), "\n")
  print(summary(modelo))

  # ── 3b. PRONÓSTICO ──────────────────────────────────────────
  pron <- forecast(modelo, h = H, level = c(80, 95))

  real      <- as.numeric(ts_test)
  predicho  <- as.numeric(pron$mean)

  # ── 3c. MÉTRICAS ────────────────────────────────────────────
  metricas <- data.frame(
    Metrica = c("RMSE", "MSE", "MAE", "MAPE (%)", "SMAPE (%)"),
    ARIMA   = c(
      Metrics::rmse(real, predicho),
      Metrics::mse(real,  predicho),
      Metrics::mae(real,  predicho),
      Metrics::mape(real, predicho) * 100,
      smape_fn(real, predicho)
    )
  ) %>% mutate(across(where(is.numeric), ~round(.x, 5)))

  cat("\n── Métricas de predicción ──────────────────\n")
  print(metricas)

  # ── 3d. GRÁFICA 1: Serie completa + pronóstico ──────────────
  df_hist <- data.frame(
    Fecha = as.Date(as.yearmon(time(serie))),
    Valor = as.numeric(serie),
    Tipo  = "Histórico"
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
                 "IC 80 %" = "#457b9d")) +
    geom_line(data = df_pred,
              aes(x = Fecha, y = Valor), colour = "#e63946",
              linewidth = 0.9, linetype = "dashed") +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(
      title    = paste("ARIMA –", etiqueta),
      subtitle = paste("Modelo:", as.character(modelo),
                       "| IC: 80 % y 95 %"),
      x = "Fecha", y = unidad_y
    ) +
    guides(
      colour = guide_legend(order = 1, title = NULL),
      fill   = guide_legend(order = 2)
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "bottom",
      plot.title      = element_text(face = "bold"),
      axis.text.x     = element_text(angle = 45, hjust = 1)
    )
  
  print(p1)

  # ── 3e. GRÁFICA 2: Real vs Predicho (conjunto de prueba) ────
  df_comp <- data.frame(
    Fecha    = df_pred$Fecha,
    Real     = real,
    Predicho = predicho
  ) %>%
    pivot_longer(cols = c(Real, Predicho),
                 names_to = "Serie", values_to = "Valor")

  p2 <- ggplot(df_comp, aes(x = Fecha, y = Valor, colour = Serie)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    scale_colour_manual(values = c("Real" = "#1d3557", "Predicho" = "#e63946")) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
    labs(
      title    = paste("Real vs. Predicho – ARIMA –", etiqueta),
      subtitle = paste("RMSE =", round(Metrics::rmse(real, predicho), 4),
                       "| MAE =", round(Metrics::mae(real, predicho), 4)),
      x = "Fecha", y = unidad_y, colour = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))

  print(p2)

  # ── 3f. GRÁFICA 3: Diagnóstico de residuos ──────────────────
  checkresiduals(modelo)

  invisible(list(modelo = modelo, forecast = pron, metricas = metricas))
}

# ============================================================
# ── 4. APLICAR ARIMA A AMBAS VARIABLES ──────────────────────
# ============================================================

res_arima_vel <- ajustar_arima(
  serie    = ts_vel,
  etiqueta = "Velocidad del Viento",
  unidad_y = "Vel. Viento (m/s)"
)

res_arima_dir <- ajustar_arima(
  serie    = ts_dir,
  etiqueta = "Dirección del Viento",
  unidad_y = "Dir. Viento (sector)"
)

# ============================================================
# ── 5. TABLA COMPARATIVA FINAL ──────────────────────────────
# ============================================================

cat("\n══════════════════════════════════════════════\n")
cat("  RESUMEN COMPARATIVO – ARIMA\n")
cat("══════════════════════════════════════════════\n")

tabla_final <- data.frame(
  Metrica   = res_arima_vel$metricas$Metrica,
  Velocidad = res_arima_vel$metricas$ARIMA,
  Direccion = res_arima_dir$metricas$ARIMA
)

print(tabla_final)
