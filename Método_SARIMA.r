# ============================================================
#  MODELO SARIMA – PREDICCIÓN DE VELOCIDAD Y DIRECCIÓN DEL VIENTO
#  Frecuencia: Mensual  |  Horizonte de predicción: 12 meses
#  SARIMA = ARIMA con componente estacional (P,D,Q)[12]
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
# ── 3. FUNCIÓN PRINCIPAL: AJUSTAR Y EVALUAR SARIMA ──────────
# ============================================================

ajustar_sarima <- function(serie, etiqueta, unidad_y) {
  
  n        <- length(serie)
  n_train  <- n - H
  ts_train <- window(serie, end   = time(serie)[n_train])
  ts_test  <- window(serie, start = time(serie)[n_train + 1])
  
  cat("\n══════════════════════════════════════════════\n")
  cat(" SARIMA –", etiqueta, "\n")
  cat("══════════════════════════════════════════════\n")
  cat(" Observaciones totales :", n, "\n")
  cat(" Entrenamiento         :", n_train, "meses\n")
  cat(" Prueba                :", H, "meses\n")
  
  # ── 3a. AJUSTE: auto.arima con componente estacional (S=12) ─
  # D = diferenciación estacional, P/Q = AR/MA estacionales
  modelo <- auto.arima(
    ts_train,
    seasonal      = TRUE,    # SARIMA — componente estacional activo
    # D no se especifica: auto.arima lo detecta automáticamente
    # (pasar D = NULL explícitamente rompe su lógica interna)
    stepwise      = FALSE,   # búsqueda exhaustiva del espacio de órdenes
    approximation = FALSE,
    ic            = "aicc",
    trace         = FALSE
  )
  
  cat("\nModelo seleccionado:", as.character(modelo), "\n")
  cat("Parámetros:\n")
  print(coef(modelo))
  print(summary(modelo))
  
  # ── 3b. PRONÓSTICO ──────────────────────────────────────────
  pron <- forecast(modelo, h = H, level = c(80, 95))
  
  real     <- as.numeric(ts_test)
  predicho <- as.numeric(pron$mean)
  
  # ── 3c. MÉTRICAS ────────────────────────────────────────────
  metricas <- data.frame(
    Metrica = c("RMSE", "MSE", "MAE", "MAPE (%)", "SMAPE (%)"),
    SARIMA  = c(
      Metrics::rmse(real, predicho),
      Metrics::mse(real,  predicho),
      Metrics::mae(real,  predicho),
      Metrics::mape(real, predicho) * 100,
      smape_fn(real, predicho)
    )
  ) %>% mutate(across(where(is.numeric), ~round(.x, 5)))
  
  cat("\n── Métricas de predicción ──────────────────\n")
  print(metricas)
  
  # ── 3d. GRÁFICA 1: Serie histórica + pronóstico con IC ──────
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
              aes(x = Fecha, y = Valor),
              colour = "#e63946", linewidth = 0.9, linetype = "dashed") +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(
      title    = paste("SARIMA –", etiqueta),
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
    theme(plot.title = element_text(face = "bold"))
  
  print(p1)
  
  # ── 3e. GRÁFICA 2: Real vs Predicho ─────────────────────────
  meses_en <- c("Jan","Feb","Mar","Apr","May","Jun",
                "Jul","Aug","Sep","Oct","Nov","Dec")
  
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
      title    = paste("Real vs. Predicho – SARIMA –", etiqueta),
      subtitle = paste("RMSE =", round(Metrics::rmse(real, predicho), 4),
                       "| MAE =", round(Metrics::mae(real, predicho), 4)),
      x = "Fecha", y = unidad_y, colour = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))
  
  print(p2)
  
  # ── 3f. GRÁFICA 3: Descomposición STL para visualizar estacionalidad
  stl_desc <- stl(ts_train, s.window = "periodic")
  plot(stl_desc,
       main = paste("Descomposición STL –", etiqueta))
  
  # ── 3g. DIAGNÓSTICO DE RESIDUOS ─────────────────────────────
  checkresiduals(modelo)
  
  invisible(list(modelo = modelo, forecast = pron, metricas = metricas))
}

# ============================================================
# ── 4. APLICAR SARIMA A AMBAS VARIABLES ─────────────────────
# ============================================================

res_sarima_vel <- ajustar_sarima(
  serie    = ts_vel,
  etiqueta = "Velocidad del Viento",
  unidad_y = "Vel. Viento (m/s)"
)

res_sarima_dir <- ajustar_sarima(
  serie    = ts_dir,
  etiqueta = "Dirección del Viento",
  unidad_y = "Dir. Viento (sector)"
)

# ============================================================
# ── 5. TABLA COMPARATIVA FINAL ──────────────────────────────
# ============================================================

cat("\n══════════════════════════════════════════════\n")
cat("  RESUMEN COMPARATIVO – SARIMA\n")
cat("══════════════════════════════════════════════\n")

tabla_final <- data.frame(
  Metrica   = res_sarima_vel$metricas$Metrica,
  Velocidad = res_sarima_vel$metricas$SARIMA,
  Direccion = res_sarima_dir$metricas$SARIMA
)

print(tabla_final)
