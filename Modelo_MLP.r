# ============================================================
#  MODELO MLP – PERCEPTRÓN MULTICAPA PARA SERIES DE TIEMPO
#  Predicción de Velocidad y Dirección del Viento
#  Frecuencia: Mensual  |  Horizonte de predicción: 12 meses
#
#  Usa el paquete {nnfor} que implementa MLP con:
#    – Selección automática de lags mediante validación cruzada
#    – Múltiples inicializaciones aleatorias (ensemble)
#    – Normalización automática de entradas/salidas
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

H       <- 12    # horizonte de predicción (meses)
REPS    <- 20    # redes en el ensemble (más = más estable, más lento)
set.seed(42)     # reproducibilidad

# ============================================================
# ── 3. FUNCIÓN PRINCIPAL: AJUSTAR Y EVALUAR MLP ─────────────
# ============================================================

ajustar_mlp <- function(serie, etiqueta, unidad_y) {
  
  n        <- length(serie)
  n_train  <- n - H
  ts_train <- window(serie, end   = time(serie)[n_train])
  ts_test  <- window(serie, start = time(serie)[n_train + 1])
  
  cat("\n══════════════════════════════════════════════\n")
  cat(" MLP –", etiqueta, "\n")
  cat("══════════════════════════════════════════════\n")
  cat(" Observaciones totales :", n, "\n")
  cat(" Entrenamiento         :", n_train, "meses\n")
  cat(" Prueba                :", H, "meses\n")
  cat(" Redes en el ensemble  :", REPS, "\n")
  
  # ── 3a. AJUSTE: mlp() ───────────────────────────────────────
  # hd        = NULL → nnfor elige automáticamente el nº de nodos ocultos
  # lags      = NULL → selección automática de retardos (lags)
  # reps      = REPS → ensemble de redes para reducir varianza
  # difforder = NULL → diferenciación automática si la serie no es estacionaria
  # sel.lag   = TRUE → selección estadística de lags relevantes
  
  modelo <- mlp(
    y         = ts_train,
    hd        = NULL,      # neuronas ocultas: auto
    lags      = NULL,      # retardos: auto
    reps      = REPS,
    sel.lag   = TRUE,
    difforder = NULL,
    allow.det.season = TRUE   # modela estacionalidad determinista
  )
  
  cat("\nResumen del modelo MLP:\n")
  print(modelo)
  
  # ── 3b. PRONÓSTICO ──────────────────────────────────────────
  pron <- forecast(modelo, h = H)
  
  real     <- as.numeric(ts_test)
  Predicho <- as.numeric(pron$mean)
  
  # ── 3c. MÉTRICAS ────────────────────────────────────────────
  metricas <- data.frame(
    Metrica = c("RMSE", "MSE", "MAE", "MAPE (%)", "SMAPE (%)"),
    MLP     = c(
      Metrics::rmse(real, Predicho),
      Metrics::mse(real,  Predicho),
      Metrics::mae(real,  Predicho),
      Metrics::mape(real, Predicho) * 100,
      smape_fn(real, Predicho)
    )
  ) %>% mutate(across(where(is.numeric), ~round(.x, 5)))
  
  cat("\n── Métricas de predicción ──────────────────\n")
  print(metricas)
  
  # ── 3d. GRÁFICA 1: Ajuste sobre entrenamiento + pronóstico ──
  df_hist <- data.frame(
    Fecha = as.Date(as.yearmon(time(serie))),
    Valor = as.numeric(serie)
  )
  
  # fitted(modelo) tiene menos filas que ts_train porque los primeros
  # p = lag_max meses no tienen historia suficiente para ser ajustados.
  # Se calcula el offset y se recortan las fechas de ts_train en consecuencia.
  fitted_vals <- as.numeric(fitted(modelo))
  offset_lags <- length(ts_train) - length(fitted_vals)
  fechas_fit  <- as.Date(as.yearmon(time(ts_train)))[(offset_lags + 1):length(ts_train)]
  
  df_fitted <- data.frame(
    Fecha  = fechas_fit,
    Ajuste = fitted_vals
  )
  
  df_pred <- data.frame(
    Fecha    = as.Date(as.yearmon(time(pron$mean))),
    Predicho = as.numeric(pron$mean)
  )
  
  p1 <- ggplot() +
    geom_line(data = df_hist,
              aes(x = Fecha, y = Valor, colour = "Historic"),
              linewidth = 0.7) +
    geom_line(data = df_fitted,
              aes(x = Fecha, y = Ajuste, colour = "Adjust (train)"),
              linewidth = 0.7, linetype = "dotted") +
    geom_line(data = df_pred,
              aes(x = Fecha, y = Predicho, colour = "MLP Pronostic"),
              linewidth = 1.0, linetype = "dashed") +
    geom_point(data = df_pred,
               aes(x = Fecha, y = Predicho, colour = "MLP Pronostic"),
               size = 2.5) +
    scale_colour_manual(
      values = c("Historic"     = "#1d3557",
                 "Adjust (train)"= "#2a9d8f",
                 "MLP Pronostic"    = "#e63946")
    ) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
labs(
      title    = paste("MLP –", etiqueta),
      subtitle = paste("Ensemble of", REPS, "self-selected networks | lags and neurons"),
      x = "Years", y = unidad_y, colour = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))
  
  print(p1)
  
  # ── 3e. GRÁFICA 2: Real vs Predicho en el período de prueba ─
  df_comp <- data.frame(
    Fecha    = df_pred$Fecha,
    Real     = real,
    Predicho = Predicho
  ) %>%
    pivot_longer(cols = c(Real, Predicho),
                 names_to = "Serie", values_to = "Valor")
  
  p2 <- ggplot(df_comp, aes(x = Fecha, y = Valor, colour = Serie)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    scale_colour_manual(values = c("Real" = "#1d3557", "Predicho" = "#e63946")) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
    labs(
      title    = paste("Real vs. Predicho – MLP –", etiqueta),
      subtitle = paste("RMSE =", round(Metrics::rmse(real, Predicho), 4),
                       "| MAE =", round(Metrics::mae(real, Predicho), 4)),
      x = "Years", y = unidad_y, colour = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))
  
  print(p2)
  
  # ── 3f. GRÁFICA 3: Dispersión Real vs Predicho ──────────────
  df_sc <- data.frame(Real = real, Predicho = Predicho)
  
  p3 <- ggplot(df_sc, aes(x = Real, y = Predicho)) +
    geom_point(colour = "#e63946", size = 3, alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "#1d3557", linewidth = 0.8) +
    labs(
      title    = paste("Dispersión Real vs. Predicho – MLP –", etiqueta),
      subtitle = "Línea punteada = predicción perfecta (y = x)",
      x = paste("Real –", unidad_y),
      y = paste("Predicho –", unidad_y)
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
  
  print(p3)
  
  # ── 3g. GRÁFICA 4: Importancia de los lags usados ───────────
  lags_usados <- modelo$lags
  if (!is.null(lags_usados) && length(lags_usados) > 0) {
    df_lags <- data.frame(Lag = factor(paste0("t-", lags_usados),
                                       levels = paste0("t-", sort(lags_usados))),
                          Importancia = 1)
    p4 <- ggplot(df_lags, aes(x = Lag, y = Importancia)) +
      geom_col(fill = "#457b9d", colour = "white") +
      labs(
        title    = paste("Lags seleccionados – MLP –", etiqueta),
        subtitle = "Retardos incluidos automáticamente por validación cruzada",
        x = "Retardo (lag)", y = "Incluido en el modelo"
      ) +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1))
    print(p4)
  }
  
  invisible(list(modelo = modelo, forecast = pron, metricas = metricas))
}

# ============================================================
# ── 4. APLICAR MLP A AMBAS VARIABLES ────────────────────────
# ============================================================

res_mlp_vel <- ajustar_mlp(
  serie    = ts_vel,
  etiqueta = "Velocidad del Viento",
  unidad_y = "Vel. Viento (m/s)"
)

res_mlp_dir <- ajustar_mlp(
  serie    = ts_dir,
  etiqueta = "Dirección del Viento",
  unidad_y = "Dir. Viento (sector)"
)

# ============================================================
# ── 5. TABLA COMPARATIVA FINAL ──────────────────────────────
# ============================================================

cat("\n══════════════════════════════════════════════\n")
cat("  RESUMEN COMPARATIVO – MLP\n")
cat("══════════════════════════════════════════════\n")

tabla_final <- data.frame(
  Metrica   = res_mlp_vel$metricas$Metrica,
  Velocidad = res_mlp_vel$metricas$MLP,
  Direccion = res_mlp_dir$metricas$MLP
)

print(tabla_final)
