# ============================================================
#  MODELO BNN – RED NEURONAL BAYESIANA PARA SERIES DE TIEMPO
#  Predicción de Velocidad y Dirección del Viento
#  Frecuencia: Mensual  |  Horizonte de predicción: 12 meses
#
#  Usa el paquete {brnn} (Bayesian Regularized Neural Networks):
#    – Regularización bayesiana de MacKay para los pesos
#    – Evita sobreajuste sin necesidad de conjunto de validación
#    – Los hiperparámetros α (regularización) y β (precisión de ruido)
#      se estiman automáticamente desde los datos
#    – Arquitectura: red feedforward con 1 capa oculta
#
#  Prerrequisito: ejecutar datos_originales.R,
#                 Interpolacion_vel_viento.R e
#                 Interpolacion_dir_viento.R
# ============================================================

# ============================================================
# ── 1. FUNCIÓN: CONSTRUIR MATRIZ DE LAGS ────────────────────
# ============================================================
# Los modelos de ML no consumen directamente objetos ts.
# Transformamos la serie en una matriz de características:
#   X = [y(t-1), y(t-2), ..., y(t-p)] → Y = y(t)

construir_lags <- function(serie, p) {
  x_vec <- as.numeric(serie)
  n     <- length(x_vec)

  X <- matrix(NA, nrow = n - p, ncol = p)
  for (i in seq_len(p)) {
    X[, i] <- x_vec[(p - i + 1):(n - i)]
  }
  colnames(X) <- paste0("lag", seq_len(p))

  Y <- x_vec[(p + 1):n]
  list(X = X, Y = Y)
}

# ============================================================
# ── 2. RECONSTRUIR SERIES INTERPOLADAS ──────────────────────
# ============================================================

inicio_vv <- c(year(min(vv_completo$Fecha)), month(min(vv_completo$Fecha)))
inicio_dv <- c(year(min(dv_completo$Fecha)), month(min(dv_completo$Fecha)))

ts_vel <- na.interp(ts(vv_completo$Valor, start = inicio_vv, frequency = 12))
ts_dir <- na.interp(ts(dv_completo$Valor, start = inicio_dv, frequency = 12))

cat("Serie velocidad – longitud:", length(ts_vel), "meses\n")
cat("Serie dirección – longitud:", length(ts_dir), "meses\n")

# ============================================================
# ── 3. PARÁMETROS GLOBALES ──────────────────────────────────
# ============================================================

H       <- 12    # horizonte de predicción (meses)
P_LAGS  <- 12    # retardos usados como predictores (1 ciclo anual)
N_NEURO <- 5     # neuronas en la capa oculta
set.seed(42)

# ============================================================
# ── 4. FUNCIÓN AUXILIAR: NORMALIZACIÓN ──────────────────────
# ============================================================

normalizar <- function(x) {
  mu  <- mean(x, na.rm = TRUE)
  sig <- sd(x,   na.rm = TRUE)
  list(z = (x - mu) / sig, mu = mu, sig = sig)
}

desnormalizar <- function(z, mu, sig) z * sig + mu

# ============================================================
# ── 5. FUNCIÓN PRINCIPAL: AJUSTAR Y EVALUAR BNN ─────────────
# ============================================================

ajustar_bnn <- function(serie, etiqueta, unidad_y) {

  n        <- length(serie)
  n_train  <- n - H

  # ── 5a. DIVIDIR EN ENTRENAMIENTO Y PRUEBA ───────────────────
  ts_train_raw <- as.numeric(window(serie, end   = time(serie)[n_train]))
  ts_test_raw  <- as.numeric(window(serie, start = time(serie)[n_train + 1]))

  cat("\n══════════════════════════════════════════════\n")
  cat(" BNN –", etiqueta, "\n")
  cat("══════════════════════════════════════════════\n")
  cat(" Observaciones totales :", n, "\n")
  cat(" Entrenamiento         :", n_train, "meses\n")
  cat(" Prueba                :", H, "meses\n")
  cat(" Lags usados           :", P_LAGS, "\n")
  cat(" Neuronas ocultas      :", N_NEURO, "\n")

  # ── 5b. NORMALIZACIÓN (se normaliza solo sobre el train) ────
  norm_obj   <- normalizar(ts_train_raw)
  ts_train_n <- norm_obj$z

  # Construir lags sobre la serie de entrenamiento normalizada
  datos_lags  <- construir_lags(ts_train_n, P_LAGS)
  X_train     <- datos_lags$X
  Y_train     <- datos_lags$Y

  # ── 5c. AJUSTE BNN ──────────────────────────────────────────
  # neurons: neuronas en la capa oculta (1 capa oculta + tangente hiperbólica)
  # epochs: iteraciones máximas del algoritmo de Levenberg-Marquardt
  # verbose: mostrar progreso de convergencia

  modelo <- brnn(
    x       = X_train,
    y       = Y_train,
    neurons = N_NEURO,
    epochs  = 1000,
    verbose = FALSE
  )

  cat("\nHiperparámetros estimados bayesianamente:\n")
  cat("  α (regularización de pesos) :", round(modelo$alpha, 6), "\n")
  cat("  β (precisión del ruido)     :", round(modelo$beta,  6), "\n")
  cat("  γ (parámetros efectivos)    :", round(modelo$gamma, 4), "\n")

  # ── 5d. PRONÓSTICO ITERATIVO (h pasos adelante) ──────────────
  # Dado que brnn predice un paso, aplicamos predicción recursiva:
  # y(t+1) = f(y(t), y(t-1), ..., y(t-p+1))
  # y(t+2) = f(ŷ(t+1), y(t), ..., y(t-p+2)), etc.

  # Semilla de predicción: últimos P_LAGS valores del train normalizado
  lag_window <- tail(ts_train_n, P_LAGS)

  predicho_n <- numeric(H)
  for (h_i in seq_len(H)) {
    x_new         <- matrix(rev(lag_window), nrow = 1)
    colnames(x_new) <- paste0("lag", seq_len(P_LAGS))
    pred_n        <- predict(modelo, x_new)
    predicho_n[h_i] <- pred_n
    lag_window    <- c(lag_window[-1], pred_n)   # deslizar ventana
  }

  # Desnormalizar predicciones
  predicho <- desnormalizar(predicho_n, norm_obj$mu, norm_obj$sig)

  # Ajuste sobre el conjunto de entrenamiento (desnormalizado)
  fitted_n <- predict(modelo, X_train)
  fitted   <- desnormalizar(fitted_n, norm_obj$mu, norm_obj$sig)

  real <- ts_test_raw

  # ── 5e. MÉTRICAS ────────────────────────────────────────────
  metricas <- data.frame(
    Metrica = c("RMSE", "MSE", "MAE", "MAPE (%)", "SMAPE (%)"),
    BNN     = c(
      Metrics::rmse(real, predicho),
      Metrics::mse(real,  predicho),
      Metrics::mae(real,  predicho),
      Metrics::mape(real, predicho) * 100,
      smape_fn(real, predicho)
    )
  ) %>% mutate(across(where(is.numeric), ~round(.x, 5)))

  cat("\n── Métricas de predicción ──────────────────\n")
  print(metricas)

  # ── 5f. FECHAS PARA LAS GRÁFICAS ────────────────────────────
  fechas_hist  <- as.Date(as.yearmon(time(serie)))
  fechas_train <- fechas_hist[seq_len(n_train)]
  fechas_test  <- fechas_hist[(n_train + 1):n]
  fechas_fit   <- fechas_train[(P_LAGS + 1):n_train]   # lags consumen P_LAGS obs

  # ── 5g. GRÁFICA 1: Serie + Ajuste + Pronóstico ──────────────
  df_hist <- data.frame(Fecha = fechas_hist,
                        Valor = as.numeric(serie),
                        Tipo  = "Histórico")

  df_fit  <- data.frame(Fecha = fechas_fit,
                        Valor = fitted,
                        Tipo  = "Ajuste (train)")

  df_pred <- data.frame(Fecha = fechas_test,
                        Valor = predicho,
                        Tipo  = "Pronóstico BNN")

  df_all <- bind_rows(df_hist, df_fit, df_pred)

  p1 <- ggplot(df_all, aes(x = Fecha, y = Valor, colour = Tipo)) +
    geom_line(linewidth = 0.8) +
    geom_point(data = df_pred, size = 2.5) +
    scale_colour_manual(
      values = c("Histórico"     = "#1d3557",
                 "Ajuste (train)"= "#2a9d8f",
                 "Pronóstico BNN"= "#e63946")
    ) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(
      x = "Fecha", y = unidad_y, colour = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))

  print(p1)

  # ── 5h. GRÁFICA 2: Real vs Predicho (prueba) ────────────────
  df_comp <- data.frame(
    Fecha    = fechas_test,
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
      title    = paste("Real vs. Predicho – BNN –", etiqueta),
      subtitle = paste("RMSE =", round(Metrics::rmse(real, predicho), 4),
                       "| MAE =", round(Metrics::mae(real, predicho), 4)),
      x = "Fecha", y = unidad_y, colour = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))

  print(p2)

  # ── 5i. GRÁFICA 3: Dispersión Real vs Predicho ──────────────
  df_sc <- data.frame(Real = real, Predicho = predicho)

  p3 <- ggplot(df_sc, aes(x = Real, y = Predicho)) +
    geom_point(colour = "#e63946", size = 3, alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "#1d3557", linewidth = 0.8) +
    labs(
      title    = paste("Dispersión Real vs. Predicho – BNN –", etiqueta),
      subtitle = "Línea punteada = predicción perfecta (y = x)",
      x = paste("Real –", unidad_y),
      y = paste("Predicho –", unidad_y)
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

  print(p3)

  # ── 5j. GRÁFICA 4: Residuos del entrenamiento ───────────────
  residuos <- Y_train - fitted_n
  df_res   <- data.frame(
    Fecha    = fechas_fit,
    Residuo  = desnormalizar(residuos, 0, norm_obj$sig)  # escala original
  )

  p4 <- ggplot(df_res, aes(x = Fecha, y = Residuo)) +
    geom_line(colour = "#457b9d", linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "#e63946", linewidth = 0.8) +
    labs(
      title    = paste("Residuos BNN – Entrenamiento –", etiqueta),
      subtitle = "Los residuos deben distribuirse aleatoriamente alrededor de 0",
      x = "Fecha", y = paste("Residuo –", unidad_y)
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

  print(p4)

  invisible(list(
    modelo   = modelo,
    predicho = predicho,
    metricas = metricas
  ))
}

# ============================================================
# ── 6. APLICAR BNN A AMBAS VARIABLES ────────────────────────
# ============================================================

res_bnn_vel <- ajustar_bnn(
  serie    = ts_vel,
  etiqueta = "Velocidad del Viento",
  unidad_y = "Vel. Viento (m/s)"
)

res_bnn_dir <- ajustar_bnn(
  serie    = ts_dir,
  etiqueta = "Dirección del Viento",
  unidad_y = "Dir. Viento (sector)"
)

# ============================================================
# ── 7. TABLA COMPARATIVA FINAL ──────────────────────────────
# ============================================================

cat("\n══════════════════════════════════════════════\n")
cat("  RESUMEN COMPARATIVO – BNN\n")
cat("══════════════════════════════════════════════\n")

tabla_final <- data.frame(
  Metrica   = res_bnn_vel$metricas$Metrica,
  Velocidad = res_bnn_vel$metricas$BNN,
  Direccion = res_bnn_dir$metricas$BNN
)

print(tabla_final)
