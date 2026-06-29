# ============================================================
#  EVALUACIÓN DE MÉTRICAS - INTERPOLACIÓN DIRECCIÓN VIENTO
#  Estrategia: Enmascaramiento de datos conocidos (Hold-Out)
# ============================================================

# ============================================================
# ── 1. ESTRATEGIA DE VALIDACIÓN: ENMASCARAMIENTO ────────────
# ============================================================

# Solo podemos evaluar métricas sobre datos QUE SÍ EXISTEN.
# Tomamos el 20% de los datos conocidos, los ocultamos,
# interpolamos como si fueran NA, y comparamos con el real.

set.seed(42)  # reproducibilidad

idx_conocidos_dv <- which(!is.na(dv_completo$Valor))
idx_test_dv      <- sample(idx_conocidos_dv,
                        size = floor(0.20 * length(idx_conocidos_dv)))
idx_test_dv      <- sort(idx_test_dv)

# Guardar valores reales del conjunto de prueba
valores_reales_dv <- dv_completo$Valor[idx_test_dv]

# Crear versión enmascarada (simula datos faltantes)
dv_enmascarado <- dv_completo
dv_enmascarado$Valor[idx_test_dv] <- NA

cat("Datos originales disponibles:", length(idx_conocidos_dv), "meses\n")
cat("Datos enmascarados para test:", length(idx_test_dv), "meses\n")

# ============================================================
# ── 2. CONVERTIR A ts ───────────────────────────────────────
# ============================================================

inicio_dv <- c(year(min(dv_completo$Fecha)), month(min(dv_completo$Fecha)))

ts_enmascarado_dv <- ts(dv_enmascarado$Valor, start = inicio_dv, frequency = 12)

# ============================================================
# ── 3. APLICAR LOS 3 MÉTODOS DE INTERPOLACIÓN ───────────────
# ============================================================

ts_Lineal     <- na.approx(ts_enmascarado_dv, na.rm = FALSE)
ts_spline     <- pmax(na.spline(ts_enmascarado_dv, na.rm = FALSE), 0)
ts_Estacional <- na.interp(ts_enmascarado_dv)

# Extraer solo las predicciones en los índices de test
pred_Lineal     <- as.numeric(ts_Lineal)[idx_test_dv]
pred_spline     <- as.numeric(ts_spline)[idx_test_dv]
pred_Estacional <- as.numeric(ts_Estacional)[idx_test_dv]

# Verificar que no haya NAs en predicciones
cat("\nNAs en predicciones:\n")
cat("  Lineal:    ", sum(is.na(pred_Lineal)), "\n")
cat("  Spline:    ", sum(is.na(pred_spline)), "\n")
cat("  Estacional:", sum(is.na(pred_Estacional)), "\n")

# Filtrar índices donde TODOS los métodos tienen predicción
idx_valido_dv <- !is.na(pred_Lineal) & 
  !is.na(pred_spline) & 
  !is.na(pred_Estacional)

actual          <- valores_reales_dv[idx_valido_dv]
pred_Lineal     <- pred_Lineal[idx_valido_dv]
pred_spline     <- pred_spline[idx_valido_dv]
pred_Estacional <- pred_Estacional[idx_valido_dv]

cat("\nPares válidos para evaluación:", sum(idx_valido_dv), "\n")

# ============================================================
# ── 4. SMAPE (función manual) ───────────────────────────────
# ============================================================

# Ninguna librería la implementa de forma estándar en R
smape_fn <- function(actual, predicted) {
  n <- length(actual)
  100 / n * sum(abs(predicted - actual) /
                  ((abs(actual) + abs(predicted)) / 2))
}

# ============================================================
# ── 5. CALCULAR MÉTRICAS ────────────────────────────────────
# ============================================================

# ─── 5a. Con librería {Metrics} ─────────────────────────────

metricas_Metrics_dv <- data.frame(
  Metrica = c("RMSE", "MSE", "MAE", "MAPE", "SMAPE"),
  
  Lineal = c(
    Metrics::rmse(actual, pred_Lineal),
    Metrics::mse(actual,  pred_Lineal),
    Metrics::mae(actual,  pred_Lineal),
    Metrics::mape(actual, pred_Lineal) * 100,
    smape_fn(actual,      pred_Lineal)
  ),
  
  Spline = c(
    Metrics::rmse(actual, pred_spline),
    Metrics::mse(actual,  pred_spline),
    Metrics::mae(actual,  pred_spline),
    Metrics::mape(actual, pred_spline) * 100,
    smape_fn(actual,      pred_spline)
  ),
  
  Estacional = c(
    Metrics::rmse(actual, pred_Estacional),
    Metrics::mse(actual,  pred_Estacional),
    Metrics::mae(actual,  pred_Estacional),
    Metrics::mape(actual, pred_Estacional) * 100,
    smape_fn(actual,      pred_Estacional)
  )
)

# ─── 5b. Con librería {MLmetrics} ─────────────────────────
# MLmetrics usa convención (y_pred, y_true) — orden inverso a {Metrics}

metricas_MLmetrics_dv <- data.frame(
  Metrica = c("RMSE", "MSE", "MAE", "MAPE", "SMAPE"),
  
  Lineal = c(
    MLmetrics::RMSE(pred_Lineal,     actual),
    MLmetrics::MSE(pred_Lineal,      actual),
    MLmetrics::MAE(pred_Lineal,      actual),
    MLmetrics::MAPE(pred_Lineal,     actual) * 100,
    smape_fn(actual, pred_Lineal)              # misma función manual
  ),
  
  Spline = c(
    MLmetrics::RMSE(pred_spline,     actual),
    MLmetrics::MSE(pred_spline,      actual),
    MLmetrics::MAE(pred_spline,      actual),
    MLmetrics::MAPE(pred_spline,     actual) * 100,
    smape_fn(actual, pred_spline)
  ),
  
  Estacional = c(
    MLmetrics::RMSE(pred_Estacional, actual),
    MLmetrics::MSE(pred_Estacional,  actual),
    MLmetrics::MAE(pred_Estacional,  actual),
    MLmetrics::MAPE(pred_Estacional, actual) * 100,
    smape_fn(actual, pred_Estacional)
  )
)

# ── 6. IMPRIMIR RESULTADOS ──────────────────────────────────

cat("\n══════════════════════════════════════════════\n")
cat("  MÉTRICAS — Librería {Metrics}\n")
cat("══════════════════════════════════════════════\n")
print(metricas_Metrics_dv %>% mutate(across(where(is.numeric), ~round(.x, 5))))

cat("\n══════════════════════════════════════════════\n")
cat("  MÉTRICAS — Librería {MLmetrics}\n")
cat("══════════════════════════════════════════════\n")
print(metricas_MLmetrics_dv %>% mutate(across(where(is.numeric), ~round(.x, 5))))

# Diferencia entre librerías (deben ser ~0 salvo SMAPE que es manual)
cat("\n══════════════════════════════════════════════\n")
cat("  DIFERENCIA entre librerías (debe ser ≈ 0)\n")
cat("══════════════════════════════════════════════\n")
diferencia <- metricas_Metrics_dv
diferencia[, -1] <- abs(metricas_Metrics_dv[, -1] - metricas_MLmetrics_dv[, -1])
print(diferencia %>% mutate(across(where(is.numeric), ~round(.x, 8))))

# ── 7. VISUALIZACIÓN 1: Heatmap de métricas ──────────────────

df_heatmap_dv <- metricas_Metrics_dv %>%
  pivot_longer(cols = -Metrica,
               names_to  = "Metodo",
               values_to = "Valor") %>%
  mutate(
    Metodo  = factor(Metodo,
                     levels = c("Lineal", "Spline", "Estacional")),
    Metrica = factor(Metrica,
                     levels = c("RMSE", "MSE", "MAE", "MAPE", "SMAPE"))
  )

ggplot(df_heatmap_dv, aes(x = Metodo, y = Metrica, fill = Valor)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = round(Valor, 3)), size = 4.5, fontface = "bold") +
  scale_fill_gradient(low = "#d4f1c0", high = "#d62828") +
  labs(
    title    = "Comparación de métricas por método de interpolación",
    subtitle = "Verde = menor error  |  Rojo = mayor error",
    x = "Método de interpolación",
    y = "Métrica",
    fill = "Valor"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    axis.text     = element_text(size = 12),
    legend.position = "right"
  )

# ── 8. VISUALIZACIÓN 2: Real vs Predicted por método ──────────

df_scatter_dv <- data.frame(
  Real       = actual,
  Lineal     = pred_Lineal,
  Spline     = pred_spline,
  Estacional = pred_Estacional
) %>%
  pivot_longer(cols = -Real,
               names_to  = "Metodo",
               values_to = "Predicted")

ggplot(df_scatter_dv, aes(x = Real, y = Predicted, colour = Metodo)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "black", linewidth = 0.8) +
  facet_wrap(~Metodo, ncol = 3) +
  scale_colour_manual(values = c("Lineal"     = "#457B9D",
                                 "Spline"     = "#2A9D8F",
                                 "Estacional" = "#E63946")) +
  labs(
    title    = "Real vs. Predicho — Validación por enmascaramiento",
    subtitle = "La línea punteada representa predicción perfecta (y = x)",
    x = "Valor real (sector)",
    y = "Valor predicho (sector)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 12))

# ── 9. ranking FINAL ──────────────────────────────────────────

cat("\n══════════════════════════════════════════════\n")
cat("  ranking FINAL (menor error = mejor)\n")
cat("══════════════════════════════════════════════\n")

ranking_dv <- metricas_Metrics_dv %>%
  mutate(across(where(is.numeric), ~round(.x, 5))) %>%
  rowwise() %>%
  mutate(Mejor = names(which.min(c(Lineal, Spline, Estacional)))) %>%
  ungroup()

print(ranking_dv)

ranking_dv_2 <- metricas_MLmetrics_dv %>%
  mutate(across(where(is.numeric), ~round(.x, 5))) %>%
  rowwise() %>%
  mutate(Mejor = names(which.min(c(Lineal, Spline, Estacional)))) %>%
  ungroup()

print(ranking_dv_2)
