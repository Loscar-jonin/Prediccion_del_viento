# ============================================================
#  EVALUACIÓN DE MÉTRICAS - INTERPOLACIÓN VELOCIDAD VIENTO
#  Estrategia: Enmascaramiento de datos conocidos (Hold-Out)
# ============================================================

# ============================================================
# ── 1. ESTRATEGIA DE VALIDACIÓN: ENMASCARAMIENTO ────────────
# ============================================================

# Solo podemos evaluar métricas sobre datos QUE SÍ EXISTEN.
# Tomamos el 20% de los datos conocidos, los ocultamos,
# interpolamos como si fueran NA, y comparamos con el real.

set.seed(42)  # reproducibilidad

idx_conocidos_vv <- which(!is.na(vv_completo$Valor))
idx_test_vv      <- sample(idx_conocidos_vv,
                        size = floor(0.20 * length(idx_conocidos_vv)))
idx_test_vv      <- sort(idx_test_vv)

# Guardar valores reales del conjunto de prueba
valores_reales_vv <- vv_completo$Valor[idx_test_vv]

# Crear versión enmascarada (simula datos faltantes)
vv_enmascarado <- vv_completo
vv_enmascarado$Valor[idx_test_vv] <- NA

cat("Datos originales disponibles:", length(idx_conocidos_vv), "meses\n")
cat("Datos enmascarados para test:", length(idx_test_vv), "meses\n")

# ============================================================
# ── 2. CONVERTIR A ts ───────────────────────────────────────
# ============================================================

inicio_vv <- c(year(min(vv_completo$Fecha)),
            month(min(vv_completo$Fecha)))

ts_enmascarado_vv <- ts(vv_enmascarado$Valor,
                     start     = inicio_vv,
                     frequency = 12)

# ============================================================
# ── 3. APLICAR LOS 3 MÉTODOS DE INTERPOLACIÓN ───────────────
# ============================================================

ts_lineal     <- na.approx(ts_enmascarado_vv, na.rm = FALSE)
ts_spline     <- pmax(na.spline(ts_enmascarado_vv, na.rm = FALSE), 0)
ts_estacional <- na.interp(ts_enmascarado_vv)

# Extraer solo las predicciones en los índices de test
pred_lineal     <- as.numeric(ts_lineal)[idx_test_vv]
pred_spline     <- as.numeric(ts_spline)[idx_test_vv]
pred_estacional <- as.numeric(ts_estacional)[idx_test_vv]

# Verificar que no haya NAs en predicciones
cat("\nNAs en predicciones:\n")
cat("  Lineal:    ", sum(is.na(pred_lineal)), "\n")
cat("  Spline:    ", sum(is.na(pred_spline)), "\n")
cat("  Estacional:", sum(is.na(pred_estacional)), "\n")

# Filtrar índices donde TODOS los métodos tienen predicción
idx_valido_vv <- !is.na(pred_lineal) & 
  !is.na(pred_spline) & 
  !is.na(pred_estacional)

actual          <- valores_reales_vv[idx_valido_vv]
pred_lineal     <- pred_lineal[idx_valido_vv]
pred_spline     <- pred_spline[idx_valido_vv]
pred_estacional <- pred_estacional[idx_valido_vv]

cat("\nPares válidos para evaluación:", sum(idx_valido_vv), "\n")

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

metricas_Metrics_vv <- data.frame(
  Metrica = c("RMSE", "MSE", "MAE", "MAPE", "SMAPE"),
  
  Lineal = c(
    Metrics::rmse(actual, pred_lineal),
    Metrics::mse(actual,  pred_lineal),
    Metrics::mae(actual,  pred_lineal),
    Metrics::mape(actual, pred_lineal) * 100,
    smape_fn(actual,      pred_lineal)
  ),
  
  Spline = c(
    Metrics::rmse(actual, pred_spline),
    Metrics::mse(actual,  pred_spline),
    Metrics::mae(actual,  pred_spline),
    Metrics::mape(actual, pred_spline) * 100,
    smape_fn(actual,      pred_spline)
  ),
  
  Estacional = c(
    Metrics::rmse(actual, pred_estacional),
    Metrics::mse(actual,  pred_estacional),
    Metrics::mae(actual,  pred_estacional),
    Metrics::mape(actual, pred_estacional) * 100,
    smape_fn(actual,      pred_estacional)
  )
)

# ─── 5b. Con librería {MLmetrics} ─────────────────────────
# MLmetrics usa convención (y_pred, y_true) — orden inverso a {Metrics}

metricas_MLmetrics_vv <- data.frame(
  Metrica = c("RMSE", "MSE", "MAE", "MAPE", "SMAPE"),
  
  Lineal = c(
    MLmetrics::RMSE(pred_lineal,     actual),
    MLmetrics::MSE(pred_lineal,      actual),
    MLmetrics::MAE(pred_lineal,      actual),
    MLmetrics::MAPE(pred_lineal,     actual) * 100,
    smape_fn(actual, pred_lineal)              # misma función manual
  ),
  
  Spline = c(
    MLmetrics::RMSE(pred_spline,     actual),
    MLmetrics::MSE(pred_spline,      actual),
    MLmetrics::MAE(pred_spline,      actual),
    MLmetrics::MAPE(pred_spline,     actual) * 100,
    smape_fn(actual, pred_spline)
  ),
  
  Estacional = c(
    MLmetrics::RMSE(pred_estacional, actual),
    MLmetrics::MSE(pred_estacional,  actual),
    MLmetrics::MAE(pred_estacional,  actual),
    MLmetrics::MAPE(pred_estacional, actual) * 100,
    smape_fn(actual, pred_estacional)
  )
)

# ── 6. IMPRIMIR RESULTADOS ──────────────────────────────────

cat("\n══════════════════════════════════════════════\n")
cat("  MÉTRICAS — Librería {Metrics}\n")
cat("══════════════════════════════════════════════\n")
print(metricas_Metrics_vv %>% mutate(across(where(is.numeric), ~round(.x, 5))))

cat("\n══════════════════════════════════════════════\n")
cat("  MÉTRICAS — Librería {MLmetrics}\n")
cat("══════════════════════════════════════════════\n")
print(metricas_MLmetrics_vv %>% mutate(across(where(is.numeric), ~round(.x, 5))))

# Diferencia entre librerías (deben ser ~0 salvo SMAPE que es manual)
cat("\n══════════════════════════════════════════════\n")
cat("  DIFERENCIA entre librerías (debe ser ≈ 0)\n")
cat("══════════════════════════════════════════════\n")
diferencia <- metricas_Metrics_vv
diferencia[, -1] <- abs(metricas_Metrics_vv[, -1] - metricas_MLmetrics_vv[, -1])
print(diferencia %>% mutate(across(where(is.numeric), ~round(.x, 8))))

# ── 7. VISUALIZACIÓN 1: Heatmap de métricas ──────────────────

df_heatmap_vv <- metricas_Metrics_vv %>%
  pivot_longer(cols = -Metrica,
               names_to  = "Metodo",
               values_to = "Valor") %>%
  mutate(
    Metodo  = factor(Metodo,
                     levels = c("Lineal", "Spline", "Estacional")),
    Metrica = factor(Metrica,
                     levels = c("RMSE", "MSE", "MAE", "MAPE", "SMAPE"))
  )

ggplot(df_heatmap_vv, aes(x = Metodo, y = Metrica, fill = Valor)) +
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

# ── 8. VISUALIZACIÓN 2: Real vs Predicho por método ──────────

df_scatter_vv <- data.frame(
  Real       = actual,
  Lineal     = pred_lineal,
  Spline     = pred_spline,
  Estacional = pred_estacional
) %>%
  pivot_longer(cols = -Real,
               names_to  = "Metodo",
               values_to = "Predicho")

ggplot(df_scatter_vv, aes(x = Real, y = Predicho, colour = Metodo)) +
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
    x = "Valor real (m/s)",
    y = "Valor predicho (m/s)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 12))

# ── 9. ranking FINAL ──────────────────────────────────────────

cat("\n══════════════════════════════════════════════\n")
cat("  ranking FINAL (menor error = mejor)\n")
cat("══════════════════════════════════════════════\n")

ranking_vv <- metricas_Metrics_vv %>%
  mutate(across(where(is.numeric), ~round(.x, 5))) %>%
  rowwise() %>%
  mutate(Mejor = names(which.min(c(Lineal, Spline, Estacional)))) %>%
  ungroup()

print(ranking_vv)

ranking_vv_2 <- metricas_MLmetrics_vv %>%
  mutate(across(where(is.numeric), ~round(.x, 5))) %>%
  rowwise() %>%
  mutate(Mejor = names(which.min(c(Lineal, Spline, Estacional)))) %>%
  ungroup()

print(ranking_vv_2)
