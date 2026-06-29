# ============================================================
#  INTERPOLACIÓN DE DATOS FALTANTES - DIRECCIÓN DEL VIENTO
#  Dataset: descargaDhime.csv  |  Frecuencia: Mensual
# ============================================================

# ============================================================
# ── 1. LIMPIEZA DE DATOS ────────────────────────────────────
# ============================================================

dv_data <- dir_data %>%
  select(Variable, Fecha, Unidad, Valor) %>%
  mutate(Fecha = as.Date(Fecha))

# ============================================================
# ── 2. CREAR SECUENCIA COMPLETA DE MESES ────────────────────
# ============================================================

# Sin este paso, los NAs no existen: simplemente no hay filas para los meses ausentes.
# Zoo/forecast los necesita explícitos.

fecha_completa <- data.frame(
  Fecha = seq(
    from = floor_date(min(dv_data$Fecha), "month"),
    to   = floor_date(max(dv_data$Fecha), "month"),
    by   = "month"
  )
)

# El comando "Join" es para insertar los NAs donde falta el mes
dv_completo <- fecha_completa %>%
  left_join(dv_data, by = "Fecha")

cat("Filas totales después del join:", nrow(dv_completo), "\n")
cat("Meses con NA:", sum(is.na(dv_completo$Valor)), "\n")

# ============================================================
# ── 3. CONVERTIR A OBJETO SERIE DE TIEMPO (ts) ──────────────
# ============================================================

# Necesario para métodos de forecast y para identificar
# la estacionalidad mensual (period = 12)

inicio_dv <- c(year(min(dv_completo$Fecha)),
            month(min(dv_completo$Fecha)))

ts_direccion <- ts(
  data      = dv_completo$Valor,
  start     = inicio_dv,
  frequency = 12          # mensual → estacionalidad anual
)

# ============================================================
# ── 4. MÉTODOS DE INTERPOLACIÓN ─────────────────────────────
# ============================================================

# ── 4a. Lineal (zoo) ────────────────────────────────────────
#   Bueno para gaps cortos (1-3 meses). Conecta puntos conocidos
#   con una línea recta. Rápido pero ignora la estacionalidad.
ts_lineal <- na.approx(ts_direccion, na.rm = FALSE)

# ── 4b. Spline cúbico (zoo) ─────────────────────────────────
#   Curvas suaves entre puntos. Mejor que lineal para datos
#   con tendencia no lineal, pero puede oscilar en gaps largos.
ts_spline <- na.spline(ts_direccion, na.rm = FALSE)
ts_spline <- pmax(ts_spline, 0)   # la velocidad no puede ser negativa

# ── 4c. Interpolación estacional (forecast) ─────────────────
#   Usa STL (descomposición estacional) + interpolación de Kalman.
#   Es el único método que respeta el patrón anual del viento,
#   por eso funciona bien incluso para gaps de varios años.
ts_estacional <- na.interp(ts_direccion)

# ============================================================
# ── 5. COMPARAR MÉTODOS VISUALMENTE ─────────────────────────
# ============================================================

# Convertir a data.frames para ggplot
df2_plot <- dv_completo %>%
  mutate(
    Lineal     = as.numeric(ts_lineal),
    Spline     = as.numeric(ts_spline),
    Estacional = as.numeric(ts_estacional)
  )

ggplot(df2_plot, aes(x = Fecha)) +
  geom_line(aes(y = Estacional, colour = "Estacional (STL)"), linewidth = 0.8) +
  geom_line(aes(y = Lineal, colour = "Lineal"), linewidth = 0.6, linetype = "dashed") +
  geom_line(aes(y = Spline, colour = "Spline"), linewidth = 0.6, linetype = "dotted") +
  geom_point(aes(y = Valor), colour = "black", size = 1.2, alpha = 0.7, na.rm = TRUE) +
  scale_colour_manual(
    values = c("Estacional (STL)" = "#E63946",
               "Lineal"           = "#457B9D",
               "Spline"           = "#2A9D8F")
  ) +
  scale_x_date(
    date_breaks = "2 years",
    date_labels = "%Y"
  ) +
  labs(
    title    = "Dirección del Viento – Interpolación de datos faltantes",
    subtitle = paste("Estación: Aeropuerto Simón Bolívar |",
                     sum(is.na(dv_completo$Valor)), "meses interpolados"),
    x        = "Fecha",
    y        = "Dir. Viento (sector)",
    colour   = "Método"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position  = "bottom",
    axis.text.x      = element_text(angle = 45, hjust = 1)
  )
# ============================================================
# ── 6. RESUMEN FINAL ────────────────────────────────────────
# ============================================================

resumen_dv <- summary(dv_completo$Valor)
print.table(resumen_dv)
