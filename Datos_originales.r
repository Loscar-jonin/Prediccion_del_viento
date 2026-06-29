# ============================================================
# ── 1. LIBRERÍAS ────────────────────────────────────────────
# ============================================================

library(Metrics)     # rmse, mae, mse, mape
library(MLmetrics)   # RMSE, MAE, MAPE, SMAPE (implementación ML)
library(forecast)
library(dplyr)
library(ggplot2)
library(lubridate)
library(readr)
library(zoo)
library(tidyr)
library(brnn)       # brnn(): redes neuronales con regularización bayesiana
library(nnfor)      # mlp(): redes neuronales para series de tiempo

# ============================================================
# ── 2. CARGAR DATOS ─────────────────────────────────────────
# ============================================================

vel_data <- read_csv("Vel_viento/descargaDhime.csv")
View(vel_data)
dir_data <- read_csv("Dir_viento/descargaDhime.csv")
View(dir_data)

# ============================================================
# ── 3. ELIMINAR COLUMNAS INNECESARIAS ───────────────────────
# ============================================================

vv_data <- vel_data %>% select(Variable, Fecha, Unidad, Valor)
View(vv_data)
dv_data <- dir_data %>% select(Variable, Fecha, Unidad, Valor)
View(dv_data)

# ============================================================
# ── 4. RESUMEN ESTADÍSTICO ──────────────────────────────────
# ============================================================

summary(vv_data$Valor)
summary(dv_data$Valor)

# ============================================================
# ── 5. GRÁFICAS LOS DATOS ───────────────────────────────────
# ============================================================

# ── 5a. HISTOGRAMAS ─────────────────────────────────────────

hist(vel_data$Valor)
hist(log(vel_data$Valor))

hist(dir_data$Valor)
hist(log(dir_data$Valor))

# ── 5b. SERIE DE TIEMPO ─────────────────────────────────────

# ── (DIR. VIENTO) ───────────────────────────────────────────

ggplot(dv_data, aes(x = Fecha, y = Valor)) +
  geom_line(color = "blue") +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  scale_x_date(
    date_breaks = "2 years",
    date_labels = "%Y"
  ) +
  labs(title = "Análisis de datos de la dirección",
       subtitle = "Tendencia temporal",
       x = "Year",
       y = "Wind Dir. (sector)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ── (VEL. VIENTO) ───────────────────────────────────────────

ggplot(vv_data, aes(x = Fecha, y = Valor)) +
  geom_line(color = "blue") +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  scale_x_date(
    date_breaks = "2 years",
    date_labels = "%Y"
  ) +
  labs(title = "Análisis de datos de la velocidad",
       subtitle = "Tendencia temporal",
       x = "Year",
       y = "Wind Vel. (m/s)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
