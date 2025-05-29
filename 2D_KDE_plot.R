# Load library
library(ggplot2)

# Simulated 2D data
set.seed(123)
df <- data.frame(
  X = rnorm(1000),
  Y = rnorm(1000)
)

# 2D KDE Contour Plot
ggplot(df, aes(x = X, y = Y)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", color = "black") +
  scale_fill_viridis_c() +
  labs(title = "2D KDE Contour Plot", x = "X", y = "Y") +
  theme_minimal()

ggplot(df, aes(x = X, y = Y)) +
  geom_density_2d_filled() +
  labs(title = "2D KDE (Filled)", x = "X", y = "Y") +
  theme_minimal()