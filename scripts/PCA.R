# Read data 
imputed_data <- data.frame(readRDS("/mnt/8TB/Projects/POIAZ/Donagh/fragpipe_wp/outputs/protein_normalized.Rds"))
annotation_df <- readRDS("/mnt/8TB/Projects/POIAZ/Donagh/fragpipe_wp/outputs/experimental_design.Rds")

# Running PCA 
pca <- prcomp(t(imputed_data), center = T, scale. = F)
percentVar <- pca$sdev^2/sum(pca$sdev^2)

PCAvalues <- data.frame(pca$x)
PCAloadings <- data.frame(Variables = rownames(pca$rotation),  pca$rotation)

# adding in ctrl vs stim variable 
PCAvalues$type <- annotation_df$type
PCAvalues$time <- annotation_df$time
PCAvalues$stim <- annotation_df$stim

# PCA (PD1 vs Cntr) 
pdf("/mnt/8TB/Projects/POIAZ/Donagh/fragpipe_wp/outputs/pca_type.pdf", width = 5, height = 4)
ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = type)) +
  geom_point(size = 1.3) + scale_color_viridis(discrete = T) +
  theme_bw() +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),
    y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))
dev.off()

# PCA (time) 
pdf("/mnt/8TB/Projects/POIAZ/Donagh/fragpipe_wp/outputs/pca_time.pdf", width = 5, height = 4)
ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = time)) +
  geom_point(size = 1.3) + scale_color_viridis(discrete = T) +
  theme_bw() +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),
    y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))
dev.off()

# PCA (stim) 
pdf("/mnt/8TB/Projects/POIAZ/Donagh/fragpipe_wp/outputs/pca_stim.pdf", width = 5, height = 4)
ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = stim)) +
  geom_point(size = 1.3) + scale_color_discrete() +
  theme_bw() +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),
    y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))
dev.off()
