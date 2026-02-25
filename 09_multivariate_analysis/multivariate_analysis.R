# =============================================================================
# MULTIVARIATE ANALYSIS IN R
# =============================================================================
# Topics: PCA, EFA, k-means, hierarchical clustering, LDA, MANOVA, CCA
# Author: Ebenezer Adjartey
# =============================================================================

pkgs <- c("ggplot2","dplyr","tidyr","factoextra","FactoMineR","psych",
          "cluster","MASS","CCA","car")
for (p in pkgs) if (!requireNamespace(p,quietly=TRUE)) install.packages(p)

library(ggplot2); library(factoextra); library(FactoMineR)
library(psych);   library(cluster);    library(MASS)
set.seed(42)

# ── 1. Generate Correlated Dataset ────────────────────────────────────────────
n  <- 200
f1 <- rnorm(n); f2 <- rnorm(n); f3 <- rnorm(n)
X  <- cbind(
  reading    = 0.8*f1 + 0.1*f2 + rnorm(n,.3),
  vocabulary = 0.7*f1 + 0.2*f2 + rnorm(n,.3),
  math       = 0.1*f1 + 0.8*f2 + rnorm(n,.3),
  statistics = 0.2*f1 + 0.9*f2 + rnorm(n,.3),
  memory1    = 0.1*f1 + 0.1*f2 + 0.8*f3 + rnorm(n,.3),
  memory2    = 0.2*f1 + 0.1*f2 + 0.7*f3 + rnorm(n,.3),
  noise1     = rnorm(n),
  noise2     = rnorm(n)
)
df <- as.data.frame(X)
group <- cut(f1+f2, breaks=3, labels=c(1,2,3))
df$group <- group

cat("Correlation matrix:\n")
print(round(cor(df[,1:8]), 2))

# ── 2. Principal Component Analysis ──────────────────────────────────────────
cat("\n=== PCA ===\n")
X_sc <- scale(df[,1:8])
pca  <- prcomp(X_sc, center=FALSE, scale.=FALSE)

# Eigenvalues and variance explained
eig_df <- data.frame(
  PC        = paste0("PC",1:8),
  Eigenvalue= round(pca$sdev^2, 4),
  Var_pct   = round(pca$sdev^2/sum(pca$sdev^2)*100, 2),
  Cum_pct   = round(cumsum(pca$sdev^2/sum(pca$sdev^2)*100), 2)
)
print(eig_df)
cat(sprintf("\nComponents for 80%% variance: %d\n\n",
            min(which(eig_df$Cum_pct>=80))))

# Loadings
loadings_df <- data.frame(round(pca$rotation[,1:4], 3))
print(loadings_df)

dir.create("09_multivariate_analysis", showWarnings=FALSE)

# Scree plot
p_scree <- fviz_eig(pca, addlabels=TRUE, ylim=c(0,40),
                     main="Scree Plot (PCA)")
ggsave("09_multivariate_analysis/scree_plot.png", p_scree, width=7, height=5)

# Biplot
p_biplot <- fviz_pca_biplot(pca, label="var",
                              habillage=df$group, palette=c("red","green","blue"),
                              title="PCA Biplot", repel=TRUE)
ggsave("09_multivariate_analysis/pca_biplot.png", p_biplot, width=8, height=6)

# ── 3. Exploratory Factor Analysis ────────────────────────────────────────────
cat("=== EXPLORATORY FACTOR ANALYSIS (3 factors) ===\n")
fa_result <- fa(X_sc, nfactors=3, rotate="varimax", fm="ml")
print(fa_result$loadings, cutoff=0.3)
cat("\nCommunalities:\n")
print(round(fa_result$communality, 4))
cat("\nVariance explained per factor:\n")
print(round(fa_result$Vaccounted, 3))

# ── 4. K-Means Clustering ─────────────────────────────────────────────────────
cat("\n=== K-MEANS CLUSTERING ===\n")
# Elbow + silhouette
sil_scores <- sapply(2:8, function(k) {
  km <- kmeans(X_sc, centers=k, nstart=20)
  mean(silhouette(km$cluster, dist(X_sc))[,"sil_width"])
})
best_k <- which.max(sil_scores) + 1
cat(sprintf("Silhouette scores for k=2..8: %s\n", paste(round(sil_scores,3),collapse=" ")))
cat(sprintf("Best k = %d\n", best_k))

km_final <- kmeans(X_sc, centers=best_k, nstart=20)
df$cluster <- factor(km_final$cluster)
cat("Cluster sizes:\n"); print(table(df$cluster))

# Plot clusters
p_clust <- fviz_cluster(km_final, data=X_sc, ellipse.type="convex",
                          palette="Set1", ggtheme=theme_minimal(),
                          main=sprintf("K-Means (k=%d)",best_k))
ggsave("09_multivariate_analysis/kmeans_plot.png", p_clust, width=7, height=6)

# ── 5. Hierarchical Clustering ────────────────────────────────────────────────
cat("\n=== HIERARCHICAL CLUSTERING ===\n")
hc <- hclust(dist(X_sc[1:50,]), method="ward.D2")
cat("Height (last 5 merges):\n"); print(round(tail(hc$height,5),3))

p_dend <- fviz_dend(hc, k=3, main="Ward Dendrogram",
                     k_colors=c("red","blue","green"),
                     cex=0.4, rect=TRUE)
ggsave("09_multivariate_analysis/dendrogram.png", p_dend, width=9, height=5)

# ── 6. Linear Discriminant Analysis ──────────────────────────────────────────
cat("\n=== LINEAR DISCRIMINANT ANALYSIS ===\n")
lda_fit <- lda(group ~ ., data=df[,c(1:8,9)])
print(lda_fit)

# Proportion of between-group variance explained
cat("Proportion of between-group var explained per LD:\n")
print(round(lda_fit$svd^2 / sum(lda_fit$svd^2), 4))

# Plot
lda_pred <- predict(lda_fit)
lda_df   <- data.frame(lda_pred$x, group=df$group)
p_lda <- ggplot(lda_df, aes(x=LD1, y=LD2, color=group)) +
  geom_point(alpha=0.6, size=2) +
  stat_ellipse(level=0.9) +
  labs(title="LDA: First Two Discriminant Functions") +
  theme_minimal()
ggsave("09_multivariate_analysis/lda_plot.png", p_lda, width=7, height=5)

# ── 7. MANOVA ─────────────────────────────────────────────────────────────────
cat("\n=== MANOVA ===\n")
manova_fit <- manova(cbind(reading,vocabulary,math,statistics) ~ group, data=df)
print(summary(manova_fit, test="Pillai"))
cat("\nUnivariate ANOVAs:\n")
print(summary.aov(manova_fit))

# ── 8. Canonical Correlation Analysis ────────────────────────────────────────
cat("\n=== CANONICAL CORRELATION ANALYSIS ===\n")
X1 <- as.matrix(df[,c("reading","vocabulary")])
X2 <- as.matrix(df[,c("math","statistics")])
cc  <- cancor(X1, X2)
cat("Canonical correlations:", round(cc$cor, 4), "\n")
cat("Coef set 1 (reading,vocabulary -> canonical var):\n")
print(round(cc$xcoef, 4))
cat("Coef set 2 (math,statistics -> canonical var):\n")
print(round(cc$ycoef, 4))

cat("\n=== MULTIVARIATE ANALYSIS COMPLETE ===\n")
