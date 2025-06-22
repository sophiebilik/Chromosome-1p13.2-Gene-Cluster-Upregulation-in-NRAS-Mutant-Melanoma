# Load libraries
library(readxl)
# Step 1: Read in the data
df <- read_excel("~/Downloads/NRAS.xlsx")
colnames(df)[1:2] <- c("Sample_ID", "NRAS_status")
df$NRAS_status <- as.factor(df$NRAS_status)
# Step 3: Clean column names (important if there are duplicates or spaces)
colnames(df) <- make.names(colnames(df), unique = TRUE)
# Step 4: Convert gene expression columns to numeric
df[, 3:ncol(df)] <- lapply(df[, 3:ncol(df)], function(x) as.numeric(as.character(x)))
genes <- colnames(df)[3:ncol(df)]
results <- data.frame(Gene = genes, p_value = NA, log2FC = NA)
for (i in seq_along(genes)) {
  gene <- genes[i]
  expr_mut <- df[df$NRAS_status == 1, gene]
  expr_wt  <- df[df$NRAS_status == 0, gene]
  
  if (length(na.omit(expr_mut)) > 1 && length(na.omit(expr_wt)) > 1) {
    t_result <- tryCatch(t.test(expr_mut, expr_wt), error = function(e) NULL)
    
    if (!is.null(t_result)) {
      results$p_value[i] <- t_result$p.value
      
      mean_mut <- mean(expr_mut, na.rm = TRUE)
      mean_wt  <- mean(expr_wt, na.rm = TRUE)
      
      if (!is.na(mean_mut) && !is.na(mean_wt) && mean_mut > 0 && mean_wt > 0) {
        results$log2FC[i] <- log2(mean_mut / mean_wt)
      }
    }
  }
}
results$adj_p <- p.adjust(results$p_value, method = "BH")
View(results[order(results$adj_p), ])


library(readxl)

# Step 1: Read file
df_raw <- read_excel("~/Downloads/NRAS.xlsx")

# Step 2: Rename first two columns
colnames(df_raw)[1:2] <- c("Sample_ID", "NRAS_status")

# Step 3: Clean column names
colnames(df_raw) <- make.names(colnames(df_raw), unique = TRUE)

# Step 4: Convert NRAS_status to factor
df_raw$NRAS_status <- as.factor(df_raw$NRAS_status)

# Step 5: Convert gene columns (columns 3 onward) to numeric
df <- df_raw
for (i in 3:ncol(df)) {
  df[[i]] <- suppressWarnings(as.numeric(df[[i]]))  # handles weird formatting
}

# Step 6: Define gene list
genes <- colnames(df)[3:ncol(df)]

# Step 7: Initialize results
results <- data.frame(Gene = genes, p_value = NA, log2FC = NA)

# Step 8: Loop over genes
for (i in seq_along(genes)) {
  gene <- genes[i]
  expr_mut <- df[df$NRAS_status == "1", gene]
  expr_wt  <- df[df$NRAS_status == "0", gene]
  
  if (length(na.omit(expr_mut)) > 1 && length(na.omit(expr_wt)) > 1) {
    test <- tryCatch(t.test(expr_mut, expr_wt), error = function(e) NULL)
    if (!is.null(test)) {
      results$p_value[i] <- test$p.value
      mean_mut <- mean(expr_mut, na.rm = TRUE)
      mean_wt <- mean(expr_wt, na.rm = TRUE)
      if (mean_mut > 0 && mean_wt > 0) {
        results$log2FC[i] <- log2(mean_mut / mean_wt)
      }
    }
  }
}

# Step 9: Adjust p-values
results$adj_p <- p.adjust(results$p_value, method = "BH")

# Step 10: Export results
write.csv(results, "~/Downloads/NRAS_ttest_results.csv", row.names = FALSE)

# Step 11: Optional — View top hits
View(results[order(results$adj_p), ])

str(df[1:5, 1:10])  # just a peek at first few rows and columns
str(df$AHCYL1)
head(df$AHCYL1, 10)
unique(df$NRAS_status)
t.test(df[df$NRAS_status == "1", "AHCYL1", drop = TRUE],
       df[df$NRAS_status == "0", "AHCYL1", drop = TRUE])
df[df$NRAS_status == "1", gene]  # with quotes
results <- data.frame(Gene = genes, p_value = NA, log2FC = NA)

for (i in seq_along(genes)) {
  gene <- genes[i]
  expr_mut <- df[df$NRAS_status == "1", gene]
  expr_wt  <- df[df$NRAS_status == "0", gene]
  
  if (length(na.omit(expr_mut)) > 1 && length(na.omit(expr_wt)) > 1) {
    test <- tryCatch(t.test(expr_mut, expr_wt), error = function(e) NULL)
    if (!is.null(test)) {
      results$p_value[i] <- test$p.value
      mean_mut <- mean(expr_mut, na.rm = TRUE)
      mean_wt  <- mean(expr_wt, na.rm = TRUE)
      if (mean_mut > 0 && mean_wt > 0) {
        results$log2FC[i] <- log2(mean_mut / mean_wt)
      }
    }
  }
}

results$adj_p <- p.adjust(results$p_value, method = "BH")
View(results[order(results$adj_p), ])
results <- data.frame(Gene = genes, p_value = NA, log2FC = NA)

for (i in seq_along(genes)) {
  gene <- genes[i]
  
  # Use [[ ]] to extract numeric vectors properly
  expr_mut <- df[df$NRAS_status == "1", ][[gene]]
  expr_wt  <- df[df$NRAS_status == "0", ][[gene]]
  
  if (length(na.omit(expr_mut)) > 1 && length(na.omit(expr_wt)) > 1) {
    test <- tryCatch(t.test(expr_mut, expr_wt), error = function(e) NULL)
    if (!is.null(test)) {
      results$p_value[i] <- test$p.value
      mean_mut <- mean(expr_mut, na.rm = TRUE)
      mean_wt <- mean(expr_wt, na.rm = TRUE)
      if (!is.na(mean_mut) && !is.na(mean_wt) && mean_mut > 0 && mean_wt > 0) {
        results$log2FC[i] <- log2(mean_mut / mean_wt)
      }
    }
  }
}
results$adj_p <- p.adjust(results$p_value, method = "BH")
View(results[order(results$adj_p), ])
write.csv(results, "~/Downloads/NRAS_ttest_results_FIXED.csv", row.names = FALSE)

# 1. Subset to NRAS-mutant samples
df_mut <- df[df$NRAS_status == "1", ]

# 2. Select only expression columns (genes) — assumes columns 3 onward
expr_matrix <- df_mut[, 3:ncol(df_mut)]

# 3. Compute correlation matrix
cor_matrix <- cor(expr_matrix, use = "pairwise.complete.obs", method = "pearson")

# 4. Visualize: heatmap
heatmap(cor_matrix, 
        main = "Pearson Correlation of 1p13.2 Genes (NRAS-mutant)", 
        col = colorRampPalette(c("blue", "white", "red"))(100), 
        margins = c(10, 10))
# Step 1: Ensure NRAS_status is character for correct subsetting
df$NRAS_status <- as.character(df$NRAS_status)
# Step 2: Subset to NRAS-mutant only
df_mut <- df[df$NRAS_status == "1", ]
# Step 3: Select only gene columns (assumes 3rd to last column = gene expression)
expr_matrix <- df_mut[, 3:ncol(df_mut)]
# Step 4: Convert all expression values to numeric (safe conversion)
expr_matrix <- as.data.frame(lapply(expr_matrix, function(x) as.numeric(as.character(x))))

# Step 5: Compute correlation matrix
cor_matrix <- cor(expr_matrix, use = "pairwise.complete.obs", method = "pearson")
# Step 6: Visualize as heatmap
heatmap(cor_matrix,
        main = "Pearson Correlation of 1p13.2 Genes (NRAS-mutant)",
        col = colorRampPalette(c("blue", "white", "red"))(100),
        margins = c(10, 10))


