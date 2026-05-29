library(ggplot2)
library(tidyverse)
library(ggpubr)

label_environment <- function(x) {
  x <- gsub("_", " ", x)
  x <- stringr::str_to_title(x)

  dplyr::recode(
    x,
    "Less Than Hs" = "< high school",
    "More Than Hs" = "> high school",
    "Amphetamines" = "Amphet.",
    "Previously Married" = "Prev. married",
    .default = x
  )
}

theme_enrichment <- function(base_size = 8) {
  ggpubr::theme_pubr(base_size = base_size) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title = element_text(size = base_size + 1),
      axis.text = element_text(size = base_size),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(size = base_size, face = "bold"),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size),
      legend.key.height = unit(0.35, "cm"),
      legend.key.width = unit(0.45, "cm"),
      panel.spacing = unit(0.45, "lines"),
      plot.margin = margin(4, 4, 4, 4)
    )
}

save_plot <- function(p, fn, w, h){
    for(ext in c(".pdf", ".png")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

load_env_enrichment <- function(pop){
    return(data.table::fread(paste0("vmr_enrichment_analysis-", pop, ".txt")))
}
memENRICH <- memoise::memoise(load_env_enrichment)

gen_data <- function(pop){
    err = 0.0000001
    dt <- memENRICH(pop) %>% mutate(across(where(is.character), as.factor)) %>%
        mutate(h2_Category=fct_relevel(h2_Category, rev), `-log10(FDR)`= -log10(FDR),
               `OR Percentile`= OR / (1+OR), p.fdr.sig=FDR < 0.05,
               `log2(OR)` = log2(OR+err),
               p.fdr.cat=cut(FDR, breaks=c(1,0.05,0.01,0.005,0),
                             labels=c("<= 0.005","<= 0.01","<= 0.05","> 0.05"),
                             include.lowest=TRUE))
        #mutate(across(Direction, factor, levels=c("All", "Down", "Up")))
    #levels(dt$Direction) <- c("All", "Heritable Bias", "Non-heritable Bias")
    return(dt)
}
memDF <- memoise::memoise(gen_data)

enrichment_grid <- function(pop) {
  memDF(pop) %>% distinct(Tissue, Test, h2_Category)
}

plot_tile <- function(pop, label, w, h){
  df <- memDF(pop) %>%
    filter(Env == label) %>%
    right_join(enrichment_grid(pop), 
               by = c("Tissue", "Test", "h2_Category")) %>%
    mutate(
      p.fdr.sig = replace_na(p.fdr.sig, FALSE),
      # Keep non-finite log2(OR) as NA so scale_fill maps them to na.value (light grey)
      `log2(OR)` = if_else(is.finite(`log2(OR)`), `log2(OR)`, NA_real_)
    )

  fin <- is.finite(df$`log2(OR)`)
  if (!any(fin)) {
    y0 <- -1
    y1 <- 1
  } else {
    y0 <- min(df$`log2(OR)`[fin], na.rm = TRUE) - 0.1
    y1 <- max(df$`log2(OR)`[fin], na.rm = TRUE) + 0.1
  }

  tile_plot <- df %>%
    ggplot(aes(y = Test, x = h2_Category, fill = `log2(OR)`)) +
    geom_tile(color = "grey50", linewidth = 0.2) +
    geom_text(aes(label = ifelse(p.fdr.sig,
                                 format(round(-log10(FDR),1), nsmall=1), "")),
              color = "black", size = 5) +
    scale_fill_gradientn(colors = c("#3B6EA8", "white", "#C65146"),
                         values = scales::rescale(c(y0, 0, y1)),
                         limits = c(y0, y1),
                         na.value = "grey90",
                         name = "log2(OR)") +
    facet_grid(. ~ Tissue) +
    labs(x = "Heritability Category", y = "Test") +
    theme_minimal(base_size = 20) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      axis.title = element_text(face = "bold", size = 28),
      strip.text = element_text(face = "bold", size = 22)
    )
  
  save_plot(tile_plot, paste0(tolower(label), "_tileplot_enrichment-", pop), w, h)
}

plot_all_env_heatmap <- function(pop, w = 8, h = 10){

  env_order <- c(
    "smoking",
    "nicotine",
    "cocaine",
    "morphine",
    "codeine",
    "amphetamines",
    "ethanol",
    "less_than_hs",
    "more_than_hs",
    "previously_married",
    "single"
  )

  df <- memDF(pop) %>%
    filter(as.character(h2_Category) != "Low prediction") %>%
    mutate(
      sig = !is.na(FDR) & FDR < 0.05,
      sig_level = case_when(
        is.na(FDR) ~ NA_character_,
        FDR <= 0.0001 ~ "****",
        FDR <= 0.001  ~ "***",
        FDR <= 0.01   ~ "***",
        FDR <= 0.05   ~ "*",
        TRUE ~ ""
      ),
      log2_or = `log2(OR)`,
      Env = factor(Env, levels = rev(env_order)),
      Test = factor(
        Test,
        levels = c("Linear", "DMR", "Both"),
        labels = c("VMR", "DMR", "VMR + DMR")
      ),
      h2_Category = factor(h2_Category, levels = c("Heritable", "Non-heritable"))
    )

  or_lim <- max(abs(df$log2_or), na.rm = TRUE)
  or_lim <- ceiling(or_lim * 2) / 2

  p <- ggplot(
    df,
    aes(
      x = h2_Category,
      y = Env,
      fill = log2_or
    )
  ) +

    geom_tile(
      color = "white",
      linewidth = 0.25,
      width = 0.95,
      height = 0.95
    ) +

    geom_text(
      aes(label = sig_level),
      color = "black",
      size = 4
    ) +

    facet_grid(
      Test ~ Tissue
    ) +

    scale_fill_gradient2(
      low = "#3B6EA8",
      mid = "white",
      high = "#C65146",
      midpoint = 0,
      limits = c(-or_lim, or_lim),
      na.value = "grey90",
      name = "log2(OR)"
    ) +

    scale_y_discrete(
      labels = label_environment
    ) +

    labs(
      x = "H2 Category",
      y = NULL
    ) +

    theme_enrichment(base_size = 10) +

    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      legend.position = "right",
      panel.grid = element_blank(),
      strip.text = element_text(
        face = "bold"
      ),
      plot.margin = margin(4, 4, 4, 10)
    )

  save_plot(p, paste0("all_env_heatmap-", pop), w, h)
}

## Main
envs <- c("smoking", "less_than_hs", "more_than_hs", 
          "single", "previously_married", "codeine", 
          "morphine", "ethanol", "nicotine",
          "amphetamines")

populations <- c("BA", "WA", "all")

for (label in envs){
  for (pop in populations){
    plot_tile(pop, label, 10, 6)
    plot_all_env_heatmap(pop)
  }
}

## Reproducibility infor mation
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()