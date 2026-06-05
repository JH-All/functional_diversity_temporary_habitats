# Packages ---------------------
library(vegan)
library(tidyverse)
library(readxl)
library(ade4)
library(mFD)
library(ape)
library(ggdendro)
library(FD)
library(picante)
library(funrar)
library(adiv)
library(cowplot)
library(ggrepel)
library(broom)    
library(rstatix) 
library(car)
library(purrr)
library(tibble)
library(e1071)
library(emmeans)
library(sf)
library(adespatial)
library(spdep)

# Data --------------------
data <- read_excel("community_data.xlsx")
comm <- data[, 6:23]

traits_1 <- read_excel("traits_data.xlsx")

traits <- traits_1 %>%
  column_to_rownames("Species") %>%
  mutate(across(everything(), as.numeric))

# Define trait columns --------------------
m_cols <- colnames(traits)[1:7]
t_col <- "CTMax"

stopifnot(t_col %in% colnames(traits))

# Align community and traits --------------------
spp <- intersect(rownames(traits), colnames(comm))

traits_aln <- as.data.frame(traits[spp, , drop = FALSE])
comm_aln <- as.matrix(comm[, spp, drop = FALSE])

rownames(comm_aln) <- paste0("site", seq_len(nrow(comm_aln)))

stopifnot(identical(rownames(traits_aln), colnames(comm_aln)))

# Environmental data --------------------
env <- read_excel("community_data.xlsx", sheet = "environmental")

env$`maximum_length (m)` <- as.numeric(env$`maximum_length (m)`)
env$`maximum_width (m)`  <- as.numeric(env$`maximum_width (m)`)
env$`maximum_depth (cm)` <- as.numeric(env$`maximum_depth (cm)`)

env$depth_m <- env$`maximum_depth (cm)` / 100

env$Volume <- (2/3) * pi *
  env$`maximum_length (m)` *
  env$`maximum_width (m)` *
  env$depth_m

env$stream_distance <- as.numeric(env$stream_distance)
env$Temp <- as.numeric(env$Temp)
env$DO <- as.numeric(env$DO)
env$pH <- as.numeric(env$pH)

# Functional space choice -----------
traits_pad <- decostand(traits_aln, method = "standardize")
dist <- vegdist(traits_pad, "euclidean")

is.euclid(dist)

quality <- quality.fspaces(
  dist,
  fdendro = "average",
  maxdim_pcoa = 3,
  deviation_weighting = "absolute",
  fdist_scaling = FALSE
)

quality$quality_fspaces

# Functional traits tree ----------------------
hc <- hclust(dist, method = "average")
dg <- dendro_data(hc)

y_max  <- max(segment(dg)$y)
offset <- y_max * 0.06

tree_plot <- ggplot() +
  geom_segment(
    data = segment(dg),
    aes(
      x = y_max - y,
      y = x,
      xend = y_max - yend,
      yend = xend
    ),
    linewidth = 1.6,
    color = "black"
  ) +
  geom_text(
    data = label(dg),
    aes(
      x = y_max + offset,
      y = x,
      label = paste0("italic('", label, "')")
    ),
    parse = TRUE,
    hjust = 0.1,
    size = 5,
    color = "black"
  ) +
  labs(x = "Cophenetic distance", y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(t = 10, r = 80, b = 10, l = 10)
  ) +
  scale_x_continuous(
    limits = c(0, y_max + offset * 4),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(clip = "off")

ggsave("Figure S1.png", plot = tree_plot, width = 10, height = 7, dpi = 300)

# Functional diversity indices -----------------
dend <- hclust(dist, method = "average")
tree_dend <- as.phylo(dend)

FD <- pd(comm_aln, tree_dend)$PD

comm_mat <- comm_aln

labs <- attr(dist, "Labels")
stopifnot(all(labs %in% colnames(comm_mat)))

comm_mat <- comm_mat[, labs, drop = FALSE]
dist_mat <- as.matrix(dist)

comm_use <- comm_mat[rowSums(comm_mat > 0) >= 2, , drop = FALSE]

FD_use <- pd(comm_use, tree_dend)$PD

mpd_values <- mpd(
  samp = comm_use,
  dis = dist_mat,
  abundance.weighted = TRUE
)

mntd_values <- mntd(
  samp = comm_use,
  dis = dist_mat,
  abundance.weighted = TRUE
)

indices <- data.frame(
  Site = rownames(comm_use),
  FD = FD_use,
  MPD = mpd_values,
  MNTD = mntd_values
)

data2 <- data %>%
  mutate(Site = paste0("site", row_number()))

indices <- indices %>%
  left_join(
    data2 %>% select(Site, period, category),
    by = "Site"
  ) %>%
  mutate(
    Period = period,
    category = category,
    habitat = paste0(category, " (", Period, ")"),
    habitat = as.factor(habitat)
  )

levels(indices$habitat) <- gsub("Drier", "Dry", levels(indices$habitat))

# SES Functional Diversity --------------------
indices2 <- indices %>%
  filter(Site %in% rownames(comm_use)) %>%
  mutate(Site = as.character(Site))

set.seed(123)

fd_ses <- ses.pd(
  samp = comm_use,
  tree = tree_dend,
  null.model = "taxa.labels",
  runs = 999,
  include.root = TRUE
)

fd_ses <- fd_ses %>%
  as.data.frame() %>%
  rownames_to_column("Site") %>%
  left_join(indices2 %>% select(Site, habitat), by = "Site") %>%
  rename(
    FD_obs = pd.obs,
    FD_null_mu = pd.rand.mean,
    FD_null_sd = pd.rand.sd,
    SES_FD = pd.obs.z,
    p_quantile = pd.obs.p
  ) %>%
  mutate(
    SES_class = case_when(
      SES_FD >  1.96 ~ "Overdispersed (limiting similarity)",
      SES_FD < -1.96 ~ "Underdispersed (environmental filtering)",
      TRUE ~ "Random expectation"
    )
  )

fig2_A <- ggplot(fd_ses, aes(x = habitat, y = SES_FD, fill = habitat)) +
  geom_hline(
    yintercept = c(-1.96, 0, 1.96),
    linetype = c("dashed", "solid", "dashed"),
    color = "grey45",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  geom_jitter(
    width = 0.12,
    height = 0,
    size = 5,
    alpha = 0.85,
    shape = 21,
    color = "black",
    show.legend = FALSE
  ) +
  labs(
    x = NULL,
    y = "Standard Effect Size (SES)",
    title = "Functional diversity"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = c("#D55E00", "#0072B2", "#E69F00", "#009E73"))

# SES MPD ---------------------------------
set.seed(123)

mpd_ses <- ses.mpd(
  samp = comm_use,
  dis = as.matrix(dist),
  null.model = "taxa.labels",
  runs = 999,
  abundance.weighted = TRUE
)

mpd_ses <- mpd_ses %>%
  as.data.frame() %>%
  rownames_to_column("Site") %>%
  left_join(indices2 %>% select(Site, habitat), by = "Site") %>%
  rename(
    MPD_obs = mpd.obs,
    MPD_null_mu = mpd.rand.mean,
    MPD_null_sd = mpd.rand.sd,
    SES_MPD = mpd.obs.z,
    p_quantile = mpd.obs.p
  ) %>%
  mutate(
    SES_class = case_when(
      SES_MPD >  1.96 ~ "Overdispersed (limiting similarity)",
      SES_MPD < -1.96 ~ "Underdispersed (environmental filtering)",
      TRUE ~ "Random expectation"
    )
  )

fig2_B <- ggplot(mpd_ses, aes(x = habitat, y = SES_MPD, fill = habitat)) +
  geom_hline(
    yintercept = c(-1.96, 0, 1.96),
    linetype = c("dashed", "solid", "dashed"),
    color = "grey45",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  geom_jitter(
    width = 0.12,
    height = 0,
    size = 5,
    alpha = 0.85,
    shape = 21,
    color = "black",
    show.legend = FALSE
  ) +
  labs(
    x = NULL,
    y = "Standard Effect Size (SES)",
    title = "Mean pairwise distance"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = c("#D55E00", "#0072B2", "#E69F00", "#009E73"))

# SES MNTD -----------------------------
set.seed(123)

mntd_ses <- ses.mntd(
  samp = comm_use,
  dis = as.matrix(dist),
  null.model = "taxa.labels",
  runs = 999,
  abundance.weighted = TRUE
)

mntd_ses <- mntd_ses %>%
  as.data.frame() %>%
  rownames_to_column("Site") %>%
  left_join(indices2 %>% select(Site, habitat), by = "Site") %>%
  rename(
    MNTD_obs = mntd.obs,
    MNTD_null_mu = mntd.rand.mean,
    MNTD_null_sd = mntd.rand.sd,
    SES_MNTD = mntd.obs.z,
    p_quantile = mntd.obs.p
  ) %>%
  mutate(
    SES_class = case_when(
      SES_MNTD >  1.96 ~ "Overdispersed (limiting similarity)",
      SES_MNTD < -1.96 ~ "Underdispersed (environmental filtering)",
      TRUE ~ "Random expectation"
    )
  )

fig2_C <- ggplot(mntd_ses, aes(x = habitat, y = SES_MNTD, fill = habitat)) +
  geom_hline(
    yintercept = c(-1.96, 0, 1.96),
    linetype = c("dashed", "solid", "dashed"),
    color = "grey45",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  geom_jitter(
    width = 0.12,
    height = 0,
    size = 5,
    alpha = 0.85,
    shape = 21,
    color = "black",
    show.legend = FALSE
  ) +
  labs(
    x = NULL,
    y = "Standard Effect Size (SES)",
    title = "Mean nearest taxon distance"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = c("#D55E00", "#0072B2", "#E69F00", "#009E73"))

## Figure 2 --------------------
fig2 <- plot_grid(fig2_A, fig2_B, fig2_C, labels = "AUTO", nrow = 1)

ggsave("Figure_2.jpg", fig2, width = 16, height = 6.5)

# Diversity predictors -------------------
env2 <- env %>%
  mutate(Site = paste0("site", row_number()))

df <- indices2 %>%
  left_join(env2, by = "Site") %>%
  mutate(
    habitat = factor(habitat),
    Volume = as.numeric(scale(Volume)),
    NSD = as.numeric(scale(stream_distance)),
    Temperature = as.numeric(scale(Temp)),
    pH = as.numeric(scale(pH))
  ) %>%
  select(Site, habitat, FD, MPD, MNTD, Volume, NSD, Temperature, pH)

responses <- c("FD", "MPD", "MNTD")
preds <- c("Volume", "Temperature", "NSD", "pH")

mods_all <- map(responses, function(resp) {
  df %>%
    group_by(habitat) %>%
    group_map(
      ~ glm(
        reformulate(preds, response = resp),
        data = .x,
        family = Gamma(link = "log")
      ),
      .keep = TRUE
    ) %>%
    set_names(levels(df$habitat))
}) %>%
  set_names(responses)

coef_df <- imap_dfr(mods_all, function(mods_resp, resp_name) {
  imap_dfr(mods_resp, function(m, hab_name) {
    broom::tidy(m) %>%
      filter(term %in% preds) %>%
      mutate(
        habitat = hab_name,
        metric = resp_name,
        sig = p.value < 0.05
      )
  })
})

term_order <- c("pH", "NSD", "Temperature", "Volume")

coef_df <- coef_df %>%
  mutate(
    term = factor(term, levels = term_order),
    metric = factor(metric, levels = c("FD", "MPD", "MNTD"))
  )

cols <- c("#D55E00", "#0072B2", "#E69F00", "#009E73")
names(cols) <- levels(df$habitat)

## Figure 3 --------------------------------
shapes <- c(FD = 21, MPD = 22, MNTD = 23)
pos <- position_dodge(width = 0.65)

fig3 <- ggplot(coef_df, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(
    aes(
      xmin = estimate - 1.96 * std.error,
      xmax = estimate + 1.96 * std.error,
      color = habitat,
      group = metric
    ),
    height = 0.18,
    linewidth = 0.8,
    position = pos,
    show.legend = FALSE
  ) +
  geom_point(
    aes(
      shape = metric,
      color = habitat,
      fill = ifelse(sig, as.character(habitat), "ns"),
      group = metric
    ),
    size = 3.5,
    stroke = 1.1,
    position = pos
  ) +
  facet_wrap(~ habitat, ncol = 2) +
  scale_color_manual(values = cols, guide = "none") +
  scale_fill_manual(values = c(cols, ns = "white"), guide = "none") +
  scale_shape_manual(values = shapes, name = "Metric") +
  scale_x_continuous(limits = c(-9, 9)) +
  labs(
    x = "Model coefficient",
    y = NULL
  ) +
  theme_classic(base_size = 16)+
  theme(
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 1
    )
  )

ggsave("Figure_3.jpg", fig3)

# Convert coordinates from DMS to decimal degrees -------------------------
dms_to_dd <- function(x) {
  x <- str_trim(x)
  nums <- str_extract_all(x, "[0-9]+\\.?[0-9]*")[[1]]
  nums <- as.numeric(nums)
  
  deg <- nums[1]
  min <- ifelse(length(nums) >= 2, nums[2], 0)
  sec <- ifelse(length(nums) >= 3, nums[3], 0)
  
  dd <- deg + min / 60 + sec / 3600
  
  if (str_detect(x, "[SsWw]")) dd <- -dd
  
  dd
}

coord_df <- tibble(coordinates = data$coordinates) %>%
  mutate(
    lat_dms = str_trim(str_extract(coordinates, "^[^ ]+")),
    lon_dms = str_trim(str_extract(coordinates, "[^ ]+$")),
    lat = sapply(lat_dms, dms_to_dd),
    lon = sapply(lon_dms, dms_to_dd)
  ) %>%
  select(lat, lon)

data <- data %>%
  select(-any_of(c("lat", "lon"))) %>%
  bind_cols(coord_df)

# Spatial autocorrelation for GLMs -------------------------
env2 <- env %>%
  mutate(Site = paste0("site", row_number()))

data_sp <- data %>%
  mutate(Site = paste0("site", row_number()))

env2 <- env2 %>%
  left_join(data_sp %>% select(Site, lon, lat), by = "Site")

df_sp <- indices2 %>%
  left_join(env2, by = "Site") %>%
  mutate(
    habitat = factor(habitat),
    Volume = as.numeric(scale(Volume)),
    NSD = as.numeric(scale(stream_distance)),
    Temperature = as.numeric(scale(Temp)),
    pH = as.numeric(scale(pH))
  ) %>%
  select(Site, habitat, FD, MPD, MNTD, Volume, NSD, Temperature, pH, lon, lat)

mods_all_sp <- map(responses, function(resp) {
  df_sp %>%
    group_by(habitat) %>%
    group_map(
      ~ list(
        habitat = unique(.x$habitat),
        response = resp,
        data = .x,
        model = glm(
          reformulate(preds, response = resp),
          data = .x,
          family = Gamma(link = "log")
        )
      ),
      .keep = TRUE
    ) %>%
    set_names(levels(df_sp$habitat))
}) %>%
  set_names(responses)

moran_residuals <- function(mod, dat, k = 4) {
  
  dat <- dat %>% filter(is.finite(lon), is.finite(lat))
  
  if (nrow(dat) < (k + 2)) {
    return(tibble(I = NA_real_, p = NA_real_, n = nrow(dat)))
  }
  
  coords <- as.matrix(dat[, c("lon", "lat")])
  nb <- spdep::knn2nb(spdep::knearneigh(coords, k = k))
  lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  
  r <- residuals(mod, type = "pearson")
  
  mi <- spdep::moran.test(r, lw, zero.policy = TRUE)
  
  tibble(
    I = unname(mi$estimate[["Moran I statistic"]]),
    p = mi$p.value,
    n = length(r)
  )
}

spatial_tests <- imap_dfr(mods_all_sp, function(mods_resp, resp_name) {
  imap_dfr(mods_resp, function(obj, hab_name) {
    moran_residuals(obj$model, obj$data, k = 4) %>%
      mutate(metric = resp_name, habitat = hab_name)
  })
})

spatial_tests

# CWM ----------------------------------
traits_matrix <- as.matrix(traits_aln)

comm_cwm <- comm_use[, rownames(traits_aln), drop = FALSE]

cwm <- functcomp(traits_matrix, comm_cwm)

cwm_df <- as.data.frame(cwm)

stopifnot(nrow(cwm_df) == nrow(comm_use))

cwm_df$Site <- rownames(comm_use)

cwm_df <- cwm_df %>%
  left_join(indices %>% select(Site, habitat), by = "Site")

cwm_num <- cwm_df %>% select(where(is.numeric))
group <- cwm_df$habitat

cwm_num_sc <- scale(cwm_num)

set.seed(123)

perm <- adonis2(
  cwm_num_sc ~ habitat,
  data = cwm_df,
  method = "euclidean",
  permutations = 999
)

perm

bd <- betadisper(
  vegdist(cwm_num_sc, method = "euclidean"),
  group
)

anova(bd)
permutest(bd)

# Variation partitioning of CWM -------------------

## Prepare environmental + spatial data ----------------------------------------

env_varpart <- env %>%
  mutate(Site = paste0("site", row_number())) %>%
  select(Site, Volume, stream_distance, Temp, pH)

coord_varpart <- data %>%
  mutate(Site = paste0("site", row_number())) %>%
  select(Site, lon, lat, category)

varpart_df <- cwm_df %>%
  left_join(env_varpart, by = "Site") %>%
  left_join(coord_varpart, by = "Site")

## Function for variation partitioning -----------------------------------------

run_varpart <- function(dat_group) {
  
  Y <- dat_group %>%
    select(
      Compression_index,
      Relative_depth,
      Index_of_ventral_flattening,
      Relative_eye_position,
      Eye_size,
      Fineness_coefficient,
      Relative_mouth_width,
      CTMax
    ) %>%
    as.data.frame()
  
  Y <- scale(Y)
  
  X_env <- dat_group %>%
    select(Volume, stream_distance, Temp, pH) %>%
    mutate(across(everything(), as.numeric)) %>%
    scale() %>%
    as.data.frame()
  
  pts <- st_as_sf(
    dat_group,
    coords = c("lon", "lat"),
    crs = 4326,
    remove = FALSE
  )
  
  pts_utm <- st_transform(pts, 32723)
  xy <- st_coordinates(pts_utm)
  
  mem <- adespatial::dbmem(xy, silent = TRUE)
  X_space <- as.data.frame(mem)
  
  X_space <- X_space[, apply(X_space, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
  
  vp <- vegan::varpart(Y, X_env, X_space)
  
  ind <- as.data.frame(vp$part$indfract)
  ind$fraction_id <- rownames(ind)
  
  get_adj <- function(pattern) {
    out <- ind %>%
      filter(str_detect(fraction_id, fixed(pattern))) %>%
      pull(Adj.R.squared)
    
    if (length(out) == 0) return(NA_real_)
    out[1]
  }
  
  tibble(
    Environment = get_adj("[a]"),
    Space = get_adj("[b]"),
    Shared = get_adj("[c]"),
    Residual = get_adj("[d]")
  )
}

### 1) All sites together --------------------------------------------------------

varpart_all <- run_varpart(varpart_df) %>%
  mutate(Group = "All sites")

### 2) By habitat type only ------------------------------------------------------

varpart_habitat <- varpart_df %>%
  group_by(category) %>%
  group_modify(~ run_varpart(.x)) %>%
  ungroup() %>%
  rename(Group = category)

### Combine results -------------------------------------------------------------

varpart_results <- bind_rows(
  varpart_all,
  varpart_habitat
)

varpart_results

# Distance-decay  ----------------------------

data <- data %>%
  group_by(lat, lon) %>%
  mutate(
    n_dup = n(),
    id_dup = row_number(),
    angle = runif(n(), 0, 2*pi),
    radius = if_else(
      n_dup > 1,
      (id_dup - 1) * runif(n(), 0.00003, 0.00008),
      0
    ),
    lat_j = lat + radius * sin(angle),
    lon_j = lon + radius * cos(angle)
  ) %>%
  ungroup() %>%
  select(-n_dup, -id_dup, -angle, -radius)

data_sites <- data %>%
  mutate(Site = paste0("site", row_number())) %>%
  filter(Site %in% cwm_df$Site) 

mantel_df <- cwm_df %>%
  left_join(
    data_sites %>% select(Site, category, lon_j, lat_j),
    by = "Site"
  )

## Function ----------------------------------------------------

run_mantel <- function(dat){
  
  Y <- dat %>%
    select(
      Compression_index,
      Relative_depth,
      Index_of_ventral_flattening,
      Relative_eye_position,
      Eye_size,
      Fineness_coefficient,
      Relative_mouth_width,
      CTMax
    )
  
  Yz <- scale(Y)
  
  pts <- st_as_sf(
    dat,
    coords = c("lon_j", "lat_j"),
    crs = 4326
  )
  
  pts_utm <- st_transform(pts, 32723)
  
  xy <- st_coordinates(pts_utm)
  
  dist_space <- as.matrix(dist(xy))
  
  dist_fun <- as.matrix(
    dist(Yz, method = "euclidean")
  )
  
  sim_fun <- 1 / (1 + dist_fun)
  
  dist_fun_decay <- 1 - sim_fun
  
  mantel_res <- mantel(
    dist_fun_decay,
    dist_space,
    method = "spearman",
    permutations = 9999
  )
  
  tibble(
    Mantel_r = mantel_res$statistic,
    p = mantel_res$signif,
    n_sites = nrow(dat)
  )
}


mantel_all <- run_mantel(mantel_df) %>%
  mutate(Group = "All sites")

mantel_rd <- mantel_df %>%
  filter(category == "Road ditches") %>%
  run_mantel() %>%
  mutate(Group = "Road ditches")

mantel_tp <- mantel_df %>%
  filter(category == "Temporary pools") %>%
  run_mantel() %>%
  mutate(Group = "Temporary pools")

mantel_results <- bind_rows(
  mantel_all,
  mantel_rd,
  mantel_tp
)

mantel_results

## Figure 4 --------------------------------------
mantel_df <- cwm_df %>%
  left_join(
    data %>%
      mutate(Site = paste0("site", row_number())) %>%
      select(Site, category, lon, lat),
    by = "Site"
  )

make_decay_df <- function(dat){
  
  Y <- dat %>%
    select(
      Compression_index,
      Relative_depth,
      Index_of_ventral_flattening,
      Relative_eye_position,
      Eye_size,
      Fineness_coefficient,
      Relative_mouth_width,
      CTMax
    )
  
  Yz <- scale(Y)
  
  pts <- st_as_sf(
    dat,
    coords = c("lon", "lat"),
    crs = 4326
  )
  
  pts_utm <- st_transform(pts, 32723)
  xy <- st_coordinates(pts_utm)
  
  dist_space <- as.matrix(dist(xy))
  dist_fun <- as.matrix(dist(Yz))
  
  sim_fun <- 1 / (1 + dist_fun)
  
  data.frame(
    spatial_dist_km = dist_space[upper.tri(dist_space)] / 1000,
    functional_sim = sim_fun[upper.tri(sim_fun)]
  )
}

decay_all <- make_decay_df(mantel_df) %>%
  mutate(Group = "All sites")

decay_rd <- mantel_df %>%
  filter(category == "Road ditches") %>%
  make_decay_df() %>%
  mutate(Group = "Road ditches")

decay_tp <- mantel_df %>%
  filter(category == "Temporary pools") %>%
  make_decay_df() %>%
  mutate(Group = "Temporary pools")

df_decay_all <- bind_rows(decay_all, decay_rd, decay_tp)

fig4 <- ggplot(
  df_decay_all,
  aes(
    x = spatial_dist_km,
    y = functional_sim
  )
) +
  geom_point(
    alpha = 0.25,
    size = 2.8
  ) +
  geom_smooth(
    method = "glm",
    se = TRUE,
    color = "black",
    linewidth = 1
  ) +
  facet_wrap(
    ~ Group,
    nrow = 1,
    scales = "free_x"
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = "Spatial distance (km)",
    y = "Functional similarity (CWM)"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  theme(
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 1
    )
  )

fig4

ggsave(
  "Figure_4.jpg",
  fig4,
  width = 14,
  height = 5,
  dpi = 300
)