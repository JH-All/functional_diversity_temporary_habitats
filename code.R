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
data = read_excel("community_data.xlsx")
comm = data[,6:23]
traits_1 = read_excel("traits_data.xlsx")
traits = traits_1[,2:10]
traits <- traits %>%
  mutate(across(
    .cols = Compression_index:ARR,  
    .fns  = as.numeric
  ))
rownames(traits) = traits_1$Species
env = read_excel("community_data.xlsx", sheet = "environmental")

env$`maximum_length (m)` = as.numeric(env$`maximum_length (m)`)
env$`maximum_width (m)` = as.numeric(env$`maximum_width (m)` )
env$`maximum_depth (cm)` = as.numeric(env$`maximum_depth (cm)`)
env$depth_m = env$`maximum_depth (cm)` / 100
V = (2/3) * pi * env$`maximum_length (m)`* env$`maximum_width (m)`  * env$depth_m 
env$Volume <- V
env$stream_distance = as.numeric(env$stream_distance)
env$Temp = as.numeric(env$Temp)
env$DO = as.numeric(env$DO)
env$pH = as.numeric(env$pH)
class(env$Temp)

# Functional space choice -----------
traits_pad = decostand(traits, method = "standardize")
dist = vegdist(traits_pad, "euclidean")
is.euclid(dist)

?quality.fspaces
quality <- quality.fspaces(dist,fdendro = "average", maxdim_pcoa = 3,
                           deviation_weighting = "absolute",fdist_scaling = FALSE)

quality$quality_fspaces

# Functional traits tree ----------------------
hc <- hclust(dist, method = "average")
dg <- dendro_data(hc)

y_max  <- max(segment(dg)$y)
offset <- y_max * 0.06

tree_plot = ggplot() +
  geom_segment(data = segment(dg),
               aes(x = y_max - y, y = x, xend = y_max - yend, yend = xend),
               linewidth = 1.6, color = "black") +
  geom_text(data = label(dg),
            aes(x = y_max + offset, y = x,
                label = paste0("italic('", label, "')")), 
            parse = TRUE,   
            hjust = 0.1, size = 5, color = "black") + 
  labs(x = "Cophenetic distance", y = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",       
                                  color = "brown4"),
        plot.margin = margin(t = 10, r = 80, b = 10, l = 10)) +
  scale_x_continuous(limits = c(0, y_max + offset*4),
                     expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(clip = "off")

## Figure S1 ---------------------
ggsave("Figure S1.png", plot = tree_plot, width = 10, height = 7, dpi = 300)

# Functional diversity indices -----------------
dend <- hclust(dist, "average")
tree_dend <- as.phylo(dend)
FD <- pd(comm, tree_dend)$PD
FD

comm_mat <- as.matrix(as.data.frame(comm))
if (is.null(rownames(comm_mat))) {
  rownames(comm_mat) <- paste0("site", seq_len(nrow(comm_mat)))
}

labs <- attr(dist, "Labels")
stopifnot(all(labs %in% colnames(comm_mat)))
comm_mat <- comm_mat[, labs, drop = FALSE]
dist_mat <- as.matrix(dist)
comm_use <- comm_mat[rowSums(comm_mat > 0) >= 2, , drop = FALSE]
mpd_values  <- mpd(samp = comm_use, dis = dist_mat, abundance.weighted = TRUE)
mpd_values 
mntd_values <- mntd(samp = comm_use, dis = dist_mat, abundance.weighted = TRUE)
mntd_values 


indices <- data.frame(
  Site = rownames(comm_use),
  FD = FD, 
  MPD  = mpd_values,
  MNTD = mntd_values
)

indices$Period = data$period
indices$category = data$category

indices <- indices %>%
  mutate(habitat = paste0(category, " (", Period, ")"))

indices$habitat = as.factor(indices$habitat)
levels(indices$habitat)
levels(indices$habitat) <- gsub("Drier", "Dry", levels(indices$habitat))
levels(indices$habitat)

# SES Functional Diversity --------------------
stopifnot(all(rownames(comm_use) %in% indices$Site))

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
  left_join(indices2 %>% select(Site, habitat), by = "Site")

fd_ses <- fd_ses %>%
  rename(
    FD_obs     = pd.obs,
    FD_null_mu = pd.rand.mean,
    FD_null_sd = pd.rand.sd,
    SES_FD     = pd.obs.z,
    p_quantile = pd.obs.p
  )


fd_ses <- fd_ses %>%
  mutate(
    SES_class = case_when(
      SES_FD >  1.96 ~ "Overdispersed (limiting similarity)",
      SES_FD < -1.96 ~ "Underdispersed (environmental filtering)",
      TRUE          ~ "Random expectation"
    )
  )


fig2_A = ggplot(fd_ses, aes(x = habitat, y = SES_FD, fill = habitat)) +
  geom_hline(
    yintercept = c(-1.96, 0, 1.96),
    linetype = c("dashed", "solid", "dashed"),
    color = "grey45", linewidth = 0.8,
    show.legend = FALSE
  ) +
  geom_jitter(
    width = 0.12, height = 0, size = 5, alpha = 0.85,
    shape = 21, color = "black", show.legend = FALSE
  ) +
  labs(
    x = NULL,
    y = "Standard Effect Size (SES)",
    title = "Functional diversity"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = c("#D55E00", "#0072B2",
                               "#E69F00", "#009E73"))

fig2_A

# SES MPD ---------------------------------
stopifnot(all(rownames(comm_use) %in% indices$Site))

indices2 <- indices %>%
  filter(Site %in% rownames(comm_use)) %>%
  mutate(Site = as.character(Site))

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
    MPD_obs     = mpd.obs,
    MPD_null_mu = mpd.rand.mean,
    MPD_null_sd = mpd.rand.sd,
    SES_MPD     = mpd.obs.z,
    p_quantile  = mpd.obs.p
  ) %>%
  mutate(
    SES_class = case_when(
      SES_MPD >  1.96 ~ "Overdispersed (limiting similarity)",
      SES_MPD < -1.96 ~ "Underdispersed (environmental filtering)",
      TRUE           ~ "Random expectation"
    )
  )

fig2_B = ggplot(mpd_ses, aes(x = habitat, y = SES_MPD, fill = habitat)) +
  geom_hline(
    yintercept = c(-1.96, 0, 1.96),
    linetype = c("dashed", "solid", "dashed"),
    color = "grey45", linewidth = 0.8,
    show.legend = FALSE
  ) +
  geom_jitter(
    width = 0.12, height = 0, size = 5, alpha = 0.85,
    shape = 21, color = "black", show.legend = FALSE
  ) +
  labs(
    x = NULL,
    y = "Standard Effect Size (SES)",
    title = "Mean pairwise distance"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  scale_fill_manual(values = c("#D55E00", "#0072B2", "#E69F00", "#009E73"))

fig2_B 

# SES MNTD -----------------------------
stopifnot(all(rownames(comm_use) %in% indices$Site))

indices2 <- indices %>%
  filter(Site %in% rownames(comm_use)) %>%
  mutate(Site = as.character(Site))

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
    MNTD_obs     = mntd.obs,
    MNTD_null_mu = mntd.rand.mean,
    MNTD_null_sd = mntd.rand.sd,
    SES_MNTD     = mntd.obs.z,
    p_quantile   = mntd.obs.p
  ) %>%
  mutate(
    SES_class = case_when(
      SES_MNTD >  1.96 ~ "Overdispersed (limiting similarity)",
      SES_MNTD < -1.96 ~ "Underdispersed (environmental filtering)",
      TRUE            ~ "Random expectation"
    )
  )

fig2_C = ggplot(mntd_ses, aes(x = habitat, y = SES_MNTD, fill = habitat)) +
  geom_hline(
    yintercept = c(-1.96, 0, 1.96),
    linetype = c("dashed", "solid", "dashed"),
    color = "grey45", linewidth = 0.8,
    show.legend = FALSE
  ) +
  geom_jitter(
    width = 0.12, height = 0, size = 5, alpha = 0.85,
    shape = 21, color = "black", show.legend = FALSE
  ) +
  labs(
    x = NULL,
    y = "Standard Effect Size (SES)",
    title = "Mean nearest taxon distance"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = c("#D55E00", "#0072B2", "#E69F00", "#009E73"))

fig2_C

## Figure 2 --------------------
fig2 = plot_grid(fig2_A, fig2_B, fig2_C, labels = "AUTO", nrow=3)
fig2
ggsave("Figure_2.jpg", fig2, width = 6, height = 16)

# Diversity predictors -------------------
env2 <- env %>%
  mutate(Site = paste0("site", row_number()))
indices2 <- indices %>% mutate(Site = as.character(Site))

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
    group_map(~ glm(
      reformulate(preds, response = resp),
      data = .x,
      family = Gamma(link = "log")
    ), .keep = TRUE) %>%
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

fig3 = ggplot(coef_df, aes(x = estimate, y = term)) +
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
    size = 5,
    stroke = 1.1,
    position = pos
  ) +
  
  facet_wrap(~ habitat, ncol = 2) +
  scale_color_manual(values = cols, guide = "none") +
  scale_fill_manual(values = c(cols, ns = "white"), guide = "none") +
  scale_shape_manual(values = shapes, name = "Metric") +
  scale_x_continuous(limits = c(-9, 9)) +
  labs(
    x = "Standardized coefficient (Gamma GLM)",
    y = NULL
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave("Figure_3.jpg", fig3)

## Spatial autocorrelation for GLMs -------------------------
env2  <- env  %>% mutate(Site = paste0("site", row_number()))
data2 <- data %>% mutate(Site = paste0("site", row_number()))

env2 <- env2 %>%
  left_join(data2 %>% select(Site, lon, lat), by = "Site")

indices2 <- indices %>% mutate(Site = as.character(Site))

df <- indices2 %>%
  left_join(env2, by = "Site") %>%
  mutate(
    habitat = factor(habitat),
    Volume = as.numeric(scale(Volume)),
    NSD = as.numeric(scale(stream_distance)),
    Temperature = as.numeric(scale(Temp)),
    pH = as.numeric(scale(pH))
  ) %>%
  select(Site, habitat, FD, MPD, MNTD, Volume, NSD, Temperature, pH, lon, lat)

responses <- c("FD", "MPD", "MNTD")
preds <- c("Volume", "Temperature", "NSD", "pH")

mods_all <- map(responses, function(resp) {
  df %>%
    group_by(habitat) %>%
    group_map(~ list(
      habitat = unique(.x$habitat),
      response = resp,
      data = .x,
      model = glm(
        reformulate(preds, response = resp),
        data = .x,
        family = Gamma(link = "log")
      )
    ), .keep = TRUE) %>%
    set_names(levels(df$habitat))
}) %>%
  set_names(responses)

moran_residuals <- function(mod, dat, k = 4){
  dat <- dat %>% filter(is.finite(lon), is.finite(lat))
  if(nrow(dat) < (k + 2)) return(tibble(I = NA_real_, p = NA_real_, n = nrow(dat)))
  
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

spatial_tests <- imap_dfr(mods_all, function(mods_resp, resp_name){
  imap_dfr(mods_resp, function(obj, hab_name){
    out <- moran_residuals(obj$model, obj$data, k = 4)
    out %>% mutate(metric = resp_name, habitat = hab_name)
  })
})

spatial_tests

# CWM ----------------------------------
traits_matrix = as.matrix(traits)
cwm <- functcomp(traits_matrix, as.matrix(comm_use))
cwm

cwm_df <- as.data.frame(cwm)
stopifnot(nrow(cwm_df) == length(indices$habitat))
cwm_df$habitat <- indices$habitat

cwm_num <- cwm_df %>% select(where(is.numeric))   
group    <- cwm_df$habitat                        

cwm_num_sc <- scale(cwm_num)                   

set.seed(123)
perm <- adonis2(cwm_num_sc ~ habitat, data = cwm_df, method = "euclidean", permutations = 999)
perm

bd   <- betadisper(vegdist(cwm_num_sc, method = "euclidean"), group)
anova(bd)           
permutest(bd)   

## Figure 4 ----------------------------
cwm_long <- cwm_df %>%
  pivot_longer(
    cols = -habitat,          
    names_to = "Trait", 
    values_to = "Value"
  )
cwm_long$Trait = as.factor(cwm_long$Trait)
levels(cwm_long$Trait) <- gsub("_", " ", levels(cwm_long$Trait))
levels(cwm_long$Trait)

str(cwm_long$Trait)

levels(cwm_long$Trait)[levels(cwm_long$Trait) == "Eye size"] <- "Relative eye size"

cwm_long$Trait <- factor(
  cwm_long$Trait,
  levels = c(
    "Compression index",
    "Relative depth",
    "Index of ventral flattening",
    "Relative eye position",
    "Relative eye size",
    "Fineness coefficient",
    "Relative mouth width",
    "CTMax",
    "ARR"
  )
)

fig4 = ggplot(cwm_long, aes(x = habitat, y = Value, fill = habitat)) +
  geom_boxplot(show.legend = F, width = 0.7, alpha = 0.6, size = 1.3)+
  geom_jitter(width = 0.15, shape = 21, color = "black",
              size = 3, show.legend = F, stroke =1.2, alpha = 0.7)+
  facet_wrap(~ Trait, scales = "free_x") +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Community Weighted Mean")+
  scale_fill_manual(values = c("#D55E00", "#0072B2",
                               "#E69F00", "#009E73"))

fig4

ggsave("Figure_4.jpg", plot = fig4, width = 15, height = 12, 
       dpi = 300)

# Distance decay ------------------------
dms_to_dd <- function(x) {
  x <- str_trim(x)
  nums <- str_extract_all(x, "[0-9]+\\.?[0-9]*")[[1]]
  nums <- as.numeric(nums)
  
  deg <- nums[1]
  min <- ifelse(length(nums) >= 2, nums[2], 0)
  sec <- ifelse(length(nums) >= 3, nums[3], 0)
  
  dd <- deg + min/60 + sec/3600
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
  mutate(
    lat = coord_df$lat,
    lon = coord_df$lon
  )

set.seed(123)

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


Y <- cwm_df %>%
  select(
    Compression_index,
    Relative_depth,
    Index_of_ventral_flattening,
    Relative_eye_position,
    Eye_size,
    Fineness_coefficient,
    Relative_mouth_width,
    CTMax,
    ARR
  )

names(Y)[names(Y) == "Eye_size"] <- "Relative_eye_size"

Yz <- scale(Y)

coords <- data %>%
  select(lon_j, lat_j)

pts <- st_as_sf(coords, coords = c("lon_j", "lat_j"), crs = 4326)
pts_utm <- st_transform(pts, 32723)

xy <- st_coordinates(pts_utm)

dist_space <- as.matrix(dist(xy))
dist_fun   <- as.matrix(dist(Yz, method = "euclidean"))

sim_fun <- 1 / (1 + dist_fun)

get_upper <- function(mat) mat[upper.tri(mat)]

df_decay <- data.frame(
  spatial_dist_km = get_upper(dist_space) / 1000,
  functional_sim  = get_upper(sim_fun)
)

dist_fun <- 1 - sim_fun

mantel_res <- mantel(
  dist_fun,    
  dist_space,  
  method = "spearman", 
  permutations = 9999
)

mantel_res

## Figure 5 -------------------------------
fig5 = ggplot(df_decay, aes(x = spatial_dist_km, y = functional_sim)) +
  geom_point(alpha = 0.3, size = 3) +
  geom_smooth(
    method = "glm",
    se = TRUE,
    color = "black"
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = "Spatial distance (km)",
    y = "Functional similarity (CWM)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25))

ggsave("Figure_5.jpg", fig5)