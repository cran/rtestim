## ----setup, message = FALSE, warning = FALSE, collapse=TRUE-------------------
library(rtestim)
library(ggplot2)
library(dplyr)
library(nnet)
library(forcats)
library(tidyr)
theme_set(theme_bw())

## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  dpi = 300,
  collapse = FALSE,
  comment = "#>",
  fig.asp = 0.618,
  fig.width = 6,
  out.width = "80%",
  fig.align = "center"
)

## ----gamma-pdf, echo = FALSE--------------------------------------------------
ggplot(data.frame(x = seq(.001, 20, length.out = 100)), aes(x)) +
  stat_function(fun = dgamma, args = list(shape = 2.5, scale = 2.5), color = "orange") +
  stat_function(fun = dgamma, args = list(shape = 1, scale = 2), color = "cornflowerblue") +
  stat_function(fun = dgamma, args = list(shape = 2, scale = 1), color = "forestgreen") +
  scale_y_continuous(expand = expansion(c(0, 0.05)))

## ----can-default--------------------------------------------------------------
can_default <- estimate_rt(cancovid$incident_cases, x = cancovid$date, nsol = 20L)
plot(can_default) + coord_cartesian(ylim = c(0.5, 2))

## ----cancov-nonpar------------------------------------------------------------
# Data from Backer et al.
delay <- read.csv("backer.csv") |>
  filter(delay > 0) |>
  select(-delay)
delay <- rowSums(delay)
delay <- delay / sum(delay)

## ----cancov-nonpar-plot, echo=FALSE-------------------------------------------
ggplot(
  data.frame(delay = 1:length(delay), probability = delay), 
  aes(delay, probability)
) +
  geom_point(colour = "cornflowerblue") +
  geom_segment(yend = 0, colour = "cornflowerblue") +
  scale_y_continuous(expand = expansion(c(0, 0.05)))

## ----can-nonpar---------------------------------------------------------------
can_nonpar <- estimate_rt(
  cancovid$incident_cases, 
  x = cancovid$date, 
  delay_distn = delay,
  nsol = 20L)
plot(can_nonpar) + coord_cartesian(ylim = c(0.5, 2))

## ----cdfs, echo = FALSE-------------------------------------------------------
cdfs <- data.frame(
  x = 0:20,
  default = pgamma(0:20, shape = 2.5, scale = 2.5),
  backer = cumsum(c(0, delay, rep(0, 20 - length(delay))))
)
ggplot(cdfs |> pivot_longer(-x), aes(x, value, color = name)) +
  geom_step() +
  scale_y_continuous(expand = expansion(0)) +
  scale_color_manual(values = c("orange", "cornflowerblue"), name = "")

## ----backer-matrix, eval=FALSE------------------------------------------------
# # library(Matrix)
# n <- nrow(cancovid)
# backer_delay <- c(0, delay, rep(0, n - length(delay) - 1))
# delay_mat <- matrix(0, n, n)
# delay_mat[1,1] <- 1
# for (iter in 2:n) delay_mat[iter, 1:iter] <- rev(backer_delay[1:iter])
# delay_mat <- drop0(as(delay_mat, "CsparseMatrix")) # make it sparse, not necessary
# delay_mat <- delay_mat / rowSums(delay_mat) # renormalize

## ----duotang-processing, eval=FALSE, echo=TRUE--------------------------------
# # Run on 19 April 2024
# duotang <- read_tsv("https://github.com/CoVaRR-NET/duotang/raw/main/data_needed/virusseq.metadata.csv.gz")
# columnlist <- c(
#   "fasta_header_name", "province", "host_gender", "host_age_bin",
#   "sample_collection_date", "sample_collected_by",
#   "purpose_of_sampling", "purpose_of_sequencing", "lineage",
#   "raw_lineage", "gisaid_accession", "isolate"
# )
# unknown.str <- c(
#   "Undeclared", "Not Provided", "Restricted Access", "Missing",
#   "Not Applicable", "", "NA", "unknow"
# )
# duotang <- duotang |>
#   rename(province = geo_loc_name_state_province_territory) |>
#   select(all_of(columnlist))
# 
# meta <- duotang |>
#   mutate(
#     week = cut(sample_collection_date, "week"),
#     month = gsub("-..$", "", as.character(cut(sample_collection_date, "month")))
#   )
# source("https://github.com/CoVaRR-NET/duotang/raw/main/scripts/scanlineages.R")
# meta <- meta |>
#   mutate(gisaid_accession = str_replace(gisaid_accession, "EPI_ISL_", "")) |>
#   rename(GID = gisaid_accession) |>
#   rowwise() |>
#   mutate(raw_lineage = ifelse(
#     grepl("^X", raw_lineage),
#     str_replace_all(paste0(
#       realtorawlineage(substr(
#         raw_lineage, 1, str_locate(raw_lineage, "\\.") - 1
#       )),
#       ".",
#       substr(raw_lineage, str_locate(raw_lineage, "\\.") + 1, nchar(raw_lineage))
#     ), "[\r\n]", ""),
#     raw_lineage
#   )) |>
#   ungroup()
# dico <- makepangolindico() # rebuild the lineage dictionary so the correct names gets assigned for XBB descedants not named XBB
# 
# VOCVOI <- read_csv("https://raw.githubusercontent.com/CoVaRR-NET/duotang/main/resources/vocvoi.csv")
# meta$pango_group <- create.pango.group(VOCVOI, meta)
# meta <- select(meta, province, week, pango_group) |>
#   mutate(week = as.Date(week))
# 
# counts <- group_by(meta, province, week, pango_group) |>
#   count() |>
#   ungroup() |>
#   arrange(province, week, pango_group)
# can_counts <- group_by(meta, week, pango_group) |>
#   count() |>
#   ungroup() |>
#   arrange(week, pango_group) |>
#   mutate(province = "Canada")
# counts <- bind_rows(can_counts, counts)
# saveRDS(counts, "duotang-counts.rds")

## ----duotang-counts, echo=FALSE, message=FALSE--------------------------------
props <- readRDS("duotang-counts.rds") |>
  pivot_wider(names_from = pango_group, values_from = n, values_fill = 0) |>
  mutate(total = rowSums(across(-c(week, province)))) |>
  mutate(across(-c(week, province, total), ~ .x / total)) |>
  select(-total)

smooth_it <- function(props_group) {
  z <- props_group |> select(-week, -province)
  n <- names(z)
  nn <- gsub(" ", "_", n)
  names(z) <- nn
  form_resp <- paste0("cbind(", paste0(names(z), collapse = ",") ,") ~ ")
  z$time <- as.numeric(props_group$week)
  form <- as.formula(paste0(form_resp, "poly(time, degree = 3)"))
  fits <- multinom(form, z, trace = FALSE)
  rng <- range(props_group$week)
  alltime <- as.numeric(seq(rng[1], rng[2], by = 1))
  z <- as_tibble(predict(fits, data.frame(time = alltime), type = "probs")) |>
    mutate(Date = as.Date(alltime))
  z
}

can_props_smoothed <- smooth_it(props |> filter(province == "Canada"))

can_props_smoothed |>
  pivot_longer(-Date) |>
  ggplot(aes(Date, y = value, fill = name)) +
  geom_area(position = "stack") +
  ylab("Variants in circulation") +
  xlab("") + 
  theme_bw() +
  scale_x_date(name = "", date_breaks = "1 year", date_labels = "%Y", expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_viridis_d(name = "")

can_pred_class <- can_props_smoothed |>
  filter(Date <= max(cancovid$date)) |>
  pivot_longer(-Date) |> 
  group_by(Date) |> 
  summarise(var = name[which.max(value)], .groups = "drop") |>
  mutate(var = case_when(
    var == "other" ~ "Ancestral lineage",
    var == "Alpha" ~ "Alpha",
    var == "Beta" ~ "Beta",
    var == "Delta" ~ "Delta",
    TRUE ~ "Omicron"
  ))

boot <- data.frame(
  Date = seq(min(cancovid$date), min(can_props_smoothed$Date) - 1, by = "day"),
  var = "Ancestral lineage"
)
can_pred_class <- bind_rows(boot, can_pred_class)

## ----echo=FALSE, fig.height=1, dev='png'--------------------------------------
ggplot(can_pred_class, aes(x = Date, fill = fct_relevel(var, "Ancestral lineage"))) + 
  geom_ribbon(aes(ymax = 1, ymin = 0)) +
  scale_y_continuous(expand = expansion(0)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_fill_brewer(palette = "Dark2", name = "")

## ----meta-analysis, eval=FALSE------------------------------------------------
# data_raw <- readxl::read_excel("xu-etal-DATA_RAW.xlsx") |>
#   select(type, para, n = Sample_size, mean, sd, se, median) |>
#   filter(!is.na(type)) |>
#   mutate(across(-c(type, para), as.numeric))
# 
# bonehead_meta <- data_raw |>
#   group_by(type, para) |>
#   mutate(
#     no_n = all(is.na(n)),
#     n = case_when(!is.na(n) ~ n, no_n ~ 1, TRUE ~ median(n, na.rm = TRUE))
#   ) |>
#   ungroup() |>
#   mutate(
#     mean = case_when(!is.na(mean) ~ mean, TRUE ~ median),
#     sd = case_when(!is.na(sd) ~ sd, TRUE ~ se * sqrt(n))
#   ) |>
#   group_by(type, para) |>
#   summarise(
#     mean = median(mean, na.rm = TRUE),
#     sd = median(sd, na.rm = TRUE),
#     .groups = "drop"
#   )
# ## There's only one Beta and only IP.
# ## We use the corresponding sd for Alpha IP,
# ## and duplicate Alpha for GT / ST
# 
# Beta_IP <- bonehead_meta |> filter(type == "Beta")
# Beta_IP$sd = bonehead_meta |>
#   filter(type == "Alpha", para == "IP") |> pull(sd)
# Beta <- bind_rows(
#   Beta_IP,
#   bonehead_meta |>
#     filter(type == "Alpha", para != "IP") |>
#     mutate(type = "Beta")
# )
# 
# delay_dstns_byvar <- bonehead_meta |>
#   filter(type != "Beta") |>
#   bind_rows(Beta) |>
#   arrange(type, para) |>
#   mutate(shape = mean^2 / sd^2, scale = mean / shape)
# saveRDS(delay_dstns_byvar, "delay-distns-byvar.rds")

## ----variant-si-plot, echo=FALSE----------------------------------------------
delay_dstns_can <- readRDS("delay-distns-byvar.rds") |>
  filter(para == "SI", type %in% unique(can_pred_class$var)) 

delay_dstns_can |>
  rowwise() |>
  mutate(probability = list(discretize_gamma(0:20, shape, scale))) |>
  ungroup() |>
  select(type, probability) |>
  unnest(probability) |>
  group_by(type) |>
  mutate(delay = row_number() - 1) |>
  ggplot(aes(delay, probability, colour = type)) +
  geom_line() + 
  scale_color_brewer(palette = "Dark2", name = "") +
  scale_y_continuous(expand = expansion(c(0, 0.05)))

## ----tvar-matrix, message=FALSE-----------------------------------------------
library(Matrix)
n <- nrow(cancovid)
delay_mat <- matrix(0, n, n)
delay_mat[1,1] <- 1
for (iter in 2:n) {
  current_var <- can_pred_class$var[iter]
  current_pars <- delay_dstns_can |> filter(type == current_var)
  delay <- discretize_gamma(0:(iter - 1), current_pars$shape, current_pars$scale)
  delay_mat[iter, 1:iter] <- rev(delay)
}
delay_mat <- drop0(as(delay_mat, "CsparseMatrix")) # make it sparse, not necessary
delay_mat <- delay_mat / rowSums(delay_mat) # renormalize

## ----tvar-Rt------------------------------------------------------------------
can_tvar <- estimate_rt(
  cancovid$incident_cases, 
  x = cancovid$date, 
  delay_distn = delay_mat,
  nsol = 20L)
plot(can_tvar) + coord_cartesian(ylim = c(0.5, 2))

