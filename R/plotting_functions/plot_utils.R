my_paletteS <- list(
  "red" = "#D55E00",
  "blue" = "#0072B2", 
  "green" = "#009E73",
  "yellow" = "#E69F00",
  "pink" = "#CC79A7",
  "light_blue" = "#56B4E9",
  "grey" = "#999999",
  "background" = "#FAFAFA"
)

my_paletteCB <- list(
  "red" = '#e41a1c',
  "blue" = '#377eb8',
  "green" = '#4daf4a',
  "violet" = '#984ea3',
  "orange" = '#ff7f00',
  "yellow" = '#ffff33',
  "brown" = '#a65628',
  "pink" = '#f781bf',
  "grey" = '#999999'
)

my_paletteN2 <- list(
  "blue"   = "#00798c",
  "red"    = "#d1495b",
  "green"  = "#a2d729",
  "yellow" = "#edae49",
  "navy"   = "#2e4057", 
  "grey"   = "#8d96a3"
)

my_paletteOH <- list(
  "blue" = rgb(57/255,106/255,177/255), 
  "green" = "#69b3a2"
)

my_palette <- my_paletteS

my_palette_methods0 <- c("Uncond"=my_paletteCB$orange,
                         "Semi_cond"=my_paletteCB$violet,
                         "GRF"=my_paletteCB$green,
                         "gpd_GAM"=my_paletteOH$green,
                         "gbex"=my_paletteCB$red,
                         "EQRNN"=my_paletteCB$blue,
                         "EQRNN2"="#25567d",
                         "EQRNN1"="#377eb8")

my_palette_methodsS <- list(
  c("method" = "EQRN", "color" = my_palette$blue),
  c("method" = "EQRN2", "color" = "#004166"),
  c("method" = "EQRNN", "color" = my_palette$blue),
  c("method" = "EQRNN1", "color" = my_palette$blue),
  c("method" = "EQRNN2", "color" = "#004166"),
  c("method" = "GBEX", "color" = my_palette$red),
  c("method" = "gbex", "color" = my_palette$red),
  c("method" = "EGAM", "color" = my_palette$green),
  c("method" = "EVGAM", "color" = my_palette$green),
  c("method" = "gpd_GAM", "color" = my_palette$green),
  c("method" = "EXQAR", "color" = my_palette$light_blue),
  c("method" = "GRF", "color" = my_palette$green),
  c("method" = "Semi-cond", "color" = my_palette$pink),
  c("method" = "Semicond", "color" = my_palette$pink),
  c("method" = "Semi_cond", "color" = my_palette$pink),
  c("method" = "Uncond", "color" = my_palette$yellow),
  c("method" = "Other", "color" = my_palette$light_blue)
) %>% 
  purrr::transpose() %>% 
  tibble::as_tibble() %>% 
  tidyr::unnest(cols = c(method, color)) %>% tibble::deframe()

my_palette_methods <- my_palette_methodsS

theme_fct <- function(font_size=11, font_size_axes=10, font_size_captions=7.5){
  ggplot2::theme_set(ggplot2::theme_bw() +
                       ggplot2::theme(
                         plot.background = element_blank(),
                         panel.background = element_blank(),
                         legend.background = element_blank(),
                         strip.background = element_blank(),
                         strip.placement = "outside",
                         strip.text = element_text(size = font_size_axes),
                         strip.switch.pad.grid = unit(0, "pt"),
                         strip.switch.pad.wrap = unit(0, "pt"),
                         plot.caption=element_text(size=font_size_captions, hjust=0, 
                                                   margin=margin(t=15)),
                         text = element_text(size = font_size),
                         axis.ticks = element_blank(),
                         axis.text = element_text(size = font_size_axes),#
                         panel.grid.major = element_line(size = 0.25),
                         legend.position = "bottom",
                         legend.box.spacing = margin(0.5),
                         plot.margin=margin(0,0,0,0, unit="pt")
                       )
  )
}

theme_fct(font_size=11, font_size_axes=10, font_size_captions=7.5)

save_myplot <- function(plt, plt_nm, width, height, 
                        width_pdf = 1270, height_pdf = 1270,
                        crop = TRUE, cairo = FALSE, wh_units="mm", ...) {
  
  dir_name <- dirname(plt_nm)
  if (!file.exists(dir_name)){
    dir.create(dir_name)
  }
  
  if (cairo) {
    ggplot2::ggsave(plt_nm, 
                    egg::set_panel_size(p = plt, 
                                        width = unit(width, wh_units), 
                                        height = unit(height, wh_units)),
                    width = width_pdf, height = height_pdf,
                    limitsize = FALSE, units = c("mm"), 
                    device = cairo_pdf, family = "Arial", ...)
  } else {
    ggplot2::ggsave(plt_nm, 
                    egg::set_panel_size(p = plt, 
                                        width = unit(width, wh_units), 
                                        height = unit(height, wh_units)),
                    width = width_pdf, height = height_pdf,
                    limitsize = FALSE, units = c("mm"), ...)
  }
  
  if (crop){
    knitr::plot_crop(plt_nm)
  } 
}
