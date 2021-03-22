# figure funs: theming

theme_proj <- 
  cowplot::theme_minimal_hgrid(font_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
theme_map <- 
  cowplot::theme_map(font_size = 12)
