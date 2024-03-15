library(hexSticker)

library(showtext)
## Loading Google fonts (http://www.google.com/fonts)
font_add_google("Outfit")

# ### sticker
# sticker("inst/figures/abn_logo.png",
#         s_x = 1,
#         s_y = 0.85,
#         s_width = 0.6,
#         s_height = 0.6,
#         package = "abn",
#         p_x = 1,
#         p_y = 1.5,
#         p_size = 28,
#         p_family = "Outfit",
#         p_color = "#2A5783",
#         h_fill = "#D0DFF3",
#         h_color = "#A3C1E8",
#         # url = "http://r-bayesian-networks.org/",
#         filename = "inst/figures/logo.png")
#
# ### sticker Purple text top
# sticker("inst/figures/abn_logo_white.png",
#         s_x = 1,
#         s_y = 0.85,
#         s_width = 0.6,
#         s_height = 0.6,
#         package = "abn",
#         p_x = 1,
#         p_y = 1.5,
#         p_size = 36,
#         p_family = "Outfit",
#         p_color = "#ffffff",
#         h_fill = "#6f1c54ff",
#         h_color = "#6f1c54ff",
#         # url = "http://r-bayesian-networks.org/",
#         filename = "inst/figures/logo.png")

### sticker purple text bottom
sticker("inst/figures/abn_logo_white.png",
        s_x = 1,
        s_y = 1.2,
        s_width = 0.6,
        s_height = 0.6,
        package = "abn",
        p_x = 1,
        p_y = 0.5,
        p_size = 36,
        p_family = "Outfit",
        p_color = "#ffffff",
        h_fill = "#6f1c54ff",
        h_color = "#6f1c54ff",
        filename = "inst/figures/logo.png")

### add sticker to package
usethis::use_logo("inst/figures/logo.png")
