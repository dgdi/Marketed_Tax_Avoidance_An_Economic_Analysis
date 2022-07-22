#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@           Marketed Tax Avoidance:            @@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@             An Economic Analysis             @@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@                   Figure 1                   @@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# This script has been tested on Win 10, R-4.1.1, RStudio-1.4.1717
# the list of the required packages is at the start of the code.
#
# These files are distributed under BSD license. For more information, please check the license.txt file.
#
# Reference article:
# Li, J., Gamannossi degl'Innocenti, D. and Rablen, M.D., 2021. Marketed Tax Avoidance Schemes: An Economic Analysis. Available at SSRN 3969577.
#
# For any question, suggestion or comment,
# write to:
# mail@dgdi.me

#cleaning environment
rm(list = ls())

#loading libraries
library("here")
library("data.table")
library("ggplot2")
library("extrafont")
library("tikzDevice")

#sourcing general parameters
source(here("R/par.R"))

#sourcing functions
source(here("R/fun.R"))

#defining parameters
rho = .24
c0 = .24
c1 = 0
gamma <-
    5 #this is unnecessary theoretically but needs to be specified to avoid numerical errors
p = 0.7
tau = 2
wmax = 16
wmin = 9

#identify numerically the W0 (no analytical solution available, see Mathematica file fig1_expr)
W0_root_a <- uniroot(
    f = W0_root_fun,
    interval = c(10, 15),
    c0 = c0,
    c1 = c1,
    gamma = gamma,
    rho = rho,
    F_ = tau,
    check.conv = check_conv,
    tol = tol
)$root

#identify analytically W1 (see the Mathematica file analytic.nb)
W1_root_a <- tau / (p * c0)

#identify analytically W2 (see the Mathematica file analytic.nb)
W2_root_a <- (((1 - p) * tau) / ((1 - c0) * (1 - p - rho)))

#defining a vector with all critical W
Ws_a <- c(W0_root_a,
          W1_root_a,
          W2_root_a)

#creating vector of W
W_vec_a <- sort(c(seq(
    from = wmin,
    to = wmax,
    length.out = 100
),
Ws_a))

#creating a data.table to store A by W
A_by_W_a_dt <- data.table(W =  W_vec_a)
#initially set all A to 0
A_by_W_a_dt[, A := 0]
#set A = T just above A W0 (allows for the correct representation using fewer points - better for tikz)
A_by_W_a_dt <- rbindlist(list(A_by_W_a_dt,
                              data.table(
                                  W = W0_root_a + .000001,
                                  A = T_fun(
                                      W = W0_root_a + .000001,
                                      c0 = c0,
                                      c1 = c1,
                                      gamma = gamma
                                  )
                              )))
#set A = T just for W1 > W > W0
A_by_W_a_dt[W > W0_root_a &
                W < W1_root_a, A := T_fun(W = W,
                                          c0 = c0,
                                          c1 = c1,
                                          gamma = gamma)]
#set A = tau/p for W2 > W >= W1
A_by_W_a_dt[W >= W1_root_a &
                W < W2_root_a, A := tau / p]

#set A = A* for W >= W2
A_by_W_a_dt[W >= W2_root_a, A := A_fun(
    W = W,
    c0 = c0,
    c1 = c1,
    gamma = gamma,
    p = p,
    rho = rho
)]

#### Plots ###############################################################
#define variable to store errors in tikz compilation
tmptikz_err <- NULL

#### Figure_1(a) ###############################################################
plotname <- "Figure_1(a)"
tmptikz <- paste0(here(), "/tikz/tmp_", plotname)
setTikzDefaults()
options(
    tikzLwdUnit = 77 / 96,
    tikzMetricPackages = c(
        "\\usepackage[utf8]{inputenc}",
        "\\usepackage[T1]{fontenc}",
        "\\usetikzlibrary{calc}",
        "\\usepackage{amssymb}"
    ),
    tikzMetricsDictionary = tmptikz
)
tikz(
    paste0(here("tikz/"), "/", plotname, ".tex"),
    width = 7,
    height = 7,
    standAlone = TRUE,
    
    packages = c(
        "\\usepackage{tikz}",
        "\\usepackage[active,tightpage,psfixbb]{preview}",
        "\\PreviewEnvironment{pgfpicture}",
        "\\setlength\\PreviewBorder{0pt}",
        "\\usepackage{amssymb}",
        "\\usepackage{amsfonts}",
        "\\usepackage{eurosym}",
        "\\usepackage{amsmath}",
        "\\usepackage{amssymb}",
        "\\usepackage[T1]{fontenc}"
    )
)

ggplot(data = A_by_W_a_dt) +
    
    geom_segment(
        aes(
            x = wmin,
            y = tau / p,
            xend = W1_root_a,
            yend = tau / p
        ),
        linetype = "dotted",
        colour = "gray",
        size = .3
    )  +
    
    geom_segment(
        aes(
            x = W1_root_a,
            y = 0,
            xend = W1_root_a,
            yend = tau / p
        ),
        linetype = "dotted",
        colour = "gray",
        size = .3
    )  +
    
    geom_segment(
        aes(
            x = W2_root_a,
            y = 0,
            xend = W2_root_a,
            yend = tau / p
        ),
        linetype = "dotted",
        colour = "gray",
        size = .3
    )  +
    
    
    geom_segment(
        aes(
            x = wmin,
            y = T_fun(
                W = wmin,
                c0 = c0,
                c1 = c1,
                gamma = gamma
            ),
            xend = wmax,
            yend = T_fun(
                W = wmax,
                c0 = c0,
                c1 = c1,
                gamma = gamma
            )
        ),
        linetype = "dotted",
        colour = "gray",
        size = .3
    )  +
    
    geom_line(aes(x = W, y = A)) +
    
    annotate(
        "text",
        x = wmax - .4,
        y = T_fun(
            W = wmax - .4,
            c0 = c0,
            c1 = c1,
            gamma = gamma
        ) + .14,
        label = "$A=T$"
    ) +
    
    labs(y = "$A$",
         x = "$W$")  +
    
    theme_bw() +
    
    scale_x_continuous(breaks = Ws_a,
                       labels = as.character(c("$W_{0}$",
                                               "$W_{1}$",
                                               "$W_{2}$")),
                       expand = c(0 , 0)) +
    
    scale_y_continuous(
        breaks = c(0, tau / p),
        labels = c("0", "$A_{1}$"),
        expand = c(0, 0)
    ) +
    
    coord_cartesian(xlim = c(wmin, wmax),
                    ylim = c(-0.005, max(A_by_W_a_dt$A) + 1)) +
    theme(
        plot.title = element_text(size = rel(1.2), hjust = .5),
        legend.text = element_text(size = rel(1)),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.key.size = grid::unit(2, "lines"),
        legend.background = element_blank(),
        legend.position = c(0.7, 0.2),
        axis.line = element_line(colour = "black", size = .1),
        axis.ticks = element_line(colour = "black", size = .1),
        axis.text.y = element_text(colour = "black", margin = margin(
            t = 0,
            r = 5,
            b = 0,
            l = 0
        )),
        axis.text.x = element_text(colour = "black", margin = margin(
            t = 5,
            r = 0,
            b = 0,
            l = 0
        )),
        axis.title.x = element_text(
            vjust = 0,
            hjust = 1,
            margin = margin(
                t = 10,
                r = 0,
                b = 0,
                l = 0
            )
        ),
        axis.title.y = element_text(
            angle = 0,
            vjust = 1,
            hjust = 0,
            margin = margin(
                t = 0,
                r = 20,
                b = 0,
                l = 0
            )
        ),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank()
    )

dev.off()

#compiling the tikz into a pdf and displaying
if (!file.exists(paste0(tmptikz, "___LOCK"))) {
    setwd(here("fig/"))
    tools::texi2pdf(paste0(here("tikz/"), "/", plotname, ".tex"), clean = T)
    Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.21/bin/gswin64c.exe")
    embed_fonts(paste0(here("fig/"), "/", plotname, ".pdf"))
    system(paste(
        getOption("pdfviewer"),
        paste0(here("fig/"), "/", plotname, ".pdf")
    ))
    tmptikz_err <- setdiff(tmptikz_err, plotname)
} else{
    file.remove(paste0(tmptikz, "___LOCK"))
    tmptikz_err <- c(tmptikz_err, plotname)
}