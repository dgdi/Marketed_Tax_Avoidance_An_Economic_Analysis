#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@           Marketed Tax Avoidance:            @@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@             An Economic Analysis             @@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@                   Figure 2                   @@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# This script has been tested on Win 10, R-4.1.1, RStudio-1.4.1717
# the list of the required packages is at the start of the code.
#
# These files are distributed under BSD license. For more information, please check the license.txt file.
#
# Reference article:
# Li, J., Gamannossi degl'Innocenti, D. and Rablen, M.D., 2021. Marketed Tax Avoidance Schemes: An Economic Analysis.
#
# For any question, suggestion or comment,
# write to:
# mail@dgdi.me

#cleaning environment
rm(list = ls())

#loading libraries
library("here")
library("data.table")
library("stats")
library("VGAM")
library("ggplot2")
library("tikzDevice")
library("extrafont")

#sourcing general parameters
source(here("R/par.R"))

#sourcing functions
source(here("R/fun.R"))

#defining parameters
rho = .15

#progressive tax system parameters (similar to UK)
c0_prg <- 0
c1_prg <- .45
gamma_prg <- .00004

##defining the c0 for the linear tax system so that the tax burden is the same

#defining integrand function to compute difference in tax burden
avg_tax_num_integrand_fun <- function(W, c0, c1, gamma, mu, sigma) {
  T_fun(W, c0, c1, gamma) *
    dlnorm(x = W,
           meanlog = mu,
           sdlog = sigma)
  # dunif(x = W, min = 0, max = wbar)
  
}

#defining error function between tax burden of progressive and linear tax systems
se_diffburden_fun <- function(c0_lin,
                              c0_prg,
                              c1_prg,
                              gamma_prg,
                              mu,
                              sigma) {
  (
    integrate(
      avg_tax_num_integrand_fun,
      lower = 0,
      upper = wbar,
      c0 = c0_prg,
      c1 = c1_prg,
      gamma = gamma_prg,
      mu = mu,
      sigma = sigma
    )$value -
      integrate(
        avg_tax_num_integrand_fun,
        lower = 0,
        upper = wbar,
        c0 = c0_lin,
        c1 = 0,
        gamma = 5,
        mu = mu,
        sigma = sigma
      )$value
  ) ^ 2
}

#numerically identify c0_lin with same tax burden as prg tax system
c0_lin_opt <- optimize(
  f = se_diffburden_fun,
  interval = c(.00001, 1),
  c0_prg = c0_prg,
  c1_prg = c1_prg,
  gamma_prg = gamma_prg,
  mu = mu,
  sigma = sigma,
  tol = tol
)

#defining the parameters of the linear tax system with same tax burden of the progressive one
c0_lin <- c0_lin_opt$minimum
c1_lin <- 0
gamma_lin <-
  5 #this is unnecessary theoretically but needs to be specified to avoid numerical errors

#creating the list of tax systems to loop over
c0_c1_list <-
  list(
    "lin" = list(
      "c0" = c0_lin,
      "c1" = c1_lin,
      "gamma" = gamma_lin
    ),
    "prg" = list(
      "c0" = c0_prg,
      "c1" = c1_prg,
      "gamma" = gamma_prg
    )
  )

#identify numerically the W1 in linear/progressive tax system
W1_root_tot <- lapply(names(c0_c1_list), function(tax_system_temp) {
  #defining the tax system parameters for the loop
  c0_temp <- c0_c1_list[[tax_system_temp]]$c0
  c1_temp <- c0_c1_list[[tax_system_temp]]$c1
  gamma_temp <- c0_c1_list[[tax_system_temp]]$gamma
  
  #computing F_ for the different c0
  F__by_t_lin_root <- optimize(
    f = se_F__tot_fun,
    interval = c(tau, 2000000),
    c0 = c0_temp,
    c1 = c1_temp,
    gamma = gamma_temp,
    mu = mu,
    sigma = sigma,
    rho = rho,
    tau = tau,
    tol = tol
  )
  
  #computing the minimum W to have a non-negative post-audit disposable income
  min_W_temp <-
    Wmin_analyt_fun(
      c0 = c0_temp,
      c1 = c1_temp,
      gamma = gamma_temp,
      F_ = F__by_t_lin_root$minimum
    )
  
  #finding W1 as root
  W1_root_temp <- uniroot(
    f = W1_root_fun,
    interval = c(min_W_temp + 1, 3000000),
    c0 = c0_temp,
    c1 = c1_temp,
    gamma = gamma_temp,
    F_ = F__by_t_lin_root$minimum,
    rho = rho,
    check.conv = check_conv,
    tol = tol
  )
  
  #returning F_ and W1
  cbind(as.data.frame(F__by_t_lin_root),
        as.data.frame(W1_root_temp))
  
})

#naming the list W1
names(W1_root_tot) <- names(c0_c1_list)

#setting a convenient upper bound for W
wtop <- 3 * max(W1_root_tot[[1]]$root, W1_root_tot[[2]]$root)

#creating a dt with p, F_, Xi by W for: 1. fettered/unfettered taxpayers, 2. linear/progressive tax systems
vars_dt_tot_nostar <-
  rbindlist(lapply(names(c0_c1_list), function(tax_system_temp) {
    #defining the tax system parameters for the loop
    c0_temp <- c0_c1_list[[tax_system_temp]]$c0
    c1_temp <- c0_c1_list[[tax_system_temp]]$c1
    gamma_temp <- c0_c1_list[[tax_system_temp]]$gamma
    
    #retrieving the W1
    W1_root_temp <- W1_root_tot[[tax_system_temp]]
    
    #defining the W interval to be considered below W1
    Wi_below_temp <-
      c(
        seq(0, W1_root_temp$root - W1_root_temp$root / 5, length.out = 25),
        seq(
          W1_root_temp$root - W1_root_temp$root / 10,
          W1_root_temp$root,
          length.out = 75
        )
      )
    
    #define the data.table to store the relevant variables for W below W1
    vars_below_dt_temp <- data.table(W = Wi_below_temp)
    
    #computing F_
    vars_below_dt_temp[, F_ := W1_root_tot[[tax_system_temp]]$minimum]
    
    #computing pW
    vars_below_dt_temp[, pW := F_ / T_fun(W,
                                          c0 = c0_temp,
                                          c1 = c1_temp,
                                          gamma = gamma_temp)]
    
    #computing Xi
    vars_below_dt_temp[, Xi := Xi_fun(
      W = W,
      c0 = c0_temp,
      c1 = c1_temp,
      gamma = gamma_temp
    )]
    
    #defining the W interval to be considered above W1
    Wi_above_temp <- c(
      seq(
        W1_root_temp$root,
        W1_root_temp$root + W1_root_temp$root / 10,
        length.out = 50
      ),
      seq(W1_root_temp$root + W1_root_temp$root / 10,
          wtop, length.out = 50)
    )
    
    #define the data.table to store the relevant variables
    vars_above_dt_temp <- data.table(W = Wi_above_temp)
    
    #computing pW
    vars_above_dt_temp[, pW := pstar_funW(
      W = W,
      c0 = c0_temp,
      c1 = c1_temp,
      gamma = gamma_temp,
      rho = rho
    )]
    
    #computing F_
    vars_above_dt_temp[, F_ := (1 - (rho / (1 - pW))) * Xi_fun(
      W = W,
      c0 = c0_temp,
      c1 = c1_temp,
      gamma = gamma_temp
    )]
    
    #computing Xi
    vars_above_dt_temp[, Xi := Xi_fun(
      W = W,
      c0 = c0_temp,
      c1 = c1_temp,
      gamma = gamma_temp
    )]
    
    #binding the dts for above and below W1tilda1
    vars_tot_dt_temp <- rbindlist(list(vars_below_dt_temp,
                                       vars_above_dt_temp), use.names = TRUE)
    
    #set the tax system variable
    vars_tot_dt_temp[, tax_system := tax_system_temp]
    
  }))

#creating a dt with p, F_, Xi by W for: 1. fettered/unfettered taxpayers, 2. linear/progressive tax systems at W1
vars_dt_star <-
  rbindlist(lapply(names(c0_c1_list), function(tax_system_temp) {
    c0_temp <- c0_c1_list[[tax_system_temp]]$c0
    c1_temp <- c0_c1_list[[tax_system_temp]]$c1
    gamma_temp <- c0_c1_list[[tax_system_temp]]$gamma
    
    #retrieving the W1
    W1_root_temp <- W1_root_tot[[tax_system_temp]]
    
    #creating dt for W1
    data.table(
      W = W1_root_temp$root,
      F_ = W1_root_tot[[tax_system_temp]]$minimum,
      pW = W1_root_tot[[tax_system_temp]]$minimum / T_fun(
        W1_root_temp$root,
        c0 = c0_temp,
        c1 = c1_temp,
        gamma = gamma_temp
      ),
      Xi = Xi_fun(
        W = W1_root_temp$root,
        c0 = c0_temp,
        c1 = c1_temp,
        gamma = gamma_temp
      ),
      tax_system = tax_system_temp
    )
  }))

#binding the nostar and star dts
vars_dt_tot <- rbindlist(list(vars_dt_tot_nostar,
                              vars_dt_star))

#setting as upper bound for A the highest in the two tax systems
Atop <-
  max(
    T_fun(
      W = 187565.2,
      c0 = 0.2302917,
      c1 = 0,
      gamma = gamma_lin
    ),
    T_fun(
      W = 187565.2,
      c0 = 0,
      c1 = 0.45,
      gamma = 4e-05
    )
  )

#creating a dt with p by A for: 1. fettered/unfettered taxpayers, 2. linear/progressive tax systems
vars_dt_A <-
  rbindlist(lapply(names(c0_c1_list), function(tax_system_temp) {
    #defining the tax system parameters
    c0_temp <- c0_c1_list[[tax_system_temp]]$c0
    c1_temp <- c0_c1_list[[tax_system_temp]]$c1
    gamma_temp <- c0_c1_list[[tax_system_temp]]$gamma
    
    #define the data.table to store the relevant variables for W below W1
    vars_dt_temp <-
      data.table(Ai = seq(0, Atop, length.out = 100))
    
    #Compute pA the data.table to store the relevant variables for W below W1
    vars_dt_temp[, pA := p_funA(Ai,
                                c0 = c0_temp,
                                c1 = c1_temp,
                                gamma = gamma_temp,
                                rho = rho)]
    
    #adding manually the values for the limit A -> 0
    if (tax_system_temp == "lin") {
      vars_dt_temp[Ai == 0,
                   pA := vars_dt_temp[Ai > 0, pA][1]]
      
    } else{
      vars_dt_temp[Ai == 0, pA := (1 - rho)]
    }
    
    #set the tax system variable
    vars_dt_temp[, tax_system := tax_system_temp]
  }))


#### Plots ###############################################################
#define variable to store errors in tikz compilation
tmptikz_err <- NULL

#### Figure_2(a) - pW_by_W ###############################################################
plotname <- "Figure_2(a)"
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

ggplot(data = vars_dt_tot, aes(x = W, y = pW)) +
  
  geom_segment(
    aes(
      x = vars_dt_star[tax_system == "lin"]$W,
      y = .72,
      xend = vars_dt_star[tax_system == "lin"]$W,
      yend = vars_dt_star[tax_system == "lin"]$pW
    ),
    linetype = "dotted",
    colour = "gray",
    size = .3
  )  +
  
  geom_segment(
    aes(
      x = vars_dt_star[tax_system == "prg"]$W,
      y = .72,
      xend = vars_dt_star[tax_system == "prg"]$W,
      yend = vars_dt_star[tax_system == "prg"]$pW
    ),
    linetype = "dotted",
    colour = "gray",
    size = .3
  )  +
  
  geom_line(aes(linetype = tax_system)) +
  
  labs(y = "$p^{\\ast}(W)$",
       x = "$W$")  +
  
  theme_bw() +
  
  scale_x_continuous(
    breaks = c(0, W1_root_tot[["lin"]]$root, W1_root_tot[["prg"]]$root),
    labels = as.character(
      c(
        "$0$",
        "$\\hspace{-.1pt}W_{2}^{fl}$",
        "$\\hspace{.1pt}W_{2}^{pr}$"
      )
    ),
    expand = c(0 , 0)
  ) +
  
  scale_y_continuous(
    breaks = c(.72, .75, 1 - rho),
    labels = c("$0$",  "$.75$", "$1-\\phi$"),
    expand = c(0, 0)
  ) +
  
  coord_cartesian(xlim = c(0, wtop),
                  ylim = c(.72, 1 - rho)) +
  
  annotate(
    geom = "segment",
    x = 0,
    xend = 0,
    y = -Inf,
    yend = Inf,
    size = .1
  ) +
  annotate(
    geom = "segment",
    x = 0,
    xend = 0,
    y =  .72,
    yend = .75,
    linetype = "dotted",
    color = "white",
    size = .2
  ) +
  
  scale_linetype(name = "",
                 labels = c("Flat tax", "Progressive tax")) +
  
  theme(
    plot.title = element_text(size = rel(1.2), hjust = .5),
    legend.text = element_text(size = rel(1)),
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.key.size = grid::unit(2, "lines"),
    legend.background = element_blank(),
    legend.position = c(0.7, 0.2),
    axis.line.x = element_line(colour = "black", size = .1),
    axis.line.y = element_blank(),
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

#compiling the tikz into a pdf and displaying ***the compiler path needs to be set by the user***
if (!file.exists(paste0(tmptikz, "___LOCK"))) {
  #set the compiler path
  Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.21/bin/gswin64c.exe")
  setwd(here("fig/"))
  tools::texi2pdf(paste0(here("tikz/"), "/", plotname, ".tex"), clean = T)
  embed_fonts(paste0(here("fig/"), "/", plotname, ".pdf"))
  system(paste(getOption("pdfviewer"), paste0(here("fig/"), "/", plotname, ".pdf")))
  tmptikz_err <- setdiff(tmptikz_err, plotname)
} else{
  file.remove(paste0(tmptikz, "___LOCK"))
  tmptikz_err <- c(tmptikz_err, plotname)
}

#### Figure_2(b) ###############################################################
plotname <- "Figure_2(b)"
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

ggplot(data = vars_dt_A, aes(x = Ai, y = pA)) +
  
  geom_line(aes(linetype = tax_system)) +
  
  labs(y = "$p^{\\ast}(A)$",
       x = "$A$")  +
  
  theme_bw() +
  
  scale_x_continuous(breaks = c(0),
                     labels = c("$0$"),
                     expand = c(0 , 0)) +
  
  scale_y_continuous(
    breaks = c(.67, .7, 1 - rho),
    labels = c("$0$",  "$0.7$", "$1-\\phi$"),
    expand = c(0, 0)
  ) +
  
  coord_cartesian(xlim = c(0, 80000),
                  ylim = c(.67, 1 - rho)) +
  
  annotate(
    geom = "segment",
    x = 0,
    xend = 0,
    y = -Inf,
    yend = Inf,
    size = .1
  ) +
  annotate(
    geom = "segment",
    x = 0,
    xend = 0,
    y =  .67,
    yend = .7,
    linetype = "dotted",
    color = "white",
    size = .2
  ) +
  
  scale_linetype(name = "",
                 labels = c("Flat tax", "Progressive tax")) +
  
  theme(
    plot.title = element_text(size = rel(1.2), hjust = .5),
    legend.text = element_text(size = rel(1)),
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.key.size = grid::unit(2, "lines"),
    legend.background = element_blank(),
    legend.position = c(0.7, 0.2),
    axis.line.x = element_line(colour = "black", size = .1),
    axis.line.y = element_blank(),
    axis.ticks = element_line(colour = "black", size = .1),
    axis.text = element_text(colour = "black", margin = margin(
      t = 0,
      r = 5,
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

#compiling the tikz into a pdf and displaying ***the compiler path needs to be set by the user***
if (!file.exists(paste0(tmptikz, "___LOCK"))) {
  #set the compiler path
  Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.21/bin/gswin64c.exe")
  setwd(here("fig/"))
  tools::texi2pdf(paste0(here("tikz/"), "/", plotname, ".tex"), clean = T)
  embed_fonts(paste0(here("fig/"), "/", plotname, ".pdf"))
  system(paste(getOption("pdfviewer"), paste0(here("fig/"), "/", plotname, ".pdf")))
  tmptikz_err <- setdiff(tmptikz_err, plotname)
} else{
  file.remove(paste0(tmptikz, "___LOCK"))
  tmptikz_err <- c(tmptikz_err, plotname)
}
