#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@           Marketed Tax Avoidance:            @@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@             An Economic Analysis             @@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@                   Figure_3                   @@@@@@@@@@@@@@@@@@@@@@@@@
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

#defining the average tax rates to be considered in the plots
c0_loop <- seq(from = .07 ,
               to = .9,
               length.out = 250)

# Linear tax system ###############################################################

#defining the c1 and gamma for the linear tax system (we will loop on c0 to adjust average tax rate)
c1_lin <- 0
gamma_lin <- 5 #this is unnecessary theoretically but needs to be specified to avoid numerical errors

#computing F_ for the different c0
F__by_t_lin_root <- sapply(c0_loop, function(c0_temp) {
  optimize(
    f = se_F__tot_fun,
    interval = c(tau, 2000000),
    c0 = c0_temp,
    c1 = c1_lin,
    gamma = gamma_lin,
    mu = mu,
    sigma = sigma,
    rho = rho,
    tau = tau,
    tol = tol
  )
})


#computing W0 for the different c0 using F_
W0_by_t_lin_root <- sapply(seq_along(c0_loop), function(id_temp) {
  min_W_temp <-
    Wmin_analyt_fun(
      c0 = c0_loop[id_temp],
      c1 = c1_lin,
      gamma = gamma_lin,
      F_ = F__by_t_lin_root[1, id_temp]$minimum
    )
  
  optimize(
    f = W0_se_fun,
    interval = c(min_W_temp + .1, 2000000),
    c0 = c0_loop[id_temp],
    c1 = c1_lin,
    gamma = gamma_lin,
    rho = rho,
    F_ = F__by_t_lin_root[1, id_temp]$minimum,
    tol = tol
  )
})

#creating dt for linear tax system variables
var_by_t_lin_dt <- data.table(
  "c0" = c0_loop,
  "c1" = 0,
  "W0" = unlist(W0_by_t_lin_root[1,]),
  "F_" = unlist(F__by_t_lin_root[1,]),
  "tax_system" = "lin"
  
)

#computing total avoidance in the linear tax setting
var_by_t_lin_dt[, A := totpop * A_lin_fun(W0 = W0, c0 = c0)]
#computing total revenues from avoidance enforcement in the linear tax setting
var_by_t_lin_dt[, RA := totpop * RA_lin_fun(
  W0 = W0,
  c0 = c0,
  rho = rho
)]
#computing total revenues from compliance in the linear tax setting
var_by_t_lin_dt[, RC := totpop * RC_lin_fun(W0 = W0, c0 = c0)]
#computing total revenues in the linear tax setting
var_by_t_lin_dt[, R := totpop * R_lin_fun(
  W0 = W0,
  c0 = c0,
  rho = rho
)]



# Progressive tax system ###############################################################

#defining the c0 and gamma for the progressive tax system (we will loop on c1 to adjust average tax rate)
c0_prg <- 0
gamma_prg <- .00004

#defining integrand function to compute cumulative tax liabilities
avg_tax_num_integrand_fun <- function(W, c0, c1, gamma, mu, sigma) {
  T_fun(W, c0, c1, gamma) *
    dlnorm(x = W,
           meanlog = mu,
           sdlog = sigma)
}

#defining integrand function to compute cumulative income
avg_tax_den_integrand_fun <- function(W, mu, sigma) {
  W * dlnorm(x = W,
             meanlog = mu,
             sdlog = sigma)
}

#defining error function to identify c1 with a given average tax rate
root_avg_t_integrand_fun <-
  function(c1, value, c0, gamma, mu, sigma) {
    (
      integrate(
        avg_tax_num_integrand_fun,
        lower = 0,
        upper = wbar,
        c0 = c0,
        c1 = c1,
        gamma = gamma,
        mu = mu,
        sigma = sigma
      )$value /
        integrate(
          avg_tax_den_integrand_fun,
          lower = 0,
          upper = wbar,
          mu = mu,
          sigma = sigma
        )$value
    ) -
      value
    
  }

#computing c1 that entails a given average t_fun
c1_by_t <- sapply(c0_loop, function(c0_temp) {
  uniroot(
    f = root_avg_t_integrand_fun,
    interval = c(.01, 10),
    value = c0_temp,
    c0 = c0_prg,
    gamma = gamma_prg,
    mu = mu,
    sigma = sigma,
    check.conv = check_conv,
    tol = tol
  )
})

#identifying the progressive tax systems with average tax rate lower than 1
c1_ok_id <- which(unlist(c1_by_t[1,]) < 1)
#selecting the c1_by_t with marginal tax rate <1
c1_by_t_ok <- unlist(c1_by_t[1,])[c1_ok_id]
#selecting the c0_loop that have a marginal tax rate <1
c0_loop_ok <- c0_loop[c1_ok_id]

#computing F_ for the different c0
F__by_t_prg_root <- sapply(c1_by_t_ok, function(c1_temp) {
  optimize(
    f = se_F__tot_fun,
    interval = c(tau, 200000),
    c0 = c0_prg,
    c1 = c1_temp,
    gamma = gamma_prg,
    mu = mu,
    sigma = sigma,
    rho = rho,
    tau = tau,
    tol = tol
  )
})

#computing W0 for the different c0 using F_
W0_by_t_prg_root <- sapply(seq_along(c1_by_t_ok), function(id_temp) {
  #computing the minimum W to have a non-negative post-audit disposable income
  min_W_temp <-
    Wmin_analyt_fun(
      c0 = c0_prg,
      c1 = c1_by_t_ok[id_temp],
      gamma = gamma_prg,
      F_ = F__by_t_prg_root[1, id_temp]$minimum
    )
  
  optimize(
    f = W0_se_fun,
    interval = c(min_W_temp + .1, 200000),
    c0 = c0_prg,
    c1 = c1_by_t_ok[id_temp],
    gamma = gamma_prg,
    rho = rho,
    F_ = F__by_t_prg_root[1, id_temp]$minimum,
    tol = tol
  )
})

#creating dt for progressive tax system variables
var_by_t_prg_dt <- data.table(
  "c0" = c0_loop_ok,
  "c1" = c1_by_t_ok,
  "W0" = unlist(W0_by_t_prg_root[1,]),
  "F_" = unlist(F__by_t_prg_root[1,]),
  "tax_system" = "prg"
)



# Compute Avoidance, Revenues (from  both Enforcement and Compliance) for the linear and progressive tax systems ###############################################################


#defining a vectorized function (that can take vectors as input for the variables "Wdwn", "Wup", "c1")
#to compute the cumulative tax liabilities (T(W)g(W)) between two levels of income ("Wdwn", "Wup")
tax_liabilities_Vfun <- Vectorize(function(Wdwn, Wup, c0, c1, gamma, mu, sigma) {
  integrate(
    avg_tax_num_integrand_fun,
    lower = Wdwn,
    upper = Wup,
    c0 = c0,
    c1 = c1,
    gamma = gamma,
    mu = mu,
    sigma = sigma
  )$value
}, vectorize.args = c("Wdwn", "Wup", "c1"))

#computing total avoidance in the progressive tax system
var_by_t_prg_dt[, A := totpop * tax_liabilities_Vfun(
  Wdwn = W0,
  Wup = wbar,
  c0 = c0_prg,
  c1 = c1,
  gamma = gamma_prg,
  mu = mu,
  sigma = sigma
)]

#computing total revenues from avoidance enforcement in the progressive tax system
var_by_t_prg_dt[, RA := A * rho]

#computing total revenues from compliance in the progressive tax system
var_by_t_prg_dt[, RC := totpop * tax_liabilities_Vfun(
  Wdwn = 0,
  Wup = W0,
  c0 = c0_prg,
  c1 = c1,
  gamma = gamma_prg,
  mu = mu,
  sigma = sigma
)]

#computing total revenues in the progressive tax system
var_by_t_prg_dt[, R := RA + RC]

#creating dt for revenues in a perfect compliance setting
var_by_t_pc_dt <- data.table(
  "c0" = c0_loop_ok,
  "c1" = c1_by_t_ok,
  "W0" = rep(wbar, length(c0_loop_ok)),
  "F_" = 0,
  "tax_system" = "pc"
)

#compute revenues collected in the linear tax system with perfect compliance
var_by_t_pc_dt[, R_lin_pc := totpop * R_lin_fun(wbar,
                                                c0 = c0,
                                                rho = rho)]

#compute revenues collected in the progressive tax system with perfect compliance
var_by_t_pc_dt[, R_prg_pc := totpop * tax_liabilities_Vfun(
  Wdwn = 0,
  Wup = wbar,
  c0 = c0_prg,
  c1 = c1,
  gamma = gamma_prg,
  mu = mu,
  sigma = sigma
)]

#binding the two dts of linear and progressive tax system
var_by_t_tot_dt_whole <- rbindlist(list(var_by_t_lin_dt,
                                        var_by_t_prg_dt))

#computing the cumulative distribution function of W0
var_by_t_tot_dt_whole[, GW0 := plnorm(q = W0, meanlog = mu, sdlog = sigma)]



#### Plots ###############################################################
#considering only the average tax rates available in both tax systems
lower_max_c0 <-
  min(var_by_t_tot_dt_whole[, max(c0), by = tax_system]$V1)
higher_min_c0 <-
  max(var_by_t_tot_dt_whole[, min(c0), by = tax_system]$V1)
var_by_t_tot_dt <-
  var_by_t_tot_dt_whole[!(c0 > lower_max_c0 | c0 < higher_min_c0)]

#expressing total avoidance in billions
var_by_t_tot_dt[, A := A / 1000000000]
#expressing total revenues from avoidance enforcement in billions
var_by_t_tot_dt[, RA := RA / 1000000000]
#expressing total revenues from compliance in billions
var_by_t_tot_dt[, RC := RC / 1000000000]
#expressing total revenues in billions
var_by_t_tot_dt[, R := R / 1000000000]
#expressing total revenues of perfect compliance in billions
var_by_t_pc_dt[, R_lin_pc := R_lin_pc / 1000000000]
var_by_t_pc_dt[, R_prg_pc := R_prg_pc / 1000000000]

#define variable to store errors in tikz compilation
tmptikz_err <- NULL

#### Figure_3(a) ###############################################################
plotname <- "Figure_3(a)"
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

ggplot(data = var_by_t_tot_dt, aes(x = c0, y = R)) +
  
  geom_line(aes(linetype = tax_system)) +
  
  geom_line(
    data = var_by_t_pc_dt,
    aes(x = c0, y = R_lin_pc),
    linetype = "dotted",
    size = .3,
    color = "gray"
  ) +
  
  annotate(
    "text",
    x = var_by_t_pc_dt[R_lin_pc > (ceiling(max(var_by_t_tot_dt$R) /
                                             100) * 100) * .95]$c0[1] * 1.35,
    y = (ceiling(max(var_by_t_tot_dt$R) / 100) * 100) * .96,
    label = "Full compliance"
  ) +
  
  labs(y = "$\\mathbf{E}(R)$",
       x = "$\\mathbf{T}/\\mathbf{W}$")  +
  
  theme_bw() +
  
  scale_x_continuous(
    breaks = seq(.1, .5, by = .1),
    labels = as.character(seq(.1, .5, by = .1)),
    expand = c(0 , 0)
  ) +
  
  scale_y_continuous(breaks = 0,
                     labels = "0",
                     expand = c(0, 0)) +
  
  coord_cartesian(xlim = c(0.05, .53),
                  ylim = c(0, (ceiling(
                    max(var_by_t_tot_dt$R) / 100
                  ) * 100) * .98)) +
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

#### Figure_3(b) ###############################################################
plotname <- "Figure_3(b)"
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

ggplot(data = var_by_t_tot_dt, aes(x = c0, y = GW0)) +
  
  geom_line(aes(linetype = tax_system)) +
  
  labs(y = "$G(W_{0})$",
       x = "$\\mathbf{T}/\\mathbf{W}$")  +
  
  theme_bw() +
  
  scale_x_continuous(
    breaks = seq(.1, .5, by = .1),
    labels = as.character(seq(.1, .5, by = .1)),
    expand = c(0 , 0)
  ) +
  
  scale_y_continuous(
    breaks = seq(0, 1, by = .25),
    labels = as.character(seq(0, 1, by = .25)),
    expand = c(0, 0)
  ) +
  
  coord_cartesian(xlim = c(0.05, .53),
                  ylim = c(0, 1)) +
  
  scale_linetype(name = "",
                 labels = c("Flat tax", "Progressive tax")) +
  
  theme(
    plot.title = element_text(size = rel(1.2), hjust = .5),
    legend.text = element_text(size = rel(1)),
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.key.size = grid::unit(2, "lines"),
    legend.background = element_blank(),
    legend.position = c(0.8, 0.6),
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