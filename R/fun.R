#####for W0_A_R_by_t_prog.R

#marginal tax function (see also the analytic.nb file)
#Tprime_fun <- function(Wi, c0, c1, gamma){(c0+c1*(1 - exp(-gamma*Wi)))}

#defining average tax 
t_fun <- function(Wi, c0, c1, gamma){(c0 + c1*(1 - ((1-exp(-gamma*Wi))/(gamma*Wi))))}

T_fun <- function(Wi, c0, c1, gamma){Wi*t_fun(Wi, c0, c1, gamma)}

#defining inverse of average tax 
Tinv_fun <- function(Tax, c0, c1, gamma){
  (1/((c0 + c1)*gamma)) * (c1 + Tax*gamma + 
                             c0*lambertW(-((c1*exp(((-Tax - c1/gamma)*gamma)/(c0 + c1)))/(
                               c0 + c1))) + 
                             c1*lambertW(-((c1*exp(((-Tax - c1/gamma)*gamma)/(c0 + c1)))/(
                               c0 + c1)))
  )
  }
# T_fun(Wi = 1000, c0 = .25, c1 = 0)

# plot(
#   seq(from = 0, to = 30000, by = 100),
#   T_fun(Wi = seq(from = 0, to = 30000, by = 100), c0 = 0, c1 = .45, gamma = .00004)
# )
# plot(
#   seq(from = 0, to = 30000, by = 100),
#   t_fun(Wi = seq(from = 0, to = 30000, by = 100), c0 = 0, c1 = .45, gamma = .00004)
# )


#defining truthful after-tax income (Xi) 
Xi_fun <- function(Wi, c0, c1, gamma){Wi - T_fun(Wi, c0, c1, gamma)}
# Xi_fun(Wi = 1000, c0 = .25, c1 = 0, gamma = .00003)

# #analytic/numeric (derivation in Mathematica file analytic.nb)
# #to be used when t_fun=c0+c1*(1 - ((1-exp(-gamma*Wi))/(gamma*Wi))))
Wmin_analyt_fun <- function(c0, c1, gamma, F_){
  (c1 - F_*gamma + (-1 + c0 + c1)*lambertW(-((
    c1*exp((-c1 + F_*gamma)/(-1 + c0 + c1)))/(-1 + c0 + c1))))/((-1 + c0 + c1)*gamma)
}
#defining EU for truthful reporters 
EU_A_0_fun <- function(W0, c0, c1, gamma){log(Xi_fun(W0, c0, c1, gamma))}
# EU_A_0_fun(Wi = 1000, c0 = .25, c1 = 0)-log(750)

###defining EU for wealth constrained taxpayer
EU_A_WC_fun <- function(W0, c0, c1, gamma, rhol, rhov, F_){
  ((rhol*rhov)*
         log(Xi_fun(W0, c0, c1, gamma) - F_ )) + ############## 
  ((1-(rhol*rhov))*
     (log(Xi_fun(W0, c0, c1, gamma) - F_ +####################
          T_fun(W0, c0, c1, gamma))))}
# EU_A_WC_fun(Wi = 10000, c0 = .25, c1 = 0, tau = tau)

#defining D''_0
W0_root_fun <- function(W0, W2, c0, c1, gamma, rhol, rhov, F_ ){
  (EU_A_0_fun(W0, c0, c1, gamma) - EU_A_WC_fun(W0, c0, c1, gamma, rhol, rhov, F_ ))
}


#defining function to find W0 as a numerical root
W0_rootfinding_fun <- function(c0, c1, gamma, rhol, rhov, F_){
  #computing the minimum W to have a non-negative W^u
  min_W_temp <-  Wmin_analyt_fun(c0 = c0, c1 = c1, gamma = gamma, F_ = F_)
  #finding W1 as root 
  uniroot(f = W0_root_fun, 
          interval = c(min_W_temp+.1, 2000000), 
          c0 = c0, 
          c1 = c1,
          gamma = gamma,
          rhol = rhol, rhov = rhov,
          F_ = F_,
          check.conv = check_conv,
          tol = tol)
}

#defining quadratic loss function to perform optimization (instead of root searching)
W0_se_fun <- function(Wi, c0, c1, gamma, rhol, rhov, F_){
  (EU_A_0_fun(Wi, c0, c1, gamma) - EU_A_WC_fun(Wi, c0, c1, gamma, rhol, rhov, F_))^2
}

#defining function to find W0 as a numerical min of se
W0_minsefinding_fun <- function(c0, c1, gamma, rhol, rhov, F_){
  #computing the minimum W to have a non-negative W^u
  min_W_temp <-  Wmin_analyt_fun(c0 = c0, c1 = c1, gamma = gamma, F_ = F_)
  #finding W1 as root 
  optimize(f = W0_se_fun, 
           interval = c(min_W_temp+.1, 2000000), 
           c0 = c0, 
           c1 = c1,
           gamma = gamma,
           rhol = rhol, rhov = rhov,
           F_ = F_,
           tol = tol)
}

#defining optimal avoidance from the demand side only (out of equilibrium)
Ai_outeq_fun <- function(Wi, c0, c1, gamma, rhol, rhov){
  (1/pstar_funW(Wi, c0, c1, gamma, rhol, rhov))*(1-((rhol*rhov)/(1/pstar_funW(Wi, c0, c1, gamma, rhol, rhov))))*Xi_fun(Wi, c0, c1, gamma)
}


#defining optimal avoidance (in equilibrium)
Ai_eq_fun <- function(Wi, c0, c1, gamma, rhol, rhov){
  ((1 - (1-((Wi-2*Xi_fun(Wi, c0, c1, gamma) +
               sqrt(((Wi-2*Xi_fun(Wi, c0, c1, gamma))^2) + (4*(rhol*rhov)*(Wi-Xi_fun(Wi, c0, c1, gamma))*Xi_fun(Wi, c0, c1, gamma))))/
              (2*(Wi-Xi_fun(Wi, c0, c1, gamma)))
  )) - (rhol*rhov))/((1-((Wi-2*Xi_fun(Wi, c0, c1, gamma) +
                            sqrt(((Wi-2*Xi_fun(Wi, c0, c1, gamma))^2) + (4*(rhol*rhov)*(Wi-Xi_fun(Wi, c0, c1, gamma))*Xi_fun(Wi, c0, c1, gamma))))/
                           (2*(Wi-Xi_fun(Wi, c0, c1, gamma)))
  ))*(1 - (1-((Wi-2*Xi_fun(Wi, c0, c1, gamma) +
                 sqrt(((Wi-2*Xi_fun(Wi, c0, c1, gamma))^2) + (4*(rhol*rhov)*(Wi-Xi_fun(Wi, c0, c1, gamma))*Xi_fun(Wi, c0, c1, gamma))))/
                (2*(Wi-Xi_fun(Wi, c0, c1, gamma)))
  )))))*Xi_fun(Wi, c0, c1, gamma)}
# Ai_fun(Wi = 1000, c0 = .25, c1 = 0, gamma = .00003, p = .19, rhol = rhol, rhov = rhov)



#defining the error function to identify W2
W2_root_fun <- function(Wi, c0, c1, gamma, rhol, rhov, F_){
  (
    Ai_eq_fun(Wi, c0, c1, gamma, rhol, rhov) *
      pstar_funW(Wi, c0, c1, gamma, rhol, rhov)
  ) - F_ #############
}

#defining function to find W2 as a numerical root
W2_rootfinding_fun <- function(c0, c1, gamma, rhol, rhov, F_){
  #computing the minimum W to have a non-negative W^u
  min_W_temp <-  Wmin_analyt_fun(c0 = c0, c1 = c1, gamma = gamma, F_ = F_)
  #finding W1 as root 
  uniroot(f = W2_root_fun, 
          interval = c(min_W_temp+.1, 2000000), 
          c0 = c0, 
          c1 = c1,
          gamma = gamma,
          rhol = rhol, rhov = rhov,
          F_ = F_,
          check.conv = check_conv,
          tol = tol)
}
# W2_rootfinding_fun(c0 = .28021974531250876,
#                          c1 = 0,
#                          gamma = 10,
#                          rhol = .35,
#                          rhov = .85,
#                         F_ = 4800
#                         )$root

#defining the error function to identify W2
W2_se_fun <- function(Wi, c0, c1, gamma, rhol, rhov, F_){
  ((
    Ai_eq_fun(Wi, c0, c1, gamma, rhol, rhov) *
      pstar_funW(Wi, c0, c1, gamma, rhol, rhov)
  ) - F_ )^2 #######################################
}

#defining function to find W2 as a numerical min of se
W2_minsefinding_fun <- function(c0, c1, gamma, rhol, rhov, F_){
  #computing the minimum W to have a non-negative W^u
  min_W_temp <-  Wmin_analyt_fun(c0 = c0, c1 = c1, gamma = gamma, F_ = F_)
  #finding W1 as root 
  optimize(f = W2_se_fun, 
           interval = c(min_W_temp+.1, 2000000), 
           c0 = c0, 
           c1 = c1,
           gamma = gamma,
           rhol = rhol, rhov = rhov,
           F_ = F_,
           tol = tol)
}

# #definition of dWdF
# dD0dW <- function(W0, F_, c0, c1, gamma, rhol, rhov){
#   (F_*(1-pstar_funW(W0, c0, c1, gamma, rhol, rhov))*(
#     (F_/pstar_funW(W0, c0, c1, gamma, rhol, rhov))-Ai_eq_fun(W0, c0, c1, gamma, rhol, rhov))*(
#       1-t_fun(W0, c0, c1, gamma)))/
#   (Xi_fun(W0, c0, c1, gamma)*(
#     Xi_fun(W0, c0, c1, gamma)-F_)*(
#       Xi_fun(W0, c0, c1, gamma)+(F_/pstar_funW(W0, c0, c1, gamma, rhol, rhov))*(1-pstar_funW(W0, c0, c1, gamma, rhol, rhov)))
#    )
# }
# dD0dF <- function(W0, F_, c0, c1, gamma, rhol, rhov){
#   - ((1-pstar_funW(W0, c0, c1, gamma, rhol, rhov))*((F_/pstar_funW(W0, c0, c1, gamma, rhol, rhov))-Ai_eq_fun(W0, c0, c1, gamma, rhol, rhov)))/((Xi_fun(W0, c0, c1, gamma)-F_)*(Xi_fun(W0, c0, c1, gamma)+((F_/pstar_funW(W0, c0, c1, gamma, rhol, rhov))*(1-pstar_funW(W0, c0, c1, gamma, rhol, rhov)))))
# }
# dWdF  <- function(W0, F_, c0, c1, gamma, rhol, rhov){
#   -dD0dF(W0, F_, c0, c1, gamma, rhol, rhov)/dD0dW(W0, F_, c0, c1, gamma, rhol, rhov)
# }
#second derivative of D0'' wrt W (from mathematica file analytic)
dDppdW_fun  <- function(Wi, F_, c0, c1, gamma, rhol, rhov){
  if(c1 > 0){
    (-1 + (rhol*rhov))/(F_ - Wi) + 
      ((c1 - (-1 + c0 + c1)*exp(Wi*gamma))*gamma)/
      (c1 + exp(Wi*gamma)*(-c1 + (-1 + c0 + c1)*
                               Wi*gamma)) + 
      ((-c1 + (-1 + c0 + c1)*exp(Wi*gamma))*(rhol*rhov)*
         gamma)/(c1 + exp(Wi*gamma)*
                      (-c1 + (F_ + (-1 + c0 + c1)*Wi)*gamma))
  } else if(c1 == 0){
    (-1 + (rhol*rhov))/(F_ - Wi) - 1/Wi + ((-1 + c0)*(rhol*rhov))/(F_ + (-1 + c0)*Wi)  
    }
}
##check
#second derivative of D0'' wrt F_ (A10 in the paper)
dDppdF__fun  <- function(Wi, F_, c0, c1, gamma, rhol, rhov){
  -(
    ((rhov*rhol)/(Xi_fun(Wi, c0, c1, gamma)-F_)) + ((1-(rhov*rhol))/(Wi-F_))
  )
}
#second derivative of W wrt F_ using IFT
dWdF <- function(Wi, F_, c0, c1, gamma, rhol, rhov){
  -dDppdF__fun(Wi, F_, c0, c1, gamma, rhol, rhov)/dDppdW_fun(Wi, F_, c0, c1, gamma, rhol, rhov)
}

# dDppdW_fun(Wi = 6000, F_ = 5500, c0 = .01, c1 = .45, gamma = .00004, rhol = 1, rhov = .25)
# dDppdF__fun(Wi = 6000, F_ = 5500, c0 = .01, c1 = .45, gamma = .00004, rhol = 1, rhov = .25)
# dWdF(Wi = 6000, F_ = 5500, c0 = .01, c1 = .45, gamma = .00004, rhol = 1, rhov = .25)
# 
# dDppdW_fun(Wi = 24969.99, F_ = 4600, c0 = .28022, c1 = 0, gamma = 5, rhol = 1, rhov = .2975)
# dDppdF__fun(Wi = 24969.99, F_ = 4600, c0 = .28022, c1 = 0, gamma = 5, rhol = 1, rhov = .2975)
# dWdF(Wi = 24969.99, F_ = 4600, c0 = .28022, c1 = 0, gamma = 5, rhol = 1, rhov = .2975)


#error  function F_
e_F__fun <- function(F_, W0, W2, tau, mu, sigma, c0, c1, gamma, rhol, rhov){
  ((tau + (
    (plnorm(q = W2, meanlog = mu, sdlog = sigma) - plnorm(q = W0, meanlog = mu, sdlog = sigma))/
      (dWdF(W0, F_, c0, c1, gamma, rhol, rhov)*dlnorm(x = W0, meanlog = mu, sdlog = sigma))
  )
  ) - F_)
}

#squared error  function F_
se_F__fun <- function(F_, W0, W2, tau, mu, sigma, c0, c1, gamma, rhol, rhov){
    ((tau + (
    (plnorm(q = W2, meanlog = mu, sdlog = sigma) - plnorm(q = W0, meanlog = mu, sdlog = sigma))/
      (dWdF(W0, F_, c0, c1, gamma, rhol, rhov)*dlnorm(x = W0, meanlog = mu, sdlog = sigma))
  )
  ) - F_)^2
}

#error function F_ computing W0 and W2
e_F__tot_fun <- function(F_, c0, c1, gamma, rhol, rhov, mu, sigma, tau){
  #computing the minimum W to have a non-negative W^u
  min_W_temp <-  Wmin_analyt_fun(c0 = c0, c1 = c1, gamma = gamma, F_ = F_)    
  #finding W1 as root 
  W0_by_rho_lin_root <- optimize(f = W0_se_fun, 
                                 interval = c(min_W_temp, 2000000), 
                                 c0 = c0, 
                                 c1 = c1,
                                 gamma = gamma,
                                 rhol = rhol, rhov = rhov,
                                 F_ = F_,
                                 tol = tol)$minimum
  W2_by_rho_lin_root <- optimize(f = W2_se_fun, 
                                 interval = c(min_W_temp, 2000000), 
                                 c0 = c0, 
                                 c1 = c1,
                                 gamma = gamma,
                                 rhol = rhol, rhov = rhov,
                                 F_ = F_,
                                 tol = tol)$minimum
  e_F__fun(F_,
            W0 = W0_by_rho_lin_root,
            W2 = W2_by_rho_lin_root,
            c0 = c0, c1 = c1, gamma = gamma, 
            rhol = rhol, rhov = rhov,
            mu = mu, sigma = sigma,
            tau = tau)
}


#squared error function F_ computing W0 and W2
se_F__tot_fun <- function(F_, c0, c1, gamma, rhol, rhov, mu, sigma, tau){
  #computing the minimum W to have a non-negative W^u
  min_W_temp <-  Wmin_analyt_fun(c0 = c0, c1 = c1, gamma = gamma, F_ = F_)    
  #finding W1 as root 
  W0_by_rho_lin_root <- optimize(f = W0_se_fun, 
                                 interval = c(min_W_temp, 2000000), 
                                 c0 = c0, 
                                 c1 = c1,
                                 gamma = gamma,
                                 rhol = rhol, rhov = rhov,
                                 F_ = F_,
                                 tol = tol)$minimum
  W2_by_rho_lin_root <- optimize(f = W2_se_fun, 
                                 interval = c(min_W_temp, 2000000), 
                                 c0 = c0, 
                                 c1 = c1,
                                 gamma = gamma,
                                 rhol = rhol, rhov = rhov,
                                 F_ = F_,
                                 tol = tol)$minimum
  se_F__fun(F_,
            W0 = W0_by_rho_lin_root,
            W2 = W2_by_rho_lin_root,
            c0 = c0, c1 = c1, gamma = gamma, 
            rhol = rhol, rhov = rhov,
            mu = mu, sigma = sigma,
            tau = tau)
}


#function to compute total avoidance in the economic system for a linear tax system
#(derivation in Mathematica file analytic.nb)
A_lin_fun <- function(W0, c0){
  17105.3769312732*c0*(1 + erf(12.33572806744155 - (1.1609233137739046*log(W0))))
}
#function to compute revenues collected from avoidance enforcement for a linear tax system
#(derivation in Mathematica file analytic.nb)
RA_lin_fun <- function(W0, c0, rhol, rhov){
  17105.3769312732*c0*(rhol*rhov)*(1 + erf(12.33572806744155 - (1.1609233137739046*log(W0))))
}
#function to compute revenues collected from compliance (honest reporting) for a linear tax system
#(derivation in Mathematica file analytic.nb)
RC_lin_fun <- function(W0, c0){
  17105.3769312732*c0*erfc(
    12.33572806744155 - 1.1609233137739046*log(W0))
}
#function to compute total revenues collected for a linear tax system
#(derivation in Mathematica file analytic.nb)
R_lin_fun <- function(W0, c0, rhol, rhov){
  RC_lin_fun(W0, c0) + RA_lin_fun(W0, c0, rhol, rhov)
}

#####for p_by_W_by_F
#defining the error function to identify W1
W1_root_fun <- function(Wi, c0, c1, gamma, F_, rhol, rhov){
  T_fun(Wi, c0, c1, gamma) - (F_/pstar_funW(Wi, c0, c1, gamma, rhol, rhov)) 
}

#####for A_not_eq
# #defining optimal avoidance (outside equilibrium) 
Ai_fun <- function(Wi, c0, c1, gamma, p, rhol, rhov){((1 - p - (rhol*rhov))/(p*(1 - p)))*Xi_fun(Wi, c0, c1, gamma)}




sciW_tikz_fun <- function(tex_file, tex_path, tikz_file, parent_file, caption = "caption", label = "", scale = 1){
  
  tikz_figure_lines <- readLines(tikz_file)
  tikz_figure_lines_start <- which(tikz_figure_lines=="\\begin{tikzpicture}[x=1pt,y=1pt]")
  tikz_figure_lines_end <- which(tikz_figure_lines=="\\end{tikzpicture}")
  tikz_figure <- paste(tikz_figure_lines[tikz_figure_lines_start:tikz_figure_lines_end],
                       collapse = "\n"
  )
  tmp <- paste0("%TCIDATA{LaTeXparent=0,0,", parent_file, "}
              \\begin{figure}[h]
              \\begin{center}
              \\scalebox{", scale,"}{",
                tikz_figure,
                "}
\\captionsetup{labelformat=empty}
\\caption{", caption, "}
\\label{", label, "}
\\end{center}
\\end{figure}")
  
  write(tmp, file = paste0(tex_path, "/", tex_file))
  
  
  
}


# #checks
# R_lin_fun_check <- function(W0, c0, rhol, rhov){
#   17105.3769312732*c0*(rhol*rhov)*(1 +
#     erf(12.33572806744155 - 1.1609233137739046*log(W0))) +
#   17105.3769312732*c0*erfc(
#     12.33572806744155 - 1.1609233137739046*log(W0))
# }
# A_lin_fun(W0 = 20000, c0 = .2)
# RA_lin_fun(W0 = 20000, c0 = .2, rhol = rhol, rhov = rhov)
# RC_lin_fun(W0 = 20000, c0 = .2)
# R_lin_fun(W0 = 20000, c0 = .2, rhol = rhol, rhov = rhov)
# R_lin_fun_check(W0 = 20000, c0 = .2, rhol = rhol, rhov = rhov)
# #there are small differences in A2/A_lin wrt M, but the plots are the same
# plot(100:100000, A_lin_fun(100:100000, c0 = .2, rhol = rhol, rhov = rhov))


# #marginal tax function (see alse the analytic.nb file)
# #Tprime_fun <- function(Wi, c0, c1, gamma){(c0+c1*(1 - exp(-gamma*Wi)))}
# 
# #defining average tax 
# t_fun <- function(Wi, c0, c1, gamma){(c0+c1*(1 - ((1-exp(-gamma*Wi))/(gamma*Wi))))}
# 
# #defining tax liabilities 
# T_fun <- function(Wi, c0, c1, gamma){Wi*t_fun(Wi, c0, c1, gamma)}
# # T_fun(Wi = 1000, c0 = .25, c1 = 0)
# 
# #defining truthful after-tax income (Xi) 
# Xi_fun <- function(Wi, c0, c1, gamma){Wi - T_fun(Wi, c0, c1, gamma)}
# # Xi_fun(Wi = 1000, c0 = .25, c1 = 0, gamma = .00003)
# 
# #@@@@@@@@@@@ functions to compute minimum W s.t. W^u is positive
# #analytic/numeric (derivation in Mathematica file analytic.nb)
# #to be used when t_fun=c0+c1*(1 - ((1-exp(-gamma*Wi))/(gamma*Wi))))
# Wmin_analyt_fun <- function(c0, c1, gamma, tau){
#   (c1 - tau*gamma + (-1 + c0 + c1)*lambertW(-((
#   c1*exp((-c1 + tau*gamma)/(-1 + c0 + c1)))/(-1 + c0 + c1))))/((-1 + c0 + c1)*gamma)
# }
# 
# #numeric (more generic, allows for arbitrary tax function, but likely less precise)
# #set error  function to find minimum W s.t. W^u is positive
# root_Wmin_fun <- function(Wi, c0, c1, gamma, tau){ 
#   Xi_fun(Wi, c0, c1, gamma) - (tau)
# }
# Wmin_num_fun <- function(interval, c0, c1, gamma, tau, tol, check.conv){
#     uniroot(f = root_Wmin_fun, 
#                  interval = interval, 
#                  c0 = c0, 
#                  c1 = c1, 
#                  gamma = gamma,
#                  tau = tau,
#                  check.conv = TRUE,
#                  tol = tol)
# }
# 
# #defining optimal avoidance (outside equilibrium) 
# Ai_fun <- function(Wi, c0, c1, gamma, p, rhol, rhov){((1 - p - (rhol*rhov))/(p*(1 - p)))*Xi_fun(Wi, c0, c1, gamma)}
# # Ai_fun(Wi = 1000, c0 = .25, c1 = 0, gamma = .00003, p = .19, rhol = rhol, rhov = rhov)
# 

#defining pstar as a function of W after-tax income (Xi)
pstar_funW <- function(Wi, c0, c1, gamma, rhol, rhov){
  1-((Wi-2*Xi_fun(Wi, c0, c1, gamma) +
        sqrt(((Wi-2*Xi_fun(Wi, c0, c1, gamma))^2) + (4*(rhol*rhov)*(Wi-Xi_fun(Wi, c0, c1, gamma))*Xi_fun(Wi, c0, c1, gamma))))/
       (2*(Wi-Xi_fun(Wi, c0, c1, gamma)))
  )
}

# #defining pstar as a function of A 
p_funA <- function(A, c0, c1, gamma, rhol, rhov){
1-((A-Xi_fun(Tinv_fun(A, c0, c1, gamma), c0, c1, gamma) + 
      sqrt(((A-Xi_fun(Tinv_fun(A, c0, c1, gamma), c0, c1, gamma))^2)+4*(rhol*rhov)*A*Xi_fun(Tinv_fun(A, c0, c1, gamma), c0, c1, gamma)))/(2*A))

}

# #defining optimal avoidance (in equilibrium) 
# Ai_eq_fun <- function(Wi, c0, c1, gamma, rhol, rhov){
#   ((1 - (1-((Wi-2*Xi_fun(Wi, c0, c1, gamma) +
#                sqrt(((Wi-2*Xi_fun(Wi, c0, c1, gamma))^2) + (4*(rhol*rhov)*(Wi-Xi_fun(Wi, c0, c1, gamma))*Xi_fun(Wi, c0, c1, gamma))))/
#               (2*(Wi-Xi_fun(Wi, c0, c1, gamma)))
#   )) - (rhol*rhov))/((1-((Wi-2*Xi_fun(Wi, c0, c1, gamma) +
#                             sqrt(((Wi-2*Xi_fun(Wi, c0, c1, gamma))^2) + (4*(rhol*rhov)*(Wi-Xi_fun(Wi, c0, c1, gamma))*Xi_fun(Wi, c0, c1, gamma))))/
#                            (2*(Wi-Xi_fun(Wi, c0, c1, gamma)))
#   ))*(1 - (1-((Wi-2*Xi_fun(Wi, c0, c1, gamma) +
#                  sqrt(((Wi-2*Xi_fun(Wi, c0, c1, gamma))^2) + (4*(rhol*rhov)*(Wi-Xi_fun(Wi, c0, c1, gamma))*Xi_fun(Wi, c0, c1, gamma))))/
#                 (2*(Wi-Xi_fun(Wi, c0, c1, gamma)))
#   )))))*Xi_fun(Wi, c0, c1, gamma)}
# # Ai_fun(Wi = 1000, c0 = .25, c1 = 0, gamma = .00003, p = .19, rhol = rhol, rhov = rhov)
# 
# 
# #defining EU for truthful reporters 
# EU_A_0_fun <- function(Wi, c0, c1, gamma){log(Xi_fun(Wi, c0, c1, gamma))}
# # EU_A_0_fun(Wi = 1000, c0 = .25, c1 = 0)-log(750)
# 
# ###defining EU for wealth constrained taxpayer
# EU_A_WC_fun <- function(Wi, c0, c1, gamma, rhol, rhov, tau){((rhol*rhov)*
#                                         log(Xi_fun(Wi, c0, c1, gamma) - tau)) + 
#     ((1-(rhol*rhov))*
#        (log(Xi_fun(Wi, c0, c1, gamma) - 
#               tau +
#               T_fun(Wi, c0, c1, gamma))))}
# # EU_A_WC_fun(Wi = 10000, c0 = .25, c1 = 0, tau = tau)
# 
# #defining D''_0
# Dpp0_fun <- function(Wi, c0, c1, gamma, rhol, rhov, tau){
#   (EU_A_0_fun(Wi, c0, c1, gamma) - EU_A_WC_fun(Wi, c0, c1, gamma, rhol, rhov, tau))
# }
# 
# #defining quadratic loss function to perform optimization (instead of root searching)
# Dpp0_se_fun <- function(Wi, c0, c1, gamma, rhol, rhov, tau){
#   (EU_A_0_fun(Wi, c0, c1, gamma) - EU_A_WC_fun(Wi, c0, c1, gamma, rhol, rhov, tau))^2
# }
# 
# #defining the error function to identify W1
# W1_root_fun <- function(Wi, c0, c1, gamma, tau, rhol, rhov){
#   T_fun(Wi, c0, c1, gamma) - (tau/pstar_funW(Wi, c0, c1, gamma, rhol, rhov)) 
# }
# 
# #defining the error function to identify W1
# W1_se_fun <- function(Wi, c0, c1, gamma, tau, rhol, rhov){
#   (T_fun(Wi, c0, c1, gamma) - (tau/pstar_funW(Wi, c0, c1, gamma, rhol, rhov)))^2
# }
# 
# #function to compute total avoidance in the economic system for a linear tax system
# #(derivation in Mathematica file analytic.nb)
# A_lin_fun <- function(W0, c0){
#   17105.3769312732*c0*(1 + erf(12.33572806744155 - (1.1609233137739046*log(W0))))
# }
# #function to compute revenues collected from avoidance enforcement for a linear tax system
# #(derivation in Mathematica file analytic.nb)
# RA_lin_fun <- function(W0, c0, rhol, rhov){
#   17105.3769312732*c0*(rhol*rhov)*(1 + erf(12.33572806744155 - (1.1609233137739046*log(W0))))
# }
# #function to compute revenues collected from compliance (honest reporting) for a linear tax system
# #(derivation in Mathematica file analytic.nb)
# RC_lin_fun <- function(W0, c0){
#   17105.3769312732*c0*erfc(
#     12.33572806744155 - 1.1609233137739046*log(W0))
# }
# #function to compute total revenues collected for a linear tax system
# #(derivation in Mathematica file analytic.nb)
# R_lin_fun <- function(W0, c0, rhol, rhov){
#   RC_lin_fun(W0, c0) + RA_lin_fun(W0, c0, rhol, rhov)
# }
# # #checks
# # R_lin_fun_check <- function(W0, c0, rhol, rhov){
# #   17105.3769312732*c0*(rhol*rhov)*(1 + 
# #     erf(12.33572806744155 - 1.1609233137739046*log(W0))) + 
# #   17105.3769312732*c0*erfc(
# #     12.33572806744155 - 1.1609233137739046*log(W0))
# # }
# # A_lin_fun(W0 = 20000, c0 = .2)
# # RA_lin_fun(W0 = 20000, c0 = .2, rhol = rhol, rhov = rhov)
# # RC_lin_fun(W0 = 20000, c0 = .2)
# # R_lin_fun(W0 = 20000, c0 = .2, rhol = rhol, rhov = rhov)
# # R_lin_fun_check(W0 = 20000, c0 = .2, rhol = rhol, rhov = rhov)
# # #there are small differences in A2/A_lin wrt M, but the plots are the same
# # plot(100:100000, A_lin_fun(100:100000, c0 = .2, rhol = rhol, rhov = rhov))
# #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# #@@@ defining firms'expected profits @@@
# 
# #defining the error function to identify W2
# W2_root_fun <- function(Wi, c0, c1, gamma, rhol, rhov, tau){
#   (
#     Ai_eq_fun(Wi, c0, c1, gamma, rhol, rhov) *
#       pstar_funW(Wi, c0, c1, gamma, rhol, rhov)
#   ) - tau 
# }
# 
# #defining function to find W2 as a numerical root
# W2_rootfinding_fun <- function(c0, c1, gamma, rhol, rhov, tau){
#   #computing the minimum W to have a non-negative W^u
#   min_W_temp <-  Wmin_analyt_fun(c0 = c0, c1 = c1, gamma = gamma, tau = tau)
#   #finding W1 as root 
#   uniroot(f = W2_root_fun, 
#           interval = c(min_W_temp+.1, 2000000), 
#           c0 = c0, 
#           c1 = c1,
#           gamma = gamma,
#           rhol = rhol, rhov = rhov,
#           tau = tau,
#           check.conv = check_conv,
#           tol = tol)
# }
# # W2_rootfinding_fun(c0 = .28021974531250876,
# #                          c1 = 0,
# #                          gamma = 10,
# #                          rhol = .35,
# #                          rhov = .85,
# #                         tau = 4800
# #                         )$root
# 
# #defining the error function to identify W2
# W2_se_fun <- function(Wi, c0, c1, gamma, rhol, rhov, tau){
#   ((
#     Ai_eq_fun(Wi, c0, c1, gamma, rhol, rhov) *
#       pstar_funW(Wi, c0, c1, gamma, rhol, rhov)
#   ) - tau )^2
# }
# 
# #defining function to find W2 as a numerical min of se
# W2_minsefinding_fun <- function(c0, c1, gamma, rhol, rhov, tau){
#   #computing the minimum W to have a non-negative W^u
#   min_W_temp <-  Wmin_analyt_fun(c0 = c0, c1 = c1, gamma = gamma, tau = tau)
#   #finding W1 as root 
#   optimize(f = W2_se_fun, 
#           interval = c(min_W_temp+.1, 2000000), 
#           c0 = c0, 
#           c1 = c1,
#           gamma = gamma,
#           rhol = rhol, rhov = rhov,
#           tau = tau,
#           tol = tol)
# }
# # W2_rootfinding_fun(c0 = .28021974531250876,
# #                          c1 = 0,
# #                          gamma = 10,
# #                          rhol = .35,
# #                          rhov = .85,
# #                         tau = 4800
# #                         )$root
# 
# #defining integrand function of firms' expected profits
# pi_integrand_fun <- function(Wi, c0, c1, rhol, rhov, gamma, mu, sigma){
#   pstar_funW(Wi, c0, c1, gamma, rhol, rhov) * T_fun(Wi, c0, c1, gamma) * dlnorm(x = Wi, meanlog = mu, sdlog = sigma)
# }
# #defining function of firms' expected profits as a function of W2
# pi_of_W2_fun <- function(c0, c1, gamma, rhol, rhov, mu, sigma, tau, n, W2, cl, vi){
#   ((1/n)*integrate(pi_integrand_fun,
#                    lower = W2, upper = wbar,
#                    c0 = c0,
#                    c1 = c1,
#                    gamma = gamma, 
#                    rhol = rhol,
#                    rhov = rhov,
#                    mu = mu, 
#                    sigma = sigma)$value) - 
#     ((tau/n)*(1-plnorm(q = W2, meanlog = mu, sdlog = sigma))) -
#     (rhol*cl) -
#     vi
# }
# # pi_of_W2_fun(c0 = .28021974531250876,
# #        c1 = 0,
# #        gamma = 5,
# #        rhol = .35,
# #        rhov = .85,
# #        mu = 10.2548,
# #        sigma = 0.60909,
# #        tau = 4800,
# #        n = 5,
# #        W2 = 72144.9,
# #        cl = 10000,
# #        vi = 1000)
# 
# #defining function of firms' expected profits 
# pi_fun_rootW2 <- function(c0, c1, gamma, rhol, rhov, mu, sigma, tau, n, cl, vi, lean = TRUE){
#   
#   W2_uniroot <- W2_rootfinding_fun(c0, c1, gamma, rhol, rhov, tau)
#   if(lean){
#   pi_of_W2_fun(W2 = W2_uniroot$root, c0, c1, gamma, rhol, rhov, mu, sigma, tau, n, cl, vi)
#   }else{
#     c(pi_of_W2_fun(W2 = W2_uniroot$root, c0, c1, gamma, rhol, rhov, mu, sigma, tau, n, cl, vi),
#       W2_uniroot$f.root)
#   }
# }
# 
# 
# # pi_fun_root(c0 = .28021974531250876,
# #        c1 = 0,
# #        gamma = 5,
# #        rhol = .35,
# #        rhov = .85,
# #        mu = 10.2548,
# #        sigma = 0.60909,
# #        tau = 4800,
# #        n = 1,
# #        cl = 3000,
# #        vi = 500)
# 
# #defining function of difference of expected profits with respect to a given value as a function of rhov 
# pi_root_fun_rootW2 <- function(rhov, rhol, value, c0, c1, gamma, mu, sigma, tau, n, cl, vi, lean = TRUE){
#   pi_fun_rootW2(c0, c1, gamma, rhol, rhov, mu, sigma, tau, n, cl, vi, lean = lean)-value
# }
# pi_root_Vfun_rootW2 <- Vectorize(pi_root_fun_rootW2, vectorize.args = "rhov")
# rhov_loop <- seq(.01, .99, by =.01)
# # plot(rhov_loop, 
# #      pi_root_Vfun_rootW2(c0 = 0,
# #                   c1 = .45,
# #                   gamma = .00004,
# #                   rhol = .35,
# #                   value = 0,
# #                   rhov = rhov_loop,
# #                   mu = 10.2548,
# #                   sigma = 0.60909,
# #                   tau = 4800,
# #                   n = 1,
# #                   cl = 3000,
# #                   vi = 500)
# # )
# 
# #defining function to find rhov such that the expected profits are equal to a given value
# pi_rootfinding_fun_rootW2 <- function(rhov, rhol, value, c0, c1, gamma, mu, sigma, tau, n, cl, vi, lean = TRUE){
# 
#   #finding rhov as root 
#   uniroot(f = pi_root_fun_rootW2, 
#           interval = c(.000001, 2), 
#           rhol = rhol,
#           value = value,
#           c0 = c0, 
#           c1 = c1,
#           gamma = gamma,
#           mu = mu,
#           sigma = sigma,
#           tau = tau,
#           n = n,
#           cl = cl,
#           vi = vi,
#           check.conv = check_conv,
#           tol = tol)
# }
# 
# 
# 
# 
# 
# #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# 
# pi_fun_minseW2 <- function(c0, c1, gamma, rhol, rhov, mu, sigma, tau, n, cl, vi, lean = TRUE){
#   
#   W2_minse <- W2_minsefinding_fun(c0, c1, gamma, rhol, rhov, tau)
#   if(lean){
#     pi_of_W2_fun(W2 = W2_minse$minimum, c0, c1, gamma, rhol, rhov, mu, sigma, tau, n, cl, vi)
#   }else{
#     c(pi_of_W2_fun(W2 = W2_minse$minimum, c0, c1, gamma, rhol, rhov, mu, sigma, tau, n, cl, vi),
#       W2_minse$f.root)
#   }
# }
# 
# 
# # pi_fun_root(c0 = .28021974531250876,
# #        c1 = 0,
# #        gamma = 5,
# #        rhol = .35,
# #        rhov = .85,
# #        mu = 10.2548,
# #        sigma = 0.60909,
# #        tau = 4800,
# #        n = 1,
# #        cl = 3000,
# #        vi = 500)
# 
# #defining function of difference of expected profits with respect to a given value as a function of rhov 
# pi_minse_fun_minseW2 <- function(rhov, rhol, value, c0, c1, gamma, mu, sigma, tau, n, cl, vi, lean = TRUE){
#   (pi_fun_minseW2(c0, c1, gamma, rhol, rhov, mu, sigma, tau, n, cl, vi, lean = lean)-value)^2
# }
# pi_minse_Vfun_minseW2 <- Vectorize(pi_minse_fun_minseW2, vectorize.args = "rhov")
# rhov_loop <- seq(.01, .99, by =.01)
# # plot(rhov_loop, 
# #      pi_minse_Vfun_minseW2(c0 = 0,
# #                   c1 = .45,
# #                   gamma = .00004,
# #                   rhol = .35,
# #                   value = 0,
# #                   rhov = rhov_loop,
# #                   mu = 10.2548,
# #                   sigma = 0.60909,
# #                   tau = 4800,
# #                   n = 1,
# #                   cl = 3000,
# #                   vi = 500)
# # )
# 
# #defining function to find rhov such that the expected profits are equal to a given value
# pi_rootfinding_fun_rootW2 <- function(rhov, rhol, value, c0, c1, gamma, mu, sigma, tau, n, cl, vi, lean = TRUE){
#   
#   #finding rhov as root 
#   optimize(f = pi_minse_fun_minseW2, 
#           interval = c(.000001, 2), 
#           rhol = rhol,
#           value = value,
#           c0 = c0, 
#           c1 = c1,
#           gamma = gamma,
#           mu = mu,
#           sigma = sigma,
#           tau = tau,
#           n = n,
#           cl = cl,
#           vi = vi,
#           tol = tol)$minimum
# }
# 
# 
# 
# 
# 
# 



