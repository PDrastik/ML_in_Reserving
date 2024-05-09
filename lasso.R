library(data.table)
library(glmnet)
library(ggplot2)
library(patchwork)
library(DT)


# other setup
options("scipen"=99) # turn off scientific notation

# colour palette for some plots
# IFoA colours
# primary colours
dblue <- "#113458"
  mblue <- "#4096b8"
    gold <- "#d9ab16"
      lgrey <- "#dcddd9"
        dgrey <- "#3f4548"
          black <- "#3F4548"
            #secondary colours
          red <- "#d01e45"
            purple <- "#8f4693"
              orange <- "#ee741d"
                fuscia <- "#e9458c"
                  violet <- "#8076cf"
                    
                  LinearSpline <- function(var, start, stop){
                    pmin(stop - start, pmax(0, var - start))
                  }
                  
                  
CreateSyntheticData<-function(whichsim, numperiods)
{
  
  # whichsim determins which data set to simulate
  # numperiods determines the size of the triangle
  
  # create the acc/dev/cal parameters
  kk <- rep(1:numperiods, each = numperiods) #AQ
  jj <- rep(1:numperiods, times= numperiods) #DQ
  tt <- kk+jj-1 # PQ
  
  # set alpha/beta/gamma - hard-code up the sim values
  if (whichsim == 1){
    alpha <- log(100000)+0.1*LinearSpline(kk,1,15)+0.2*LinearSpline(kk,15,20) - 0.05*LinearSpline(kk,30,40) 
    beta  <- (16/3 - 1)*log(jj)- (1/3)*jj
    gamma <- 0
    mu <- exp( alpha + beta + gamma)  
  }
  else if (whichsim == 2){
    alpha <- log(100000)+0.1*LinearSpline(kk,1,15)+0.2*LinearSpline(kk,15,20) - 0.05*LinearSpline(kk,30,40) 
    beta  <- (16/3 - 1)*log(jj)- (1/3)*jj  # a is 16/3, b is 1/3 
    gamma <- gammafunc(tt)
    mu <- exp( alpha + beta + gamma)  
  }
  else if (whichsim == 3){
    alpha <- log(100000)+0.1*LinearSpline(kk,1,15)+0.2*LinearSpline(kk,15,20) - 0.05*LinearSpline(kk,30,40) 
    beta  <- (16/3 - 1)*log(jj)- (1/3)*jj  # a is 16/3, b is 1/3 
    gamma <- gammafunc(tt)
    mu <- exp( alpha + beta + gamma + 0.3*beta*ifelse(kk>16 & jj>20,1,0))  
  }
  else if (whichsim == 4){
    alpha <- log(100000)+0.1*LinearSpline(kk,1,15)+0.2*LinearSpline(kk,15,20) - 0.05*LinearSpline(kk,30,40) 
    beta  <- (16/3 - 1)*log(jj)- (1/3)*jj  # a is 16/3, b is 1/3 
    gamma <- gammafunc(tt)
    mu <- exp( alpha + beta + gamma*((numperiods-1)-LinearSpline(jj,1,numperiods))/(numperiods-1) )  # need to check
  }
  
  varbase <- (0.3 * mu[  kk==1 & jj ==16] )^2 # can scale variance up and down here
  CC  <-  varbase / mu[  kk==1 & jj ==16]
  
  vars   <- CC*mu
  tausq  <- log (vars / (mu^2) + 1)
  
  pmts <- exp( rnorm( numperiods^2, mean = log(mu)-0.5*tausq , sd = sqrt(tausq)  ) )
  
  # indicator for past/future = traint/test
  train_ind<-(tt<=numperiods)
  
  ### data.table for output
  full<-data.table(pmts, acc=as.integer(kk), dev=as.integer(jj), cal=as.integer(tt), mu, train_ind )
  full
}


#---------------------------
# function to generate calendar period effects used in CreateSyntheticData()
# written as a seperate function for convenience

gammafunc <- function(t){
  gg <- 
    ifelse( t<=12, gg <- 0.0075*LinearSpline(t,1,12),
            ifelse(t<=24,  gg <- 0.0075*LinearSpline(12,1,12) + 0.001* (t-12)*(t-11)/2,
                   ifelse(t<=32, gg <- 0.0075*LinearSpline(12,1,12) + 0.001* (24-12)*(24-11)/2,
                          ifelse(t<=40, gg <- 0.0075*LinearSpline(12,1,12) + 0.001* (24-12)*(24-11)/2 + 0.002*(t-32)*(t-31)/2,
                                 0.0075*LinearSpline(12,1,12) + 0.001* (24-12)*(24-11)/2 + 0.002*(40-32)*(40-31)/2
                          ))))
  gg  
}                  

use_data_set <- 3       # which data set to use
use_data_set_seed <- 130  # seed to generate data
num_periods <- 40   # size of data set

set.seed(use_data_set_seed)

dat <- CreateSyntheticData(whichsim=use_data_set, numperiods=num_periods)

# head(dat,10) gives the same information, the rest is just formatting

head(dat, 10)  |> 
  datatable() |> 
  formatRound(c("pmts", "mu"), digits = 0)

# get limits for use with raster plots
data_raster_limits <- c(floor(dat[, min(log(pmts))]), ceiling(dat[, max(log(pmts))]))

ggplot(dat, aes(dev, acc)) +
  geom_raster(aes(fill = log(pmts)))+
  geom_line(aes(x=num_periods+1-acc, y=acc), colour=dgrey, size=2)+
  scale_y_reverse()+
  scale_fill_viridis_c(begin=1, end=0, limits=data_raster_limits)+
  theme_classic()+
  labs(x="Development quarter", y="Accident quarter", title="Log(payments)")+
  NULL

# function to calculate scaling factors for the basis functions
# scaling is discussed in the paper
GetScaling <- function(vec) {
  fn <- length(vec)
  fm <- mean(vec)
  fc <- vec - fm
  rho_factor <- ((sum(fc^2))/fn)^0.5
}

# function to create the ramps for a particular primary vector
GetRamps <- function(vec, vecname, np, scaling){
  
  # vec = fundamental regressor
  # vecname = name of regressor
  # np = number of periods
  # scaling = scaling factor to use
  
  # pre-allocate the matrix to hold the results for speed/efficiency
  n <- length(vec)
  nramps <- (np-1)
  
  mat <- matrix(data=NA, nrow=n, ncol=nramps)
  cnames <- vector(mode="character", length=nramps)
  
  
  col_indx <- 0
  
  for (i in 1:(np-1)){
    col_indx <- col_indx + 1
    
    mat[, col_indx] <- LinearSpline(vec, i, 999) / scaling
    cnames[col_indx] <- paste0("L_", i, "_999_", vecname)
  }
  
  colnames(mat) <- cnames
  
  return(mat)
}

# create the step (heaviside) function interactions
GetInts <- function(vec1, vec2, vecname1, vecname2, np, scaling1, scaling2) {
  
  # vec1 = fundamental regressor 1
  # vec2 = fundamental regressor 2
  # vecname1 = name of regressor 1
  # vecname2 = name of regressor 2
  # np = number of periods
  # scaling1 = scaling factor to use for regressor 1
  # scaling2 = scaling factor to use for regressor 2
  
  
  # pre-allocate the matrix to hold the results for speed/efficiency
  n <- length(vec1)
  nints <- (np-1)*(np-1)
  
  mat <- matrix(data=NA_real_, nrow=n, ncol=nints)
  cnames <- vector(mode="character", length=nints)
  
  
  col_indx <- 0
  
  for (i in 2:np){
    
    ivec <- LinearSpline(vec1, i-1, i) / scaling1
    iname <- paste0("I_", vecname1, "_ge_", i)
    
    if (length(ivec[is.na(ivec)]>0)) print(paste("NAs in ivec for", i))
    
    for (j in 2:np){
      col_indx <- col_indx + 1  
      mat[, col_indx] <- ivec * LinearSpline(vec2, j-1, j) / scaling2
      cnames[col_indx] <- paste0(iname, "xI_", vecname2, "_ge_", j)
      
      jvec <- LinearSpline(vec2, j-1, j) / scaling2
      if (length(jvec[is.na(jvec)]>0)) print(paste("NAs in jvec for", j))
      
    }
  }
  
  colnames(mat) <- cnames
  
  return(mat)
  
  
}

# get the scaling values
rho_factor_list <- vector(mode="list", length=3)
names(rho_factor_list) <- c("acc", "dev", "cal")

for (v in c("acc", "dev", "cal")){
  # NB: only calculating scaling using past data    
  rho_factor_list[[v]] <- GetScaling(dat[train_ind == TRUE, get(v)])
}


# main effects - matrix of values of Ramp functions

main_effects_acc <- GetRamps(vec = dat[, acc], vecname = "acc", np = num_periods, scaling = rho_factor_list[["acc"]])
main_effects_dev <- GetRamps(vec = dat[, dev], vecname = "dev", np = num_periods, scaling = rho_factor_list[["dev"]])
main_effects_cal <- GetRamps(vec = dat[, cal], vecname = "cal", np = num_periods, scaling = rho_factor_list[["cal"]])

main_effects <- cbind(main_effects_acc, main_effects_dev, main_effects_cal)


# interaction effects
int_effects <- cbind(
  GetInts(vec1=dat[, acc], vecname1="acc", scaling1=rho_factor_list[["acc"]], np=num_periods, 
          vec2=dat[, dev], vecname2="dev", scaling2=rho_factor_list[["dev"]]),
  
  GetInts(vec1=dat[, dev], vecname1="dev", scaling1=rho_factor_list[["dev"]], np=num_periods, 
          vec2=dat[, cal], vecname2="cal", scaling2=rho_factor_list[["cal"]]),
  
  GetInts(vec1=dat[, acc], vecname1="acc", scaling1=rho_factor_list[["acc"]], np=num_periods, 
          vec2=dat[, cal], vecname2="cal", scaling2=rho_factor_list[["cal"]])
)


varset <- cbind(main_effects, int_effects)

varset[1:10, 39:44]

rho_factor_list

# drop any constant columns over the training data set
# do this by identifying the constant columns and dropping them

# get the past data subset only using TRUE/FALSE from train_ind in dat
varset_train <- varset[dat$train_ind, ]

# identify constant columns as those with max=min
rm_cols <- varset_train[, apply(varset_train, MARGIN=2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))]

# drop these constant columns
varset <- varset[, !(colnames(varset) %in% colnames(rm_cols))]

my_pmax <- num_periods^2   # max number of variables ever to be nonzero
my_dfmax <- num_periods*10  #max number of vars in the model

time1 <- Sys.time()
cv_fit <- cv.glmnet(x = varset[dat$train_ind,], 
                    y = dat[train_ind==TRUE, pmts], 
                    family = "poisson", 
                    nlambda = 200,  # default is 100 but often not enough for this type of problem
                    nfolds = 8,
                    thresh = 1e-08, # default is 1e-7 but we found we needed to reduce it from time to time
                    lambda.min.ratio = 0, # allow very small values in lambda
                    dfmax = my_dfmax, 
                    pmax = my_pmax, 
                    alpha = 1, 
                    standardize = FALSE, 
                    maxit = 200000)  # convergence can be slow so increase max number of iterations
print("time taken for cross validation fit: ")
## [1] "time taken for cross validation fit: "
Sys.time() - time1
## Time difference of 54.57403 secs

plot(cv_fit)

# all coefficients, including those that are 0
coefs_min <- predict(cv_fit, type = "coefficients", s = cv_fit$lambda.min)  # NB this is a data.frame not a vector
coefnames <- c("Intercept", colnames(varset))  # don't forget the intercept  # coeff names - this is a vector

# get indicators for non-zero ones
ind_nz_min<-which(!(coefs_min == 0))

# make a data.table for easy viewing
nzcoefs_min <- data.table(Parameter = coefnames[ind_nz_min], Coefficient = coefs_min[ind_nz_min,])

# print the table
datatable(nzcoefs_min) |> 
  formatRound("Coefficient", digits = 4)

dat[, fitted := as.vector(predict(cv_fit, 
                                  newx = varset, 
                                  s = cv_fit$lambda.min, 
                                  type="response"))]

GraphModelVals<-function(dat, primary_predictor, secondary_predictor, secondary_predictor_val, 
                         xaxis_label, yaxis_label, var_names, log_values = TRUE, 
                         include_st=FALSE, include_legend=FALSE, font_size=6){
  
  # dat = input data. Must be in wide format
  # primary_predictor = plot values for this predictor
  # secondary_predictor = hold this predictor fixed
  # secondary_predictor_val = value to hold secondary_predictor fixed at
  # xaxis_label = label for x axis
  # yaxis_label = label for y axis
  # var_names = names of actual / fitted variables in the input data to plot.     
  #    var_names must be list with names like this: list(actual="pmts", mean="mu", fitted="fitted")
  # log_values = plot log(values) - default = TRUE
  # include_st = include a subtitle in the plot saying the grey rectangle is past data, default=FALSE    
  # include_legend = include the legend in the plot, default=FALSE    
  # font_size = size of font, default = 6
  # (these last 3 variables + default values might seem odd - this function is used in a later article
  #    and the defaults are sensible in that context)    
  
  
  # extract data we want to use
  use_dat <- dat[get(secondary_predictor) == secondary_predictor_val, ]
  
  # turn into long format (tidy format in this case) using melt.data.table since that works better with ggplot
  dat_long <- melt(dat[get(secondary_predictor) == secondary_predictor_val, ],
                   measure.vars = unlist(var_names),
                   id.vars = primary_predictor)
  
  # make the names nicer - colnames and labels
  setnames(dat_long, primary_predictor, "predictor")
  
  dat_long[variable == var_names$actual, variable := "Simulated"
  ][variable == var_names$mean, variable := "Underlying"
  ][variable == var_names$fitted, variable := "Fitted"]
  
  # get the levels of the variables right so that they are plotted in the right order
  dat_long[, variable := factor(variable, levels=c("Fitted", "Simulated", "Underlying"))]
  
  
  if (log_values) dat_long[, value := log(value)]
  
  # figure out past data rectangle coordinates
  xmin1 <- use_dat[train_ind == TRUE, min(get(primary_predictor))]
  xmax1 <- use_dat[train_ind == TRUE, max(get(primary_predictor))]
  
  ymin1 <- dat_long[, min(value)]*0.95
  ymax1 <- dat_long[, max(value)]*1.05
  
  
  # draw the tracking plots
  g <- ggplot(data=dat_long, aes(x=predictor, y=value, group=variable))+
    geom_line(aes(linetype=variable, colour=variable, size=variable, alpha=variable))+
    geom_line(aes(linetype=variable, colour=variable))+
    scale_colour_manual(name="", values=c(red, dgrey, dgrey))+
    scale_linetype_manual(name="", values=c("solid", "solid", "dotted"))+
    scale_size_manual(name="", values=c(2,1,1))+
    scale_alpha_manual(name="", values=c(0.8, 0.5, 0.5))+
    theme_classic()+
    annotate(geom="rect", xmin=xmin1, xmax=xmax1, ymin=ymin1, ymax=ymax1, alpha=0.1)+
    labs(x=xaxis_label, y=yaxis_label, title=paste(xaxis_label, "tracking for", secondary_predictor, "=", secondary_predictor_val)) +
    theme(axis.title = element_text(size = font_size), axis.text = element_text(size = font_size-1))
  
  if(include_st==TRUE) g <- g + labs(subtitle="Past data in grey rectangle") + theme(plot.subtitle = element_text (size = font_size))
  
  g <- if(include_legend==TRUE) g + theme(legend.position="bottom") else g + theme(legend.position = "none")
  
  
  
  # return the results  
  invisible(list(data=dat_long, graph=g))
}

dev_graph_list <- GraphModelVals(dat, 
                                 primary_predictor = "dev", 
                                 secondary_predictor = "acc", 
                                 secondary_predictor_val = 20, 
                                 xaxis_label = "Development quarter",
                                 yaxis_label = "Log(Payments)",
                                 var_names = list(actual="pmts", mean="mu", fitted="fitted"),
                                 include_st = TRUE,
                                 include_legend = TRUE,
                                 font_size = 10)

dev_graph_list$graph

acc_graph_list <- GraphModelVals(dat, 
                                 primary_predictor = "acc", 
                                 secondary_predictor = "dev", 
                                 secondary_predictor_val = 24, 
                                 xaxis_label = "Accident quarter",
                                 yaxis_label = "Log(Payments)",
                                 var_names = list(actual="pmts", mean="mu", fitted="fitted"),
                                 include_st = TRUE,
                                 include_legend = TRUE,
                                 font_size = 10)
acc_graph_list$graph

# heat maps

GraphHeatMap <- function(dat, x="dev", y="acc", actual, fitted, lims=c(0.25, 4),
                         xlab="Development quarter", ylab="Accident Quarter"){
  
  # copy data to avoid modifying original
  localdat <- copy(dat)
  
  # get fails if there is a variable with the same name so make local copies
  local_x <- x
  local_y <- y
  local_actual <- actual
  local_fitted <- fitted
  
  # make restricted Avs F for heatmap and set up past/future split line
  np <- max(localdat[[y]])
  
  localdat[, .avsf := get(local_actual) / get(local_fitted)
  ][, .avsf_restrict_log := log(pmax(min(lims), pmin(max(lims), .avsf)))
  ][, .past_line := np + 1 - get(local_y)]
  
  
  g <- ggplot(data=localdat, aes_string(x=local_x, y=local_y)) +
    geom_tile(aes(fill = .avsf_restrict_log))+scale_y_reverse()+
    theme_classic()+
    scale_fill_gradient2(name="AvF_min", low=mblue, mid="white", high=red, midpoint=0, 
                         space="Lab", na.value=lgrey, guide="colourbar")+
    labs(x=xlab, y=ylab)+
    geom_line(aes_string(x=".past_line", y=local_y), colour=dgrey, size=2)+
    theme(strip.text = element_text(size=8,colour=dgrey), 
          strip.background = element_rect(colour="white", fill="white"))+
    theme(axis.title.x = element_text(size=10), axis.text.x  = element_text(size=10))+
    theme(axis.title.y = element_text(size=10), axis.text.y  = element_text(size=10))+
    theme(element_line(size=0.25, colour=dgrey))+
    theme(legend.position="none", )+  
    NULL    
  
  
  invisible(list(data=localdat, graph=g))
  
  
}

g <- GraphHeatMap(dat, x="dev", y="acc", actual="pmts", fitted="fitted")
g$graph

os_acc <- dat[train_ind == FALSE, 
              .(LASSO = sum(fitted), simulated = sum(pmts), underlying = sum(mu)), 
              keyby=.(acc)]

os <- os_acc[, .(LASSO = sum(LASSO), 
                 simulated = sum(simulated), 
                 underlying = sum(underlying))]

# get cumulative payments to make it easier to calculate CL factors
dat[, cumpmts := cumsum(pmts), by=.(acc)][train_ind==FALSE, cumpmts := NA]

# 8-period average
cl_fac <- numeric(num_periods-1) # hold CL factors

for (j in 1:num_periods-1){
  
  cl_fac[j] <- 
    dat[train_ind==TRUE & dev == (j+1) & acc > (num_periods-8-j) & acc <= (num_periods-j), sum(cumpmts)] /
    dat[train_ind==TRUE & dev == (j) & acc > (num_periods-8-j) & acc <= (num_periods-j), sum(cumpmts)]
}

# accumulate the CL factors
cl_cum <- cumprod(rev(cl_fac))

# leading diagonal for projection
leading_diagonal <- dat[train_ind==TRUE & cal == num_periods & acc > 1, cumpmts]

# CL amounts now
cl_os <- cl_cum * leading_diagonal - leading_diagonal

os_acc[, Chainladder := cl_os]

os[, Chainladder := sum(cl_os)]

# make a long [tidy format] version of the data for use with ggplot2
os_acc_long <- melt(os_acc, id.vars = "acc")

# divide by Bn to make numbers more readable
os_acc_long[, value := value/1e9]

os_plot <-
  ggplot(data=os_acc_long, aes(x=acc, y=value, colour=variable, 
                               linetype=variable, size=variable, alpha=variable))+
  geom_line()+
  scale_linetype_manual(name="", values=c("solid", "dashed", "dotted", "solid" ))+
  scale_colour_manual(name="", values=c(red, dgrey, dgrey, mblue ))+
  scale_size_manual(name="", values=c(2, 1, 1, 1.5))+
  scale_alpha_manual(name="", values=c(0.8, 0.5, 0.5, 0.8))+
  coord_cartesian(ylim=c(0, 40))+
  theme_classic()+
  theme(legend.position="bottom", legend.title=element_blank())+
  labs(x="Accident quarter", y="Amount (B)")+
  annotate(geom="text", x=37, y=40, label="488B->")+
  ggtitle("Outstanding amounts")

os_plot

cols <- setdiff(names(os_acc), "acc")

os_acc[, lapply(.SD, function(x) x/1e9)][, acc := acc*1e9]  |> 
  datatable() |> 
  formatRound(cols, digits = 1)

os[, lapply(.SD, function(x) x/1e9)] |> 
  datatable() |> 
  formatRound(cols, digits = 1)

