library(tidyverse)
library(lubridate)
library(forecast)
library(vars) ## loads urca, MASS
library(tsDyn) ## loads vars
library(xtable)
library(tseries)
library(gridExtra)

companion <- function (x) {
        ## Create companion matrix of vec2var object
        
        K <- x$K
        p <- x$p
        A <- unlist(x$A)
        companion <- matrix(0, nrow = K * p, ncol = K * p)
        companion[1:K, 1:(K * p)] <- A
        if (p > 1) {
                j <- 0
                for (i in (K + 1):(K * p)) {
                        j <- j + 1
                        companion[i, j] <- 1
                }
        }
        return(companion)
}

## ---------------------------------------------------------------------------- ##

back_tran <- function(y, cap = 1000){(cap*exp(y))/(1 + exp(y))}

hosp01 <- read_csv("Hospital_Stats.csv")
hosp01 <- hosp01 %>% mutate(date = mdy(date),
                            Lcensus= log(Census/(1000 - Census)),
                            DayWeek = wday(date, label = TRUE))


i_t <- read_csv("CRI_Incidence.csv")
i_t <- i_t %>% mutate(date = mdy(Date), Lincidence = log(Incidence), .before = 1) %>% dplyr::select(date, Incidence, Lincidence)

master_data <- full_join(hosp01, i_t) %>% filter(date <= mdy("12-05-2020"))
master_data$Time <- 1:nrow(master_data)

## ---------------------------------------------------------------------------- ##

## Pearson Correlation
blast02_I <- matrix(NA, nrow = 22, ncol = 2)
colnames(blast02_I) <- c("Lag", "Pearson Correlation")
blast02_I[,1] <- 0:-21
for(i in 0:21){blast02_I[i + 1, 2] <- cor(master_data$Lcensus, lag(master_data$Lincidence, i), use = "pairwise.complete.obs", method = "pearson")}
blast02_I <- as_tibble(blast02_I)
blast02_I

## Fig 1.
fig1_data <- master_data %>% mutate(Scaled_C = 100*(Lcensus - min(Lcensus))/(max(Lcensus) - min(Lcensus)), 
                                   Scaled_i = 100*(Lincidence - min(Lincidence))/(max(Lincidence) - min(Lincidence)))
fig1_data$date <- as.Date(fig1_data$date, "%m/%d/%y") 
fig1 <- ggplot(fig1_data, aes(x = date)) 
fig1 <- fig1 + geom_line(aes(y = Scaled_C), color = "blue")
fig1 <- fig1 + geom_line(aes(y = Scaled_i), color = "red")
fig1 <- fig1 + labs(y = "Standardized Scale", x = "Date")
fig1 <- fig1 + scale_x_date(date_breaks = "1 month", limits = as.Date(c("2020-05-15", "2020-12-05")), date_labels = "%b")
fig1 <- fig1 + theme(axis.text = element_text(size = 6),
                     axis.title = element_text(size = 6),
                     panel.background = element_rect(fill = "white"),
                     panel.grid = element_line(colour = "grey85"), 
                     panel.border = element_rect(colour = "black", fill = NA, size = 1))
fig1

## ---------------------------------------------------------------------------- ##

## VECM
dat <- master_data  %>% dplyr::select("Lcensus", "Lincidence") %>% as.matrix()
colnames(dat) <- c("Census", "Incidence")
H1 <- ca.jo(dat, type = "trace", ecdet = "const", spec = "transitory", K = 7, season = 7)

## AIC scores 
VARselect(dat, lag.max = 14, type = "const", season = 7)

## Johansen trace test
summary(H1)

## Maximum eigenvalue modulus
stb <- vec2var(H1, r = 1)
comp <- companion(stb)
Mod(eigen(comp)$values) 

## Likelihood ratio test for linear trend
lttest(H1, r = 1)

## Check stationarity of error correction term
ect <- dat[, 1] - 0.8013363 * dat[, 2] + 7.8266336
plot(ect, type = "l"); abline(h = 0)
## KPSS test
kpss.test(ect, null = c("Level"), lshort = FALSE) 
## Augmented Dickey-Fuller test
trunc(12*(length(ect)/100)^(1/4)) ## 14 chosen as upper bound
tryit01 <- ur.df(ect, type = "trend", lags = 14, selectlags = "BIC")
summary(tryit01) ## Claim of 6 lags
tryit02 <- ur.df(ect, type = "trend", lags = 6)
summary(tryit02) ## tau_3 and Phi_3
tryit03 <- ur.df(ect, type = "drift", lags = 6)
summary(tryit03) ## tau_2 and Phi_1
tryit04 <- ur.df(ect, type = "none", lags = 6)
summary(tryit04) ## tau_1; reject at 10%, but not 5%

## Generate seasonal terms
week("2020-05-15") ## 20
x1 <- ts(dat[,1], start = c(20, 1), frequency = 7) ## Wlog start at 1 to be consistent with cajorls()
seas <- seasonaldummy(x1)
seas[seas == 0] <- -1/7 ## Center columns to have mean 0 within week; simple coding
seas[seas == 1] <- 6/7 ## Center columns to have mean 0 within week; simple coding
vecm2 <- VECM(dat, lag = 6, estim = "ML", r = 1, include = "const", LRinclude = "none", exogen = seas) 
vecm <- VECM(dat, lag = 6, estim = "ML", r = 1, include = "none", LRinclude = "const", exogen = seas)
AIC(vecm); AIC(vecm2) ## AIC ~ -1519, -1516

## Serial correlation test
serial.test(stb) 
## Fig 3. 
cors <- acf(residuals(stb), lag.max = 21)
cors <- tibble(cor = c(cors$acf[, , 1][, 1], cors$acf[, , 1][, 2], cors$acf[, , 2][, 1], cors$acf[, , 2][, 2]),
        lag = rep(0:-21, 4), 
        name = c(rep("A", 22), rep("B", 22), rep("C",22), rep("D",22)))

fig3 <- grid.arrange(ggplot(cors %>% filter(name=="A"), aes(x = lag, y = cor)) +
        geom_segment(aes(x = lag, y = 0, xend = lag, yend = cor)) + 
        geom_hline(yintercept = 2/sqrt(198), lty = 2, colour = "blue4") + 
        geom_hline(yintercept = -2/sqrt(198), lty = 2, colour = "blue4") + 
        labs(x = "Lag", y = "Correlation", title = "A") +
        theme(axis.text = element_text(size = 6),
                axis.title = element_text(size = 6),
                panel.background = element_rect(fill = "white"),
                panel.grid = element_line(colour = "grey85"), 
              panel.border = element_rect(colour = "black", fill=NA, size=1)),
        
        ggplot(cors %>% filter(name=="B"), aes(x = lag, y = cor)) +
                geom_segment(aes(x = lag, y = 0, xend = lag, yend = cor)) + 
                geom_hline(yintercept = 2/sqrt(198), lty = 2, colour = "blue4") + 
                geom_hline(yintercept = -2/sqrt(198), lty = 2, colour = "blue4") + 
                labs(x = "Lag", y = "Correlation", title = "B") +
                theme(axis.text = element_text(size = 6),
                      axis.title = element_text(size = 6),
                      panel.background = element_rect(fill = "white"),
                      panel.grid = element_line(colour = "grey85"), 
                      panel.border = element_rect(colour = "black", fill=NA, size=1)), 
        
        ggplot(cors %>% filter(name=="C"), aes(x = lag, y = cor)) +
                geom_segment(aes(x = lag, y = 0, xend = lag, yend = cor)) + 
                geom_hline(yintercept = 2/sqrt(198), lty = 2, colour = "blue4") + 
                geom_hline(yintercept = -2/sqrt(198), lty = 2, colour = "blue4") + 
                labs(x = "Lag", y = "Correlation", title = "C") +
                theme(axis.text = element_text(size = 6),
                      axis.title = element_text(size = 6),
                      panel.background = element_rect(fill = "white"),
                      panel.grid = element_line(colour = "grey85"), 
                      panel.border = element_rect(colour = "black", fill=NA, size=1)), 
        
        ggplot(cors %>% filter(name=="D"), aes(x = lag, y = cor)) +
                geom_segment(aes(x = lag, y = 0, xend = lag, yend = cor)) + 
                geom_hline(yintercept = 2/sqrt(198), lty = 2, colour = "blue4") + 
                geom_hline(yintercept = -2/sqrt(198), lty = 2, colour = "blue4") + 
                labs(x = "Lag", y = "Correlation", title="D") +
                theme(axis.text = element_text(size = 6),
                      axis.title = element_text(size = 6),
                      panel.background = element_rect(fill = "white"),
                      panel.grid = element_line(colour = "grey85"), 
                      panel.border = element_rect(colour = "black", fill=NA, size=1)),
        nrow=2)
fig3

## Normality test
normality.test(stb, multivariate.only = FALSE)

qqnorm(resid(vecm)[,1])
qqline(resid(vecm)[,1])
shapiro.test(resid(vecm)[,1])
qqnorm(resid(vecm)[,2])
qqline(resid(vecm)[,2])

## Table 1
summary(cajorls(H1, r = 1)$rlm)
tab11 <- summary(cajorls(H1, r = 1)$rlm)[["Response Census.d"]]$coefficients[c(2:7, 1, 8:19), c(1, 3, 4)] %>%
        cbind(
                summary(cajorls(H1, r = 1)$rlm)[["Response Incidence.d"]]$coefficients[c(2:7, 1, 8:19), c(1, 3, 4)]  
        )
colnames(tab11) <- c("Estimate1", "T-statistics1", "P-value1",
                     "Estimate2", "T-statistics2", "P-value2")

as_tibble(as.data.frame(tab11)) %>% mutate(Predictor = rownames(tab11), .before = 1) %>% 
        xtable(align = c("lccccccc"), digits = 4, caption = "Table 1. Parameter estimates and T-tests.") %>%
        print(include.rownames = F)

## MAPE (Dec 06 - Dec 12)
first_try <- predict(stb, n.ahead = 7, ci = 0.8)
Obs <- c(336, 347, 353, 348, 345, 362, 344)
Fcst <- first_try$fcst$Census[,1] %>% back_tran()
100*mean(abs(Obs - Fcst)/Obs, na.rm = TRUE)

## Bootstrapped forecast interval for Dec 06 - Dec 12
seasPred <- seasonaldummy(x1, 7)
seasPred[seasPred == 0] <- -1/7 ## Center columns to have mean 0 within week
seasPred[seasPred == 1] <- 6/7

set.seed(123)
e <- residuals(stb)
pred_path <- matrix(NA, nrow = 1000, ncol = 7)

for(rep in 1:1000){
        pred_data <- dat[199:205, ]
        for(i in 1:7){
                e_temp <- e[sample(1:198, size = 1, replace = TRUE), ]
                pred_temp <- predict(vecm, n.ahead = 1, newdata = pred_data[i:(i + 6), ], exoPred = matrix(seasPred[i, ], nrow = 1, ncol = 6))
                y_temp <- pred_temp + e_temp
                pred_data <- rbind(pred_data, y_temp)
        } 
        pred_path[rep, ] <- pred_data[8:14, 1]
}
Fcst_bstr <- pred_path %>% apply(2, quantile, p = c(0.1, 0.9)) %>% t()

## Fig 6. 
upper <- c(rep(NA, 7), fitted(vecm, level = "original")[,1] + 2*sqrt(summary(vecm)$sigma[1,1]))
upper <- back_tran(upper)
lower <- c(rep(NA, 7), fitted(vecm, level = "original")[,1] - 2*sqrt(summary(vecm)$sigma[1,1]))
lower <- back_tran(lower)
my_dates <- seq(as.Date("2020/05/15"), as.Date("2020/12/05"), "days")
poly01 <- tibble(Date = c(my_dates, rev(my_dates)), Bounds = c(upper, rev(lower)))
my_dates02 <- seq(as.Date("2020/12/06"), as.Date("2020/12/12"), "days")
poly02 <- tibble(Date = c(my_dates02, rev(my_dates02)), Bounds = c(back_tran(Fcst_bstr[, 1]), rev(back_tran(Fcst_bstr[, 2]))))
fig6_graph <- tibble(Census = c(master_data$Census[1:205], rep(NA, 7)),  
                    Date = c(my_dates, my_dates02),
                    Fit = c(rep(NA, 7), back_tran(fitted(vecm, level = "original")[,1]), 
                            back_tran(first_try$fcst$Census[,1])),
                    Obs = c(rep(NA, 205), Obs))

fig6 <- ggplot(fig6_graph, aes(x = Date))
fig6 <- fig6 + geom_polygon(data = poly01, aes(x = Date, y = Bounds), fill = "#00BFC4", alpha = 0.4)
fig6 <- fig6 + geom_polygon(data = poly02, aes(x = Date, y = Bounds), fill = "#F8766D", alpha = 0.4)
fig6 <- fig6 + geom_line(aes(y = Fit), colour = "red")
fig6 <- fig6 + geom_point(aes(y = Census), size = 0.5, color = "black")
fig6 <- fig6 + geom_point(aes(y = Obs), size = 0.5, color = "black")
fig6 <- fig6 + scale_x_date(date_breaks = "1 month", date_labels = "%b",  
                          limits = as.Date(c("2020-05-15", "2020-12-12")))  
fig6 <- fig6 + labs(x = "Date", y = "COVID-19 Hospital Census")
fig6 <- fig6 + theme(axis.text = element_text(size = 6),
                     axis.title = element_text(size = 6),
                     panel.background = element_rect(fill = "white"),
                     panel.grid = element_line(colour = "grey85"),
                     panel.border = element_rect(colour = "black", fill=NA, size=1))
fig6


## ---------------------------------------------------------------------------- ##

## Cross-validation
set.seed(123)
place_holder <- data.frame(matrix(NA, nrow = 166, ncol = 3))
colnames(place_holder) <- c("Time", "Eigen", "MAPE")
for(i in 33:198){
        place_holder[i - 32, "Time"] <- i
        dat_cv <- master_data[1:i, c("Lcensus", "Lincidence")]
        colnames(dat_cv) <- c("Census", "Incidence")
        H1_cv <- ca.jo(dat_cv, type = "trace", ecdet = "const", spec = "transitory", K = 7, season = 7)
        ## Max eigenvalue modulus
        stb_cv <- vec2var(H1_cv, r = 1)
        eigen_cv <- Mod(eigen(companion(stb_cv))$values)
        place_holder[i - 32, "Eigen"] <- ifelse(round(max(eigen_cv), 7) == 1, eigen_cv[2], eigen_cv[1])
        ## MAPE 
        preds_cv <- predict(stb_cv, n.ahead = 7, ci = 0.95)
        preds_cv <- preds_cv$fcst$Census[, "fcst"] %>% back_tran()
        Obs_cv <- unlist(master_data[(i + 1):(i+7), "Census"])
        place_holder[i - 32, "MAPE"] <- mean(abs((preds_cv - Obs_cv)/Obs_cv))
}
place_holder <- tibble(place_holder)  %>% mutate(Date = seq(as.Date("2020/06/16"), as.Date("2020/11/28"), "days"))

## Fig 4.
fig4 <- ggplot(place_holder) + 
        geom_line(aes(x = Date, y = Eigen)) + 
        labs(x = "Date", y = "Maximum Eigenvalue Modulus") + 
        scale_x_date(breaks = "1 month", date_labels = "%b") + 
        theme(axis.text = element_text(size = 6),
              axis.title = element_text(size = 6),
              panel.background = element_rect(fill = "white"),
              panel.grid = element_line(colour = "grey85"),
              panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
        geom_hline(yintercept = 1, linetype = 2)
        
fig4

## MAPE 50th, 95th percentile
quantile(place_holder$MAPE, c(0.5, 0.95))

## Fig 5.
fig5 <- ggplot(place_holder)+
        geom_histogram(aes(x = place_holder$MAPE), binwidth = 0.01, color = "black", fill = "white") + 
        xlim(0, 0.3) + 
        labs(x = "Mean Absolute Percentage Error", y = "Count") + 
        geom_vline(xintercept = quantile(place_holder$MAPE, 0.5), linetype = 2, color = "blue", size = 1.2) + 
        geom_vline(xintercept = quantile(place_holder$MAPE, 0.95), linetype = 2, color = "red", size = 1.2) + 
        theme(axis.text = element_text(size = 6),
              axis.title = element_text(size = 6),
              panel.background = element_rect(fill = "white"),
              panel.grid = element_line(colour = "grey85"),
              panel.border = element_rect(colour = "black", fill=NA, size=1))
fig5
