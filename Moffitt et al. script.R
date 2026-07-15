rm(list=ls())
library(MASS); library(car); library(dplyr); library(MuMIn); library(ggplot2); library(tidyr)

##############################################################################################################
########################################### Snail mortality ##################################################
##############################################################################################################

snail_dat<-read.csv(choose.files(), header =T, sep =",")
snail$Group<- factor(snail$Group, levels = c("A","J"))
snail$Season<- factor(snail$Season, levels = c("S","W"))

#Model selection
mod1 <- glm(cbind(Dead,Alive) ~ Temp,family = binomial(link = "logit"), data = snail)
mod2 <- glm(cbind(Dead,Alive) ~ Temp+Group+Season, family = binomial(link = "logit"), data = snail)
mod3 <- glm(cbind(Dead,Alive) ~ Temp*Group*Season, family = binomial(link = "logit"), data = snail)
mod4 <- glm(cbind(Dead,Alive) ~ Season:Group+ Temp:Season:Group -1, family = binomial(link = "logit"), data = snail)
AICc(mod1,mod2,mod3, mod4)

options(contrasts = c("contr.sum", "contr.poly"))
Anova(mod5, type="3", test.statistic = "F") # significant interaction of slope between seasons and ages
summary(mod4)

######### obtain the LT50 and LT90 with bootstrapping 95% CI ############
get_LT <- function(model, data, B = 10000){
  coef_hat <- coef(model)
  coef_sim <- MASS::mvrnorm(n = B, mu = coef_hat, Sigma = vcov(model))
  groups <- expand.grid(Season = sort(unique(data$Season)),Group = sort(unique(data$Group)))
  LT_list <- list()
  
  for(i in 1:nrow(groups)){
    season <- groups$Season[i]
    group <- groups$Group[i]
    int_name <-paste0("Season", season, ":Group", group)
    slope_name <- paste0("Season", season, ":Group", group, ":Temp")
    alpha <- coef_hat[int_name]
    beta  <- coef_hat[slope_name]
    alpha_sim <- coef_sim[, int_name]
    beta_sim  <- coef_sim[, slope_name]
    LT50 <- (qlogis(0.5) - alpha) / beta
    LT90 <- (qlogis(0.9) - alpha) / beta
    LT50_sim <- (qlogis(0.5) - alpha_sim) / beta_sim
    LT90_sim <- (qlogis(0.9) - alpha_sim) / beta_sim
    LT_list[[i]] <- data.frame(Season = season, Group = group,
                               LT50 = LT50, LT50_L = quantile(LT50_sim, 0.025),  
                               LT50_U = quantile(LT50_sim, 0.975), LT90 = LT90,
                               LT90_L = quantile(LT90_sim, 0.025),
                               LT90_U = quantile(LT90_sim, 0.975),
                               LT50_sim = I(list(LT50_sim)),
                               LT90_sim = I(list(LT90_sim)))
  }
  
  bind_rows(LT_list)
}

LT_results <- get_LT(model = mod5, data = snail, B = 10000)
LT_results


####### compare the LT50 and LT90 by season and age ###########
compare_LT_nested <- function(LT_results, metric = "LT50"){
  sims <- LT_results[[paste0(metric,"_sim")]]
  LT_results$id <- seq_len(nrow(LT_results))
  
  # Compare Groups within Season
  season_results <- list()
  seasons <- unique(LT_results$Season)
  k <- 1
  
  for(s in seasons){
    sub <- LT_results[LT_results$Season == s, ]
    if(nrow(sub) != 2) next
    diff <- sub[[paste0(metric,"_sim")]][[1]] -
      sub[[paste0(metric,"_sim")]][[2]]
    p <- 2 * min(mean(diff > 0), mean(diff < 0))
    
    season_results[[k]] <- data.frame(Metric = metric,Season = s,
                                      Group1 = sub$Group[1], 
                                      Group2 = sub$Group[2],
                                      Mean_Diff = mean(diff),
                                      CI_Low = quantile(diff,0.025),
                                      CI_High = quantile(diff,0.975),
                                      p_value = p)
    
    k <- k + 1
  }
  
  season_results <- bind_rows(season_results)
  
  # Compare Seasons within Group
  group_results <- list()
  groups <- unique(LT_results$Group)
  k <- 1
  
  for(g in groups){
    
    sub <- LT_results[LT_results$Group == g, ]
    if(nrow(sub) != 2) next
    diff <- sub[[paste0(metric,"_sim")]][[1]] -
      sub[[paste0(metric,"_sim")]][[2]]
    p <- 2 * min(mean(diff > 0), mean(diff < 0))
    
    group_results[[k]] <- data.frame(Metric = metric, Group = g,
                                     Season1 = sub$Season[1],
                                     Season2 = sub$Season[2],
                                     Mean_Diff = mean(diff),
                                     CI_Low = quantile(diff,0.025),
                                     CI_High = quantile(diff,0.975),
                                     p_value = p)
    
    k <- k + 1
  }
  
  group_results <- bind_rows(group_results)
  
  list(season_comparisons = season_results,
       group_comparisons = group_results)
}

LT50_nested <- compare_LT_nested(LT_results, metric = "LT50")
LT90_nested <- compare_LT_nested(LT_results, metric = "LT90")

LT50_nested
LT90_nested


##### create fitted values with 95% CI from bootstrapping #####
B <- 10000
coef_sim <- mvrnorm(n = B, mu = coef(mod5), Sigma = vcov(mod5))

groups <- data.frame(Season = c("S","S","W","W"), Group  = c("A","J","A","J"))

temp_grid <- seq(20,60,by=0.1)
pred_grid <- expand.grid(Temp = temp_grid, Season = c("S","W"), Group = c("A","J"))

X <- model.matrix(~ Season:Group -1 + Temp:Season:Group, data = pred_grid)
eta_hat <- X %*% coef(mod5)
pred_grid$Fit <- plogis(eta_hat)
pred_sim <- X %*% t(coef_sim)
prob_sim <- plogis(pred_sim)

pred_grid$Lower <-  apply(prob_sim,1,quantile,0.025)
pred_grid$Upper <-  apply(prob_sim,1,quantile,0.975)

fitted_df <- pred_grid %>%
  select(Temp, Season, Group, Fit, Lower, Upper)


#Plot the mortality curve with fitted and real data
snail$Prop <- snail$Dead / (snail$Dead + snail$Alive)*100

ggplot(pred_grid, aes(x = Temp, y = Fit*100,  color = Season, fill = Season)) +
  geom_ribbon(aes(ymin = Lower*100, ymax = Upper*100), alpha = 0.2, colour = NA) +
  geom_point(data = snail, aes(x = Temp, y = Prop, color = Season), alpha = 0.6, 
             position = position_jitter(width = 0.2, height = 0)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c(S = "red", W = "steelblue"), labels = c("Summer", "Winter")) +
  scale_fill_manual(values = c(S = "red", W = "steelblue"), labels = c("Summer", "Winter"))+
  facet_grid(rows = vars(Group), scales="free_y", labeller = labeller(Group = c( "A" = "Adult", "J" = "Juvenile"))) +
  labs(x = "Temperature (°C)", y = "Mortality (%)")+ geom_line(size = 1) +
  geom_vline(data = LT_results, aes(xintercept = LT50,
                                    color = Season), linetype = "solid", linewidth = 0.8, show.legend = FALSE) +
  geom_vline(data = LT_results, aes(xintercept = LT90,
                                    color = Season), linetype = "dashed", linewidth = 0.8, show.legend = FALSE) +
  theme_test()+
  theme(text = element_text(color = "black", size = 20), axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"), strip.text = element_text(color = "black"),
        legend.text = element_text(color = "black"), legend.title = element_text(color = "black"))


#Plot the LT50 and LT90 with asymmetrical 95%CI
#create dataframe for plotting
LT_plot <- LT_results %>%
  select(Season, Group, LT50, LT50_L, LT50_U, LT90, LT90_L, LT90_U) %>%
  pivot_longer(cols = c(LT50, LT90), names_to = "Metric", values_to = "Estimate") %>%
  mutate(Lower = ifelse(Metric=="LT50", LT50_L, LT90_L),
         Upper = ifelse(Metric=="LT50", LT50_U, LT90_U))

LT_plot$Season <- factor(LT_plot$Season, levels = c("S", "W"))
LT_plot$Group <- factor(LT_plot$Group, levels = c("A", "J"), labels = c("Adult", "Juvenile"))
LT_plot$int<-interaction(LT_plot$Season,LT_plot$Group)
LT_plot$int<- factor(LT_plot$int, levels = c("S.A", "S.J", "W.A", "W.J"))

#plotting the results
ggplot(LT_plot,  aes(x = Group, y = Estimate, fill = Season)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.59) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), position = position_dodge(width = 0.6),
                width = 0.15, linewidth = 0.8) +
  facet_wrap(~Metric, labeller = as_labeller(c(LT50 = "LT[50]",
                                               LT90 = "LT[90]"), default = label_parsed)) +
  scale_fill_manual(values = c(S = "#eb1427", W = "steelblue"), labels = c(S = "Summer", W = "Winter")) +
  coord_cartesian(ylim = c(30, 50)) +
  theme_test() +
  labs(x = NULL, y = "Temperature (°C)")+
  theme(text = element_text(color = "black", size = 20), axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"), strip.text = element_text(color = "black"),
        legend.text = element_text(color = "black"), legend.title = element_text(color = "black"))


##############################################################################################################
############################################### oyster mortality #############################################
##############################################################################################################
oyster<-read.csv(choose.files(), header =T, sep =",")
oyster$Season<- factor(oyster$Season, levels = c("S","W"))

#Model selection
mod1 <- glm(cbind(Dead,Alive) ~ Temp, family = binomial(link = "logit"), data = oyster)
mod2 <- glm(cbind(Dead,Alive) ~ Temp + Season-1, family = binomial(link = "logit"), data = oyster)
mod3 <- glm(cbind(Dead,Alive) ~  Temp * Season -1, family = binomial(link = "logit"), data = oyster)
AICc(mod1,mod2,mod3)
anova(mod1)

######### obtain the LT50 and LT90 with  95% CI using bootstrapping ############
get_LT <- function(model, data, B = 10000){
  
  coef_hat <- coef(model)
  
  coef_sim <- MASS::mvrnorm(n = B, mu = coef_hat, Sigma = vcov(model))
  alpha <- coef_hat["(Intercept)"]
  beta  <- coef_hat["Temp"]
  alpha_sim <- coef_sim[, "(Intercept)"]
  beta_sim  <- coef_sim[, "Temp"]
  LT50 <- (qlogis(0.5) - alpha) / beta
  LT90 <- (qlogis(0.9) - alpha) / beta
  LT50_sim <- (qlogis(0.5) - alpha_sim) / beta_sim
  LT90_sim <- (qlogis(0.9) - alpha_sim) / beta_sim
  
  data.frame(LT50 = LT50,LT50_L = quantile(LT50_sim, 0.025),LT50_U = quantile(LT50_sim, 0.975),
             LT90 = LT90, LT90_L = quantile(LT90_sim, 0.025), LT90_U = quantile(LT90_sim, 0.975),
             LT50_sim = I(list(LT50_sim)), LT90_sim = I(list(LT90_sim)))
}

LT_results <- get_LT(mod1, oyster, B = 10000)
LT_results

##### create fitted values with 95% CI from bootstrapping #####

B <- 10000
coef_sim <- MASS::mvrnorm(n = B, mu = coef(mod1), Sigma = vcov(mod1))
temp_grid <- seq(20, 64, by = 0.1)
pred_grid <- expand.grid(Temp = temp_grid)

X <- model.matrix(~ Temp, data = pred_grid)
eta_hat <- X %*% coef(mod1)
pred_grid$Fit <- plogis(eta_hat)

pred_sim <- plogis(X %*% t(coef_sim))
pred_grid$Lower <- apply(pred_sim, 1, quantile, 0.025)
pred_grid$Upper <- apply(pred_sim, 1, quantile, 0.975)

fitted_df <- pred_grid

#### plotting the fitted and real data #####
oyster$Prop <- oyster$Dead / (oyster$Dead + oyster$Alive) * 100

ggplot() +  geom_ribbon(data = fitted_df, aes(x = Temp, ymin = Lower*100, ymax = Upper*100),
                        alpha = 0.2, colour = "grey") +
  geom_line(data = fitted_df, aes(x = Temp, y = Fit*100), linewidth = 1) +
  geom_point(data = oyster, aes(x = Temp, y = Prop, colour = Season),
             alpha = 0.6, position = position_jitter(width = 0.2), size=3) +
  scale_colour_manual(values = c(S = "#eb1427", W = "steelblue"), labels = c(S = "Summer", W = "Winter")) +
  geom_vline(data = LT_results, aes(xintercept = LT50), linetype = "solid", linewidth = 0.8) +
  geom_vline(data = LT_results, aes(xintercept = LT90), linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = LT_results$LT50+3, y = 23, label = "49.5 °C", angle = 90, vjust = -0.3, size=5) +
  annotate("text", x = LT_results$LT90+3, y = 23, label = "58.5 °C", angle = 90, vjust = -0.3, size=5) +
  labs(x = "Temperature (°C)", y = "Mortality (%)") +
  theme_test(base_size = 16)+
  theme(text = element_text(color = "black", size = 20), axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"), strip.text = element_text(color = "black"),
        legend.text = element_text(color = "black"), legend.title = element_text(color = "black"))


######################################################################################################
###################################### Barnacle mortality#############################################
######################################################################################################
barnacle<-read.csv(choose.files(), header =T, sep =",")
barnacle$Season<- factor(barnacle$Season, levels = c("S","W"))

#Model selection
mod1 <- glm(cbind(Dead_count,Alive_count) ~  SA + Temp, family = binomial(link = "logit"), data = barnacle)
mod2 <- glm(cbind(Dead_count,Alive_count) ~ SA + Temp + Season -1, family = binomial(link = "logit"), data = barnacle)
mod3 <- glm(cbind(Dead_count,Alive_count) ~  Season + Temp:Season + SA-1, family = binomial(link = "logit"), data = barnacle)
AICc(mod1,mod2,mod3)

options(contrasts = c("contr.sum", "contr.poly"))
Anova(mod3, test.statistic = "F", type="3") # significant interaction of slope between seasons and ages
summary(mod3) # although mod1 lowest AICc, I want to test if there are seasonal difference, use mod2 share the same slope


######### obtain the LT50 and LT90 with  95% CI using bootstrapping ############
get_LT <- function(model, data, SA_ref = mean(data$SA, na.rm = TRUE), B = 10000){
  coef_hat <- coef(model)
  coef_sim <- MASS::mvrnorm(n = B, mu = coef_hat, Sigma = vcov(model))
  seasons <- sort(unique(data$Season))
  LT_list <- vector("list", length(seasons))
  
  for(i in seq_along(seasons)){
    s <- seasons[i]
    alpha <- coef_hat[paste0("Season", s)] + coef_hat["SA"] * SA_ref
    alpha_sim <- coef_sim[, paste0("Season", s)] + coef_sim[, "SA"] * SA_ref
    beta_name <- paste0("Season", s, ":Temp")
    beta <- coef_hat[beta_name]
    beta_sim <- coef_sim[, beta_name]
    ## LT50 and LT90
    LT50 <- (qlogis(0.5) - alpha) / beta
    LT90 <- (qlogis(0.9) - alpha) / beta
    LT50_sim <- (qlogis(0.5) - alpha_sim) / beta_sim
    LT90_sim <- (qlogis(0.9) - alpha_sim) / beta_sim
    LT_list[[i]] <- data.frame(Season = s, SA = SA_ref,LT50 = LT50,
                               LT50_L = quantile(LT50_sim, 0.025),
                               LT50_U = quantile(LT50_sim, 0.975),
                               LT90 = LT90, LT90_L = quantile(LT90_sim, 0.025),
                               LT90_U = quantile(LT90_sim, 0.975),
                               LT50_sim = I(list(LT50_sim)),
                               LT90_sim = I(list(LT90_sim)))
  }
  
  bind_rows(LT_list)
}

LT_results <- get_LT(model = mod3, data = barnacle,  SA_ref = mean(barnacle$SA), B = 10000)
LT_results


compare_LT_season <- function(LT_results, metric = "LT50") {
  

  sims <- LT_results[[paste0(metric, "_sim")]]
  diff <- sims[[1]] - sims[[2]]
  B <- length(diff)
  p <- 2 * min((sum(diff >= 0) + 1)/(B + 1),
               (sum(diff <= 0) + 1)/(B + 1))
  prob_gt <- mean(diff > 0)
  
  data.frame(Metric = metric, Season1 = LT_results$Season[1], Season2   = LT_results$Season[2],
             Mean_Diff = mean(diff), CI_Low    = quantile(diff, 0.025),
             CI_High   = quantile(diff, 0.975), Prob_GT   = prob_gt,
             p_value   = format.pval(p, digits = 3,eps = .Machine$double.xmin))
}

LT50_compare <- compare_LT_season(LT_results, "LT50")
LT90_compare <- compare_LT_season(LT_results, "LT90")

LT50_compare
LT90_compare


##### create fitted values with 95% CI from bootstrapping #####
B <- 10000 # 10000 time simulations

coef_sim <- MASS::mvrnorm(n = B, mu = coef(mod3), Sigma = vcov(mod3))
temp_grid <- seq(20,60,by=0.1)
SA_ref <- mean(barnacle$SA)
pred_grid <- expand.grid(Temp = temp_grid, Season = c("S","W"))
pred_grid$SA <- SA_ref

X <- model.matrix(~ Season + Temp:Season + SA -1, data = pred_grid)
eta_hat <- X %*% coef(mod3)
pred_grid$Fit <- plogis(eta_hat)
pred_sim <- plogis( X %*% t(coef_sim))
pred_grid$Lower <-  apply(pred_sim,1,quantile,0.025)
pred_grid$Upper <-  apply(pred_sim,1,quantile,0.975)

fitted_df <- pred_grid


#Plot the mortality curve with fitted and real data
barnacle$Prop <-  barnacle$Dead_count / (barnacle$Dead_count + barnacle$Alive_count) * 100

ggplot() +  geom_ribbon(data = fitted_df,aes(x = Temp, ymin = Lower*100, ymax = Upper*100, fill = Season),
                        alpha = 0.2, colour = NA) +
  geom_line(data = fitted_df, aes(x = Temp, y = Fit*100,colour = Season), linewidth = 1) +
  geom_point(data = barnacle, aes(x = Temp, y = Prop, colour = Season), alpha = 0.6,
             position = position_jitter(width = 0.2)) +
  geom_vline(data = LT_results, aes(xintercept = LT50, colour = Season), linetype = "solid") +
  geom_vline(data = LT_results, aes(xintercept = LT90, colour = Season), linetype = "dashed") +
  scale_colour_manual(values = c(S = "#eb1427", W = "steelblue"), 
                      labels = c(S = "Summer", W = "Winter")) +
  scale_fill_manual(values = c(S = "#eb1427", W = "steelblue"),
                    labels = c(S = "Summer", W = "Winter")) +
  labs(x = "Temperature (°C)", y = "Mortality (%)") +
  theme_test(base_size = 18)+
  theme(text = element_text(color = "black", size = 20), axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"), strip.text = element_text(color = "black"),
        legend.text = element_text(color = "black"), legend.title = element_text(color = "black"))


#Plot the LT50 and LT90 with asymmetrical 95%CI
LT_plot <- LT_results %>%
  select(Season, LT50, LT50_L, LT50_U, LT90, LT90_L, LT90_U) %>%
  pivot_longer(cols = c(LT50, LT90),
               names_to = "Metric",
               values_to = "Estimate") %>%
  mutate(Lower = ifelse(Metric == "LT50", LT50_L, LT90_L),
         Upper = ifelse(Metric == "LT50", LT50_U, LT90_U))

LT_plot$Season <- factor(LT_plot$Season, levels = c("S","W"), labels = c("Summer","Winter"))

ggplot(LT_plot, aes(x = Season, y = Estimate, fill = Season)) +
  geom_col(width = 0.6) + geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.15) +
  facet_wrap(~Metric, labeller = as_labeller(c(LT50 = "LT[50]", LT90 = "LT[90]"),
                                             default = label_parsed)) +
  scale_fill_manual(values = c(Summer = "#eb1427", Winter = "steelblue")) +
  coord_cartesian(ylim = c(30,65)) +
  labs(x = NULL, y = "Temperature (°C)") +
  theme_test(base_size = 18)



##############################################################################################################
######################################### snail temporal mortality ###########################################
##############################################################################################################
temp_m <- read.csv(choose.files(), header=T, sep=",")
temp_m$Season <- factor(temp_m$Season, levels = c("S", "W"), labels = c(S = "Summer", W = "Winter"))
temp_m$Group <- factor(temp_m$Group, levels = c("A", "J"), labels = c(A = "Adult", J = "Juvenile"))

#Model selection
mod1 <- glm(cbind(Dead,Alive) ~ Time_min, family = binomial(link = "logit"), data = temp_m)
mod2 <- glm(cbind(Dead,Alive) ~ Time_min + Season -1, family = binomial(link = "logit"), data = temp_m)
mod3 <- glm(cbind(Dead,Alive) ~ Time_min * Season - Time_min-1, family = binomial(link = "logit"), data = temp_m)
mod4 <- glm(cbind(Dead,Alive) ~ Time_min * Season -1, family = binomial(link = "logit"), data = temp_m)
mod5 <- glm(cbind(Dead,Alive) ~ Time_min * Group * Season -1, family = binomial(link = "logit"), data = temp_m)
AICc(mod1,mod2, mod3,mod4, mod5)

options(contrasts = c("contr.sum", "contr.poly"))
Anova(mod3, test.statistic = "F", type="3") # significant interaction of slope between seasons and ages


######### obtain the LTime50 and LTime90 with  95% CI using bootstrapping ############
get_LTime <- function(model, data, B = 10000){
  
  coef_hat <- coef(model)
  coef_sim <- MASS::mvrnorm(B, coef_hat, vcov(model))
  seasons <- levels(data$Season)
  out <- vector("list", length(seasons))
  
  for(i in seq_along(seasons)){
    s <- seasons[i]
    alpha_name <- paste0("Season", s)
    beta_name  <- paste0("Time_min:Season", s)
    alpha <- coef_hat[alpha_name]
    beta  <- coef_hat[beta_name]
      alpha_sim <- coef_sim[, alpha_name]
    beta_sim  <- coef_sim[, beta_name]
    LTime50 <- (qlogis(0.5)-alpha)/beta
    LTime90 <- (qlogis(0.9)-alpha)/beta
    LTime50_sim <- (qlogis(0.5)-alpha_sim)/beta_sim
    LTime90_sim <- (qlogis(0.9)-alpha_sim)/beta_sim
    
    out[[i]] <- data.frame(Season = s, LTime50 = LTime50,
                           LTime50_L = quantile(LTime50_sim,.025), LTime50_U = quantile(LTime50_sim,.975),
                           LTime90 = LTime90, LTime90_L = quantile(LTime90_sim,.025), 
                           LTime90_U = quantile(LTime90_sim,.975),
                           LTime50_sim = I(list(LTime50_sim)),
                           LTime90_sim = I(list(LTime90_sim)))
  }
  
  bind_rows(out)
}

LTime_results <- get_LTime(mod3, temp_m)
LTime_results

compare_LT_season <- function(LTime_results, metric = "LTime50"){
    sims <- LTime_results[[paste0(metric, "_sim")]]
    ord <- order(LTime_results$Season)
    sims <- sims[ord]
    labels <- LTime_results$Season[ord]
    diff <- sims[[1]] - sims[[2]]
    B <- length(diff)
    p <- 2 * min((sum(diff >= 0) + 1) / (B + 1), (sum(diff <= 0) + 1) / (B + 1))
    data.frame(Metric = metric, Season1 = labels[1], Season2 = labels[2],
               Mean_Diff = mean(diff), CI_Low = quantile(diff, 0.025),
               CI_High = quantile(diff, 0.975), Prob_GT = mean(diff > 0),
               p_value = format.pval(p, digits = 3))
}

LT50_posthoc <- compare_LT_season (LTime_results, "LTime50")
LT90_posthoc <- compare_LT_season (LTime_results, "LTime90")

LT50_posthoc
LT90_posthoc


##### create fitted values with 95% CI from bootstrapping #####
B <- 10000

coef_sim <- MASS::mvrnorm(B,coef(mod3),vcov(mod3))

time_grid <- seq(0, 63, by = 0.5)
pred_grid <- expand.grid(Time_min = time_grid, Season = c("Summer", "Winter"))

X <- model.matrix(~Time_min*Season -Time_min- 1, data = pred_grid)
eta <- X %*% coef(mod2.1)
pred_grid$Fit <- plogis(eta)
pred_sim <- plogis(X %*% t(coef_sim))

pred_grid$Lower <- apply(pred_sim, 1, quantile, 0.025)
pred_grid$Upper <- apply(pred_sim, 1, quantile, 0.975)


#Plot the mortality curve with fitted and real data
temp_m$Prop <- 100 * temp_m$Dead / (temp_m$Dead + temp_m$Alive)

ggplot() +  geom_ribbon(data = pred_grid, aes(Time_min, ymin = Lower * 100, ymax = Upper * 100, 
                                              fill = Season), alpha = 0.2) +
  geom_line(data = pred_grid, aes(Time_min, Fit * 100, colour = Season), linewidth = 1) +
  geom_point(data = temp_m, aes(Time_min, Prop, colour = Season), alpha = 0.6) +
  geom_vline(data = LTime_results, aes(xintercept = LTime50, colour = Season), linetype = "solid") +
  geom_vline(data = LTime_results, aes(xintercept = LTime90, colour = Season), linetype = "dashed") +
  scale_colour_manual(values = c(Summer = "#eb1427", Winter = "steelblue")) +
  scale_fill_manual(values = c(Summer = "#eb1427", Winter = "steelblue")) +
  labs(x = "Exposure time (min)", y = "Mortality (%)") +
  theme_test(base_size = 18) +
  theme(text = element_text(color = "black", size = 20), axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"), strip.text = element_text(color = "black"),
        legend.text = element_text(color = "black"), legend.title = element_text(color = "black"))


#Plot barchart LTime50 and LTime90 with asymeetrical 95%CI
LTime_plot <- LTime_results %>%
  pivot_longer(cols = c(LTime50, LTime90), names_to = "Metric", values_to = "Estimate") %>%
  mutate(Lower = ifelse(Metric == "LTime50", LTime50_L, LTime90_L), 
         Upper = ifelse(Metric == "LTime50", LTime50_U, LTime90_U))
pd <- position_dodge(width = 0.6)

ggplot(LTime_plot,aes(x = Season, y = Estimate, fill = Season)) +
    geom_col(position = pd, width = 0.6) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper),position = pd,width = 0.15) +
    scale_fill_manual(values = c(Summer = "#eb1427",Winter = "steelblue")) +
    facet_wrap(~ Metric, labeller = labeller(Metric = c(LTime50 = "LTime[50]",
                                                        LTime90 = "LTime[90]" ), 
                                             .default = label_parsed)) +
    labs(x = NULL, y = "Exposure time (min)") +
    theme_test(base_size = 18) +
    theme(text = element_text(color = "black", size = 20), axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"), strip.text = element_text(color = "black"),
        legend.text = element_text(color = "black"), legend.title = element_text(color = "black"))


##############################################################################################################
######################################### snail DDT ##########################################################
##############################################################################################################
ddt <- read.csv(choose.files(), header=T, sep=",")
ddt<-subset(ddt, !Group =="S")
ddt$Season <- factor(ddt$Season, levels = c("S", "W"), labels = c(S = "Summer", W = "Winter"))
ddt$Group <- factor(ddt$Group, levels = c("A", "J"), labels = c(A = "Adult", J = "Juvenile"))

shapiro.test(ddt$DDT)
hist(ddt$DDT) # right-skewed, use Gamma distribution

mods1<- glm(DDT ~ Season * Group, family = Gamma(link = "log"), data = ddt)

options("contrasts")
contrasts(ddt$Season)
contrasts(ddt$Group)

options(contrasts = c("contr.sum", "contr.poly"))
Anova(mods2, type="III", test.statistic = "F")
emmeans::emmeans(mods1, pairwise ~ Group|Season)

ggplot(ddt, aes(x = Season, y=DDT)) +
  geom_boxplot(aes(x = Season, y = DDT, fill = Season),
               staplewidth = 0.5, width = 0.7, outlier.shape = NA,median.linewidth = 1.2) +
  geom_jitter(data = ddt, aes(y= DDT, x = Season, fill = Season), shape=21, alpha=0.8,
              position = position_jitter(width = 0.1, height = 0), stroke=1.2)+
  facet_wrap(~Group) + labs(x= "Season", y="DDT (min)") +
  scale_fill_manual(values = c(Summer = "#eb1427", Winter = "steelblue")) +
  stat_summary( fun = mean, geom = "point", shape = 21, size = 4, fill = "yellow", stroke=1.2)+
  theme_test(base_size = 18) + theme(legend.position = "none")+
  theme(text = element_text(color = "black", size = 20), axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"), strip.text = element_text(color = "black"),
        legend.text = element_text(color = "black"), legend.title = element_text(color = "black"))

mean(subset(ddt, Group =="Adult")$DDT)
sd(subset(ddt, Group =="Adult")$DDT)

mean(subset(ddt, Group =="Juvenile")$DDT)
sd(subset(ddt, Group =="Juvenile")$DDT)
