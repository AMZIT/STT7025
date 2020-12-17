#-----------------------------------------------
# TRAVAIL PRATIQUE 2
# 
# Cours: STT-7120
# Professeur: M. Thierry Duchesne
# 
# Équipe 7
# Étudiants: - Alexandre Lepage 111 144 776
#            - Amedeo Zito
#-----------------------------------------------

library(tidyverse)
library(lme4)
library(nlme)
library(olsrr)

current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))

#============================= Question 1 ======================================
rm(list=ls())
data <- read.table(file='Homework.txt', header = T)
data %>%  glimpse()
data <- 
   data %>%
   select(-c(stuid, public, sex)) %>%
   mutate(white = factor(white),
          schid = factor(schid),
          ratio = (ratio-median(ratio)) / sd(ratio)
          )
data %>%  glimpse()
attach(data)


# a) ---------------------------------------------------------------------------
#--- Entraînement d'un modèle linéaire standard ----
lm_complet <- lm(
   math ~ (homework + white + ratio),
   data=data)
ols_vif_tol(lm_complet)
ols_plot_resid_fit(lm_complet)
#' Le modèle est biaisé. Les résidus ne sont pas centrés à zéro.
#' Voyons voir si des transformations sur les xi pourrait améliorer le modèle.

library(gam)
gam_mod <- gam(math ~ s(homework) + white + s(ratio))
plot(gam_mod, se=TRUE)
#' Possiblement que l'ajout d'une variable ratio2=(ratio-18)^2 pourrait 
#' améliorer le modèle. En ce qui attrait à homework, celle-ci peut demeurer 
#' telle quelle.
summary(gam_mod)
# Les deux transformations sont significatives au seuil de 1%.

detach(data)
data <- data %>% mutate(ratio2 = ratio^2)
attach(data)

gam_mod <- gam(math ~ s(homework) + white + s(ratio) + s(ratio2))
par(mfrow=c(2,2))
plot(gam_mod, se=TRUE)
par(mfrow=c(1,1))
summary(gam_mod)
# Les relations sont maintenant toutes linéaires.

lm_complet <- lm(
   math ~ homework + white + ratio + ratio2,
   data=data
)
# Analyse de multicolinéarité
ols_vif_tol(lm_complet)

# Ajout d'interactions
add1(lm_complet, .~. +.^2, test = "F")
lm_complet <- update(lm_complet, .~.+ white:ratio)
add1(lm_complet, .~. +.^2, test = "F")
lm_complet <- update(lm_complet, .~.+ white:ratio2)
add1(lm_complet, .~. +.^2, test = "F")

ols_vif_tol(lm_complet)

ols_plot_resid_fit(lm_complet)
#' Ce n'est pas bon, mais c'est le mieux que je puisse faire avec les variables
#' disponibles.

# Analyse BoxCox.
MASS::boxcox(lm_complet)
#' Comme lambda = 1 appartient à l'intervalle de confiance, on n'apporte aucune
#' transformation à Y.

summary(lm_complet)

# Calcul du MSE en entraînement.
Y_hat <- predict(lm_complet)
MSE_lm <- mean((math - Y_hat)^2)


#' Voyons voir si un modèle mixte pourrait expliquer ce biais.

#--- Analyse de corrélation ---
data_residuals <- data %>% mutate(
   residus_lm = rstudent(lm_complet),
   id = row_number())


plot_residuals_variable <- function(variable){
   #' Fonction qui affiche le graphique des résidus en fonction d'une variable
   #' explicative
   #' @param variable La variable explicative que l'on désire regarder 
   data_residuals %>% 
      ggplot() +
      geom_point(aes(
         x = variable,
         y = residus_lm,
         colour = schid
      ), alpha = 0.6) +
      geom_hline(yintercept = 0) +
      xlab(deparse(substitute(variable))) +
      ylab("Résidus")
}
plot_residuals_id <- function(n_samples, RandomSeed=FALSE) {
   #' Fonction qui affiche le graphique des résidus en fonction des indices i.
   #' @param n_samples Nombre d'observations à afficher.
   #' @param RandomSeed Valeur de l'ancrage aléatoire. Si False, aucun ancrage.
   if (!isFALSE(RandomSeed)) set.seed(RandomSeed)
   n_row = nrow(data_residuals)
   data_residuals[sample(1:n_row, n_samples), ] %>% 
      ggplot() +
      geom_point(aes(
         x = id,
         y = residus_lm,
         colour = schid
      ), alpha = 0.6) +
      geom_hline(yintercept = 0) +
      xlab("i") +
      ylab("Résidus")
}


plot_residuals_id()
#' La distribution des résidus est grandement influencée par les grappes.

plot_residuals_variable(homework)
plot_residuals_variable(white)
plot_residuals_variable(ratio)
plot_residuals_variable(ratio2)
#' Il est difficile de tirer des conclusions de ces graphiques.
#' Les tests du ratio des vraisemblances permettra de mettre les choses
#' au clair sur les effets aléatoires.


#--- Entraînement des modèles linéaires mixtes (lmm) complets ----
# des matrices de variance appropriées

lm_complet

lmm_VC_UN <- lmer(  # Modèle mixte avec structure de variance VC, UN
   math ~ homework + white + ratio + ratio2 + white:ratio + 
      white:ratio2 +
      (homework | schid),
   data=data,
   REML = TRUE
   )
lmm_VC_UN1 <- lmer(  # Modèle mixte avec structure de variance VC, UN(1)
   math ~ homework + white + ratio + ratio2 + white:ratio + 
      white:ratio2 +
      (homework || schid),
   data=data,
   REML = TRUE
)
lmm_CS_UN <- lme(# Modèle mixte avec structure de variance CS, UN
   fixed = math ~ homework + white + ratio + ratio2 + white:ratio +
      white:ratio2,
   random = ~ homework | schid,
   data = data,
   method='REML',
   correlation = corCompSymm()
)


extractAIC_lme <- function(object) {
   #' Fonction qui prend un modèle de la classe lme et qui calcule l'AIC
   #' @param object Modèle entraîné de la classe lme.
   df = sum(object$dims$ncol, object$dims$qvec)
   aic = -2 * (object$logLik - df)
   return(c(df, aic))
}


# Comparaison des AIC
tbl_AIC <- rbind(
      "VC UN"=extractAIC(lmm_VC_UN),
      "VC UN(1)"=extractAIC(lmm_VC_UN1),
      "CS UN"=extractAIC_lme(lmm_CS_UN)
   )
colnames(tbl_AIC) <- list("dl", "AIC")
xtable::xtable(tbl_AIC)
tbl_AIC

lmm_complet <- lmm_VC_UN
summary(lmm_complet)
#' On choisit la structure nonstructurée (UN) pour la matrice de variance D
#' et la structure par composantes de variance (VC) pour la matrice V.
#' 
#' En ce qui attrait à la matrice V, la structure AR(1) n'est pas idéale si on
#' se base sur le fait que les j sont des étudiants d'une classe... donc ne
#' correspondent pas à des observations pouvant être ordonnancées 
#' chronologiquement ou d'un point de vue spatial.
#' 
#' Le fait de prendre une structure CS, quant à elle, reviendrait à dire qu'il
#' existerait une corrélation entre les résultats de différents individus d'une
#' même école. À ce stade, on pourrait considérer une corrélation entre les 
#' élèves d'une même classe, mais, autrement, on peut présumer l'indépendance
#' des étudiants jusqu'à preuve du contraire. Dans ce cas, la structure 
#' par composantes de variance (VC) serait justifiée.


#--- Test des effets aléatoires ----
lrtest_lmm <- function(model_H0, model_H1) {
   #' Fonction qui calcule le test du rapport des vraisemblances pour tester 
   #' la nécessité des effets aléatoires.
   #' @param model_H0 Le modèle sous l'hypothèse nulle (modèle simple)
   #' @param model_H1 Le modèle sous H1 (modèle complet)
   test <- lmtest::lrtest(model_H0, model_H1)
   df <- test$Df[2]
   chisq <- test$Chisq[2]
   pval <- 0.5 * (2 - pchisq(chisq, df - 1) - pchisq(chisq, df))
   print(paste("chisq:", round(chisq, 2), "Df:", df, "p-value:", round(pval,3)),
         quote=F)
}


lmm_H0 <- update(
   lmm_complet,.~.- (homework | schid) + (1 | schid), REML=T
   )
lrtest_lmm(lmm_H0, lmm_complet)
#' On rejette H0 et on conserve l'effet aléatoire (homework | schid)

lmm_Q1a <- lmm_complet
summary(lmm_Q1a)


#--- Sélection des effets fixes ----
car::Anova(lmm_Q1a, type=3)
lmm_Q1a <- update(lmm_Q1a, .~.- white:ratio2)
car::Anova(lmm_Q1a, type=3)
lmm_Q1a <- update(lmm_Q1a, .~.- white:ratio)
car::Anova(lmm_Q1a, type=3)
lmm_Q1a <- update(lmm_Q1a, .~.- ratio2)
car::Anova(lmm_Q1a, type=3)
lmm_Q1a <- update(lmm_Q1a, .~.- ratio)
car::Anova(lmm_Q1a, type=3)

summary(lmm_Q1a)

sum_lmm_Qst1a <- summary(lmm_Q1a)
beta <- sum_lmm_Qst1a$coefficients[,1]
sd_beta <- sum_lmm_Qst1a$coefficients[,2]
tbl_betas <- cbind(
   "Estimateurs"=beta,
   "Écarts-types"=sd_beta,
   beta + cbind(rep(-1, 3), rep(1, 3)) * qnorm(0.975) * sd_beta
   )
xtable::xtable(tbl_betas)

sum_lmm_Qst1a$varcor$schid %*% diag(nrow=2)
sum_lmm_Qst1a$sigma^2

# Calcul du MSE.
Y_hat <- predict(lmm_Q1a)
(MSE_lmm_Q1a <- mean((data$math - Y_hat)^2))
MSE_lm
#' Belle amélioration de l'erreur quadratique

data %>% mutate(  # Visualisation des résidus pour le modèle mixte entraîné
   predictions = Y_hat,
   residus = rstudent(lmm_Q1a)
   ) %>%
   ggplot() +
   geom_point(aes(
      x=predictions,
      y=residus
   ), color='blue') +
   geom_hline(yintercept = 0)
#' Le prolème a changé: Les résidus sont linéaires, mais on a un problème
#' d'hétéroscédasticité.
data %>% arrange(schid)%>%
   mutate(  # Visualisation des résidus pour le modèle mixte entraîné
   residus = rstudent(lmm_Q1a),
   id = row_number()
   ) %>%
   ggplot() +
   geom_point(aes(
      x=id,
      y=residus
   ), color=schid) +
   geom_hline(yintercept = 0)


# b) ---------------------------------------------------------------------------
#--- Entraînement d'un modèle linéaire standard ----
lm_complet <- lm(math ~ meanses + homework + white + ratio + ratio2, data=data)
ols_vif_tol(lm_complet)

# Ajout d'interactions
add1(lm_complet, .~. +.^2, test = "F")
lm_complet <- update(lm_complet, .~.+ white:ratio2)
add1(lm_complet, .~. +.^2, test = "F")
lm_complet <- update(lm_complet, .~.+ meanses:white)
add1(lm_complet, .~. +.^2, test = "F")

ols_vif_tol(lm_complet)

ols_plot_resid_fit(lm_complet)


Y_hat <- predict(lm_complet)
MSE_lm2 <- mean((math - Y_hat)^2)


data_residuals <- data %>% mutate(
   residus_lm = rstudent(lm_complet),
   id = row_number())


plot_residuals_id()
#' Il semble y avoir quelques valeurs abérantes, mais, globalement, les résidus
#' sont bien répartis autour de zéro.
plot_residuals_variable(meanses)
plot_residuals_variable(homework)
plot_residuals_variable(white)
plot_residuals_variable(ratio)
plot_residuals_variable(ratio2)
#' On voit que l'effet aléatoire (meanses|schid) pourrait expliquer la variation
#' des résidus d'une grappe à une autre dans le modèle linéaire.


#--- Entraînement des modèles linéaires mixtes (lmm) complets ----
lm_complet

rm(lmm_VC_UN, lmm_VC_UN1, lmm_CS_UN)
lmm_VC_UN <- lmer(
   math ~ meanses + homework + white + ratio + ratio2 + 
      white:ratio2 + meanses:white +
      (homework|schid),
   REML = TRUE
   )
lmm_VC_UN1 <- lmer(
   math ~ meanses + homework + white + ratio + ratio2 + 
      white:ratio2 + meanses:white +
      (homework||schid),
   REML = TRUE
)
lmm_CS_UN <- lme(# Modèle mixte avec structure de variance CS, UN
   fixed = math ~ meanses + homework + white + ratio + ratio2 + 
      white:ratio2 + meanses:white,
   random = ~ homework | schid,
   data = data,
   method='REML',
   correlation = corCompSymm()
)

# Comparaison des AIC
tbl_AIC <- rbind(
   "VC UN" = extractAIC(lmm_VC_UN),
   "VC UN(1)" = extractAIC(lmm_VC_UN1),
   "CS UN" = extractAIC_lme(lmm_CS_UN)
)
colnames(tbl_AIC) <- list("dl", "AIC")
xtable::xtable(tbl_AIC)
tbl_AIC

lmm_complet <- lmm_VC_UN
summary(lmm_complet)
#' On a la même structure de variance qu'en a).

#--- Test des effets aléatoires ----
lmm_H0 <- update(lmm_complet, .~.- (homework | schid) + (1 | schid), REML=T)
lrtest_lmm(lmm_H0, lmm_complet)
# On rejette H0 et on conserve l'effet aléatoire (homework|schid)

lmm_Q1b <- lmm_complet
summary(lmm_Q1b)


#--- Sélection des effets fixes ----
car::Anova(lmm_Q1b, type=3)
lmm_Q1b <- update(lmm_Q1b, .~.- white:ratio2 )
car::Anova(lmm_Q1b, type=3)
lmm_Q1b <- update(lmm_Q1b, .~.- ratio2)
car::Anova(lmm_Q1b, type=3)
lmm_Q1b <- update(lmm_Q1b, .~.- ratio)
car::Anova(lmm_Q1b, type=3)
lmm_Q1b <- update(lmm_Q1b, .~.- meanses:white)
car::Anova(lmm_Q1b, type=3)
lmm_Q1b <- update(lmm_Q1b, .~.+ meanses:homework )
car::Anova(lmm_Q1b, type=3)
lmm_Q1b <- update(lmm_Q1b, .~.- meanses:homework )
summary(lmm_Q1b)

sum_lmm_Qst1b <- summary(lmm_Q1b)
beta <- sum_lmm_Qst1b$coefficients[,1]
sd_beta <- sum_lmm_Qst1b$coefficients[,2]
tbl_betas <- cbind(
   "Estimateurs"=beta,
   "Écarts-types"=sd_beta,
   beta + cbind(rep(-1, 4), rep(1, 4)) * qnorm(0.975) * sd_beta
)
xtable::xtable(tbl_betas)

sum_lmm_Qst1a$varcor$schid %*% diag(nrow=2)
sum_lmm_Qst1a$sigma^2


data %>% mutate(  # Visualisation des résidus pour le modèle mixte entraîné
   predictions = Y_hat,
   residus = rstudent(lmm_Q1a)
) %>%
   ggplot() +
   geom_point(aes(
      x=predictions,
      y=residus
   ), color=schid) +
   geom_hline(yintercept = 0)
#' Les résidus sont beaucoup mieux répartis. Les problèmes soulevés en (a) ont
#' été réglés.
data %>% arrange(schid) %>%
   mutate(  # Visualisation des résidus pour le modèle mixte entraîné
      residus = rstudent(lmm_Q1a),
      id = row_number()
   ) %>%
   ggplot() +
   geom_point(aes(
      x=id,
      y=residus
   ), color=schid) +
   geom_hline(yintercept = 0)


# Évaluation des résultats finaux pour la question 1 ---------------------------
Y_hat <- predict(lmm_Q1b)
MSE_lmm_Q1b <- mean((math - Y_hat)^2)

MSE_Qst1 <- matrix(c(MSE_lm, MSE_lm2, MSE_lmm_Q1a, MSE_lmm_Q1b), ncol=2)
dimnames(MSE_Qst1) <- list(
   list("Question 1a", "Question 1b"), list("Modèle linéaire", "Modèle mixte"))
xtable::xtable(MSE_Qst1, digits=2, align="lcc",
               caption="Écarts quadratiques moyens")
MSE_Qst1
#' L'ajout de la variable meanses diminue le besoin de faire appel aux effets
#' aléatoires.
detach(data)

#============================= Question 2 ======================================
rm(list=ls())
detach(data)
data <- read.table(file='GirlsGrowth.dat', header = T)
data %>%  glimpse()
data <- 
   data %>%
   mutate(child = factor(child),
          group = factor(group),
          temps=age-6
   )
data %>%  glimpse()
attach(data)


data %>% arrange(child, age) %>%
   ggplot() +
   geom_line(aes(
      x = age,
      y = height,
      colour = child
   ), alpha = 0.6) +
   ylim(c(110, 150)) +
   ylab('Grandeur')
#' On voit que la relation entre l'âge des petites filles et leur grandeur est
#' très linéaire.


#--- Entraînement d'un modèle linéaire standard ----
lm_complet <- lm(height ~ group + temps + temps:group,
                 data=data)
car::vif(lm_complet)

Y_hat <- predict(lm_complet)
MSE_lm <- mean((height - Y_hat)^2)


#--- Analyse des résidus de régression ---
data_residuals <- data %>% 
   mutate(
      residus_lm = rstudent(lm_complet),
      id = row_number()
      ) %>%
   arrange(child, temps)


plot_residuals_variable <- function(variable){
   #' Fonction qui affiche le graphique des résidus en fonction d'une variable
   #' explicative
   #' @param variable La variable explicative que l'on désire regarder 
   data_residuals %>% 
      ggplot() +
      geom_point(aes(
         x = variable,
         y = residus_lm,
         colour = child
      ), alpha = 0.6) +
      geom_hline(yintercept = 0) +
      xlab(deparse(substitute(variable))) +
      ylab("Résidus")
}
plot_residuals_id <- function() {
   #' Fonction qui affiche le graphique des résidus en fonction des indices i.
   n_row = nrow(data_residuals)
   data_residuals %>% 
      ggplot() +
      geom_line(aes(
         x = id,
         y = residus_lm,
         colour = child
      ), alpha = 0.6) +
      geom_hline(yintercept = 0) +
      xlab("i") +
      ylab("Résidus")
}


plot_residuals_id()
plot_residuals_variable(group)
data_residuals %>% 
   ggplot() +
   geom_line(aes(
      x = temps,
      y = residus_lm,
      colour = child
   ), alpha = 0.6) +
   geom_hline(yintercept = 0) +
   xlab('temps') +
   ylab("Résidus")
#' On remarque que les courbes des résidus varient selon l'enfant.


#--- Entraînement des modèles linéaires mixtes (lmm) complets ----
lm_complet

lmm_VC_UN <- lmer(# Modèle mixte avec structure de variance VC, UN
   height ~ group + temps + group:temps +
      (group | child),
   data = data, 
   REML = TRUE
   )
lmm_CS_UN <- lme(# Modèle mixte avec structure de variance CS, UN
   fixed = height ~ group + temps + group:temps,
   random = ~ group | child,
   data = data,
   method='REML',
   correlation = corCompSymm()
   )
lmm_AR1_UN <- lme(# Modèle mixte avec structure de variance AR(1), UN
   fixed = height ~ group + temps + group:temps,
   random = ~ group | child,
   data = data,
   method='REML',
   correlation = corAR1()
)


extractAIC_lme <- function(object) {
   #' Fonction qui prend un modèle de la classe lme et qui calcule l'AIC
   #' @param object Modèle entraîné de la classe lme.
   df = sum(object$dims$ncol, object$dims$qvec) + 1
   aic = -2 * (object$logLik - df)
   if (round(AIC(object), 4) != round(aic, 4)){
      print(paste(AIC(object), aic))
      warning('Le nombre de degrés de liberté spécifié est incorrect.')
      }
   return(c(df, aic))
}


# Comparaison des AIC
tbl_AIC <- rbind(
   "VC UN" = extractAIC(lmm_VC_UN),
   "CS UN" = extractAIC_lme(lmm_CS_UN),
   "AR(1) UN" = extractAIC_lme(lmm_AR1_UN)
)
colnames(tbl_AIC) <- list("dl", "AIC")
xtable::xtable(tbl_AIC)
tbl_AIC

summary(lmm_AR1_UN)

#' La structure de variance pour les résidus qui semble la plus adaptée selon 
#' le critère de l'AIC est la structure auto-régressive d'ordre 1 (AR1).
lmm_complet <- lmm_AR1_UN

#--- Test des effets aléatoires ----
lrtest_lmm <- function(model_H0, model_H1) {
   #' Fonction qui calcule le test du rapport des vraisemblances pour tester 
   #' la nécessité des effets aléatoires.
   #' @param model_H0 Le modèle sous l'hypothèse nulle (modèle simple)
   #' @param model_H1 Le modèle sous H1 (modèle complet)
   test <- lmtest::lrtest(model_H0, model_H1)
   df <- test$Df[2]
   chisq <- test$Chisq[2]
   pval <- 0.5 * (2 - pchisq(chisq, df - 1) - pchisq(chisq, df))
   print(paste("chisq:", round(chisq, 2), "Df:", df, "p-value:", round(pval,3)),
         quote=F)
}


lmm_H0 <- lme(
   fixed = height ~ group + temps + group:temps,
   random = ~ 1 | child,
   data = data,
   method='REML',
   correlation = corAR1()
)
lrtest_lmm(lmm_H0, lmm_complet)
lmm_complet <- lmm_H0
#' On ne peut rejeter H0; ce qui veut dire que l'on retire l'effet aléatoire
#' group | child

lmm_H1 <- lme(
   fixed = height ~ group + temps + group:temps,
   random = ~ 1 | child,
   data = data,
   method='ML',
   correlation = corAR1()
)

lmm_H0 <- lm(
   height ~ group + temps + group:temps,
   data = data
)
lrtest_lmm(lmm_H0, lmm_H1)
#' On rejette H0 et on conserve l'effet aléatoire (1 | child)
#' On en déduit que le modèle mixte est significativement plus efficace qu'un
#' modèle linéaire standard pour ce contexte.
lmm_Q2 <- lmm_complet


#--- Sélection des effets fixes ----
car::Anova(lmm_Q2, type=3)

# Résultats finaux pour la question 2 ---------------------------
sum_lmm_Qst2 <- summary(lmm_Q2)
beta <- sum_lmm_Qst2$tTable[, 1]
sd_beta <- sum_lmm_Qst2$tTable[, 2]
tbl_betas <- cbind(
   "Estimateurs"=beta,
   "Écarts-types"=sd_beta,
   beta + cbind(rep(-1, length(beta)),
                rep(1, length(beta))) * qnorm(0.975) * sd_beta
)
xtable::xtable(tbl_betas)

sum_lmm_Qst2
sum_lmm_Qst2$sigma^2

rho <- 0.9498482
sd_residus <- 2.987839
Vi <- sd_residus * matrix(
   c(
      rho^(0:4),
      rho, rho^(0:3),
      rho^(2:1), rho^(0:2),
      rho^(3:0), rho,
      rho^(4:0)
   ), ncol=5
   )
round(Vi, 3)


Y_hat <- predict(lmm_Q2)
(MSE_lmm_Q2 <- mean((height - Y_hat)^2))
MSE_lm

data %>% 
   mutate(
      residus_lm = residuals(lmm_Q2),
      id = row_number()
   ) %>%
   arrange(child, age) %>% 
   ggplot() +
   geom_line(aes(
      x = temps,
      y = residus_lm,
      colour = child
   ), alpha = 0.6) +
   geom_hline(yintercept = 0) +
   xlab('temps') +
   ylab("Résidus")

detach(data)


#============================= Question 3 ======================================
rm(list=ls())
library(gee)
# Lecture des donnees

data <- read.table("bolusdata_1.txt", header = TRUE)
data %>% glimpse()
data <- data %>% mutate(
   group = factor(group)
   )
# a) ---------------------------------------------------------------------------
#---- Structures de correlation de travail ----

# Non structuré (UN) 
#' (ici on peut l'utiliser car l'indice j signifie la meme chose pour toutes
#' les grappes)
model_UN <-
   gee(
      count ~ group + time + group:time,
      id = id,
      family = poisson(link = "log"),
      corstr = "unstructured",
      data = data
   )
summary(model_UN)
#' La corrélation semble s'atténuer au fil du temps. Une structure AR(1) serait
#' donc appropriée.

# Et avec la structure AR(1)
model_AR1 <-
   gee(
      count ~ group + time + group:time,
      id = id,
      family = poisson(link = "log"),
      corstr = "AR-M",
      Mv = 1,
      data = data
   )
summary(model_AR1)

# Structure d'indépendance
model_IND <-
   gee(
      count ~ group + time + group:time,
      id = id,
      family = poisson(link = "log"),
      corstr = "independence",
      data = data
   )
summary(model_IND)

# Structure échangeable
model_EX <-
   gee(
      count ~ group + time + group:time,
      id = id,
      family = poisson(link = "log"),
      corstr = "exchangeable",
      data = data
   )
summary(model_EX)

#' Choix du model: On suppose que la drogue un impact à court terme. 
#' De ce fait, l'impacte se degrade dans le temps. 
#' Selon cette hypothese le model avec AR(1) serait le plus approprie.
#' Les administrations de drogue les plus recents sont les plus importants et 
#' sont donc plus correllés.
model_selectionne <- model_AR1
#---- Analyse des variables exogènes ----

selection_variable_GEE <- function(model, L) {
   #' Fonction qui effectue le test chi-carré pour vérifier que les estimateurs
   #' ne sont pas égaux à zéro.
   #' H0: Le(s) paramètre(s) spécifiée(S) dans la matrice L sont égaux à zéro.
   #' H1: Le(s) paramètre(s) spécifiée(S) dans la matrice L ne sont pas égaux
   #' à zéro.
   #' @param model Un modèle de type GEE pré-entraîné.
   #' @param L Matrice de test de dimension nb_test x nb_variables
   #' @return le seuil observé du test.
   coefs <- model$coefficients
   Vs <- model$robust.variance
   r <- if (purrr::is_null(nrow(L))) 1 else nrow(L)
   d <- if (r == 1) 0 else rep(0, r)
   Lbeta <- L %*% coefs
   if (r == 1) L <- t(L)
   chi2 <- t(Lbeta - d) %*% solve(L %*% Vs %*% t(L)) %*% (Lbeta - d)
   p_value <- 1 - pchisq(chi2, r)
   print(paste("p-value:", round(p_value, 3), "Df:", r), quote=FALSE)
}


summary(model_selectionne)$coefficients
L <- matrix(c(
   0, 1, 0, 0,
   0, 0, 1, 0,
   0, 0, 0, 1), 
   nrow=3, byrow = T)
selection_variable_GEE(model_selectionne, L)
#' On rejette l'hypothèse nulle que tous les coefficients de la régression sont
#' égaux à zéro. Le modèle est donc utile.

L <- c(0,0,0,1)
selection_variable_GEE(model_selectionne, L)
#' On ne peut pas rejeter H0 et on conclu que l'interaction group:time n'est pas
#' utile.
model_selectionne <- update(model_selectionne, .~.- group:time)

L <- c(0, 0, 1)
selection_variable_GEE(model_selectionne, L)
L <- c(0, 1, 0)
selection_variable_GEE(model_selectionne, L)
#' Les deux autres variables sont significatives au seuil de 5%

model_final <- model_selectionne
summary(model_final)


# b) ---------------------------------------------------------------------------
# Estimation ponctuelle
L <- c(1, 1, 5)
c(t(L) %*% model_final$coefficients)
# I.C. a 95%
c(t(L) %*% model_final$coefficients) + c(-1, 1) *
   1.96 * sqrt(c(t(L) %*% model_final$robust.variance %*% L))
# pour E(Y):
exp(c(t(L) %*% model_final$coefficients))
exp(c(t(L) %*% model_final$coefficients) +
       c(-1, 1) * 1.96 * sqrt(c(t(L) %*% model_final$robust.variance %*% L)))

# I.C. effet moyen de la dose
L <- c(0, 1, 0)
c(t(L) %*% model_final$coefficients)
# I.C. a 95%
c(t(L) %*% model_final$coefficients) + c(-1, 1) *
   1.96 * sqrt(c(t(L) %*% model_final$robust.variance %*% L))
# pour E(Y):
exp(c(t(L) %*% model_final$coefficients))
exp(c(t(L) %*% model_final$coefficients) +
       c(-1, 1) * 1.96 * sqrt(c(t(L) %*% model_final$robust.variance %*% L)))

#' La dose de 2mg augmente en moyenne le nombre d'administrations de 35% avec
#' IC entre 6% a 74%
