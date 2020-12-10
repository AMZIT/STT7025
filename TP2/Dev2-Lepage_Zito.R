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
          schid = factor(schid)
          )
data %>%  glimpse()
attach(data)


#------------------------- Analyse graphique -----------------------------------
data %>%
   ggplot() +
   geom_boxplot(aes(
      x = factor(homework),
      y = math
   ), alpha = 0.6)
#' Il semble exister une relation non linéaire entre le nombre d'heures
#' d'étude et les résultats en math.

data %>%
   ggplot() +
   geom_boxplot(aes(
      x = factor(ratio),
      y = math
   ), alpha = 0.6)
#' La variable ratio semble avoir une relation significative, mais il est 
#' difficile de voir si celle-ci est linéaire.

data %>%
   ggplot() +
   geom_boxplot(aes(
      x = schid,
      y = homework
   ), alpha = 0.6)
#' Le nombre d'heures d'étude des étudiants varie énormément selon
#' l'école d'appartenance

data %>%
   ggplot() +
   geom_boxplot(aes(
      x = schid,
      y = math
   ), alpha = 0.6)
#' Les résultats scolaires aussi.


# a) ---------------------------------------------------------------------------
#--- Entraînement d'un modèle linéaire standard ----
lm_complet <- lm(
   math ~ (homework + white + ratio)^2,
   data=data)
ols_plot_resid_fit(lm_complet)
#' Le modèle est biaisé. Les résidus ne sont pas centrés à zéro.
#' Voyons voir si des transformations sur les xi pourrait améliorer le modèle.

library(gam)
gam_mod <- gam(math ~ (s(homework) + white + s(ratio))^2)
plot(gam_mod, se=TRUE)
#' Possiblement que transformer ratio pour (ratio-18)^2 pourrait améliorer le
#' modèle. En ce qui attrait à homework, celle-ci peut demeurer telle quelle.
summary(gam_mod)
# Les deux transformations sont significatives au seuil de 1%.
data %>% mutate(  # Visualisation des résidus pour le GAM
   Y_hat = predict(gam_mod),
   residus = rstudent(gam_mod)
   ) %>%
   ggplot() +
   geom_point(aes(
         x=Y_hat,
         y=residus
      ), color='blue') +
   geom_hline(yintercept = 0)
# Problème non réglé.


# Analyse BoxCox.
MASS::boxcox(lm_complet)
#' Comme lambda = 1 appartient à l'intervalle de confiance, on n'apporte aucune
#' transformation à Y.


lm_complet <- lm(
   math ~ homework + I((ratio - 18)^2) + white,
   data=data
)


# Analyse de multicolinéarité
ols_vif_tol(lm_complet)


# Ajout d'interactions
add1(lm_complet, .~. +.^2, test = "F")
lm_complet <- update(lm_complet, .~.+ I((ratio - 18)^2):white)
add1(lm_complet, .~. +.^2, test = "F")


ols_plot_resid_fit(lm_complet)
#' Ce n'est pas bon, mais c'est le mieux que je puisse faire avec les variables
#' disponibles.
#' 
#' Voyons voir si un modèle mixte pourrait expliquer ce biais.


# Calcul du MSE en entraînement.
Y_hat <- predict(lm_complet)
MSE_lm <- mean((math - Y_hat)^2)


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
#' Il est difficile de tirer des conclusions de ces graphiques.
#' Les tests du ratio des vraisemblances permettra de mettre les choses
#' au clair sur les effets aléatoires.


#--- Entraînement des modèles linéaires mixtes (lmm) complets ----
# des matrices de variance appropriées
detach(data)
data <- data %>% mutate(ratio2 = (ratio-18)^2)
attach(data)
lm_complet

lmm_VC_UN <- lmer(  # Modèle mixte avec structure de variance VC, UN
   math ~ homework + ratio2 + white + ratio2:white +
      (homework | schid) + (ratio2 | schid),
   data=data,
   REML = TRUE
   )
lmm_VC_UN1 <- lmer(  # Modèle mixte avec structure de variance VC, UN(1)
   math ~ homework + ratio2 + white + ratio2:white +
      (homework || schid) + (ratio2 || schid),
   data=data,
   REML = TRUE
)
lmm_CS_UN <- lme(# Modèle mixte avec structure de variance CS, UN
   fixed = math ~ homework + ratio2 + white + ratio2:white,
   random = ~ homework + ratio2 | schid,
   data = data,
   method='REML',
   correlation = corCompSymm()
)
lmm_AR1_UN <- lme(# Modèle mixte avec structure de variance AR(1), UN
   fixed = math ~ homework + ratio2 + white + ratio2:white,
   random = ~ homework + ratio2 | schid,
   data = data,
   method='REML',
   correlation = corAR1()
)


extractAIC_lme <- function(object) {
   #' Fonction qui prend un modèle de la classe lme et qui calcule l'AIC
   #' @param object Modèle entraîné de la classe lme.
   df = sum(object$dims$ncol, object$dims$qvec)
   aic = -2 * (object$logLik - df)
   return(c(df, aic))
}


# Comparaison des AIC
extractAIC(lmm_VC_UN)
extractAIC(lmm_VC_UN1)
extractAIC_lme(lmm_CS_UN)
extractAIC_lme(lmm_AR1_UN)


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

lmm_H0 <- update(lmm_complet, .~.- (ratio2 | schid), REML=T)
lrtest_lmm(lmm_H0, lmm_complet)
lmm_complet <- lmm_H0
#' On rejete H0 au seuil de 5% et on conserve l'effet aléatoire (ratio2 | schid)
#' On conserve donc tous les effets aléatoires.

lmm_H0 <- update(lmm_complet, .~.- (homework | schid) + (1 | schid), REML=T)
lrtest_lmm(lmm_H0, lmm_complet)
#' On rejette H0 et on conserve l'effet aléatoire (homework | schid)

lmm_Q1a <- lmm_complet
summary(lmm_Q1a)


#--- Sélection des effets fixes ----
car::Anova(lmm_Q1a, type=3)


# Calcul du MSE en entraînement.
Y_hat <- predict(lmm_Q1a)
(MSE_lmm_Q1a <- mean((math - Y_hat)^2))
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
   ), color=schid) +
   geom_hline(yintercept = 0)
#' Le prolème a changé: Les résidus sont linéaires, mais on a un problème
#' d'hétéroscédasticité.
data %>% mutate(  # Visualisation des résidus pour le modèle mixte entraîné
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
lm_complet <- lm(math ~ (meanses + homework + white + ratio2)^2, data=data)
ols_plot_resid_fit(lm_complet)
#' Mêmes constats qu'en a)

ols_vif_tol(lm_complet)
# Il y a présence de multicolinéarité


multicol_diagnosis <- function(modele, TOL=0.5, TeX_table=FALSE) {
   #' Fonction qui cible les variables problématiques
   #' @param modele Un modèle linéaire sur lequel diagnostiquer le problème
   #' de multicolinéarité;
   #' @param TOL Seuil de tolérance pour déterminer que la proportion de 
   #' variance soulève un problème de multicolinéarité;
   #' @param TeX_table Variable booléenne. Si TRUE, affiche un tableau synthèse
   #' en format LaTeX.
   multicol_diagnosis <- ols_eigen_cindex(modele)
   eigen_index <- rownames(multicol_diagnosis)
   multicol_diagnosis <- 
      multicol_diagnosis %>%  
      mutate(j = eigen_index, .after=1) %>%
      select(-Eigenvalue, -intercept) %>%
      top_n(2, `Condition Index`) %>%
      pivot_longer(-c(1,2), values_to='p_lj') %>%
      filter(p_lj>=TOL)
   
   if (TeX_table)
      multicol_diagnosis %>% xtable::xtable() %>% print(include.rownames=FALSE)
   return(multicol_diagnosis)
}


multicol_diagnosis(lm_complet)
lm_complet2 <- update(lm_complet, .~.-homework:white)
ols_vif_tol(lm_complet2) %>% top_n(1, VIF)

lm_complet <- lm_complet2


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
#' On voit que l'effet aléatoire (meanses|schid) pourrait expliquer la variation
#' des résidus d'une grappe à une autre dans le modèle linéaire.


#--- Entraînement des modèles linéaires mixtes (lmm) complets ----
# des matrices de variance appropriées
lm_complet

rm(lmm_VC_UN, lmm_VC_UN1)
lmm_VC_UN <- lmer(
   math ~ meanses + homework + white + ratio2 + meanses:homework + 
      meanses:white +
      (meanses|schid) + (ratio2|schid) + (homework|schid),
   REML = TRUE
   )
lmm_complet_UN1 <- lmer(
   math ~ meanses + homework + white + ratio2 + meanses:homework + 
      meanses:white +
      (meanses||schid) + (ratio2||schid) + (homework||schid),
   REML = TRUE
)

extractAIC(lmm_complet_UN)
extractAIC(lmm_complet_UN1)
#' On a la même structure de variance qu'en a).

lmm_complet <- lmm_complet_UN
summary(lmm_complet)


#--- Test des effets aléatoires ----

lmm_H0 <- update(lmm_complet, .~.- (homework | schid), REML=T)
lrtest_lmm(lmm_H0, lmm_complet)
# On rejette H0 et on conserve l'effet aléatoire (homework|schid)

lmm_H0 <- update(lmm_complet, .~.- (ratio2 | schid), REML=T)
lrtest_lmm(lmm_H0, lmm_complet)
lmm_complet <- lmm_H0
# On ne peut rejeter H0 au seuil de 5% et on retire l'effet (ratio|schid)

lmm_H0 <- update(lmm_complet, .~.- (meanses | schid), REML=T)
lrtest_lmm(lmm_H0, lmm_complet)
lmm_complet <- lmm_H0
# On ne peut rejeter H0 au seuil de 5% et on retire l'effet (meanses|schid)

lmm_Q1b <- lmm_complet
summary(lmm_Q1b)


#--- Sélection des effets fixes ----
car::Anova(lmm_Q1b, type=3)
lmm_Q1b <- update(lmm_Q1b, .~.- meanses:homework)
car::Anova(lmm_Q1b, type=3)
lmm_Q1b <- update(lmm_Q1b, .~.- ratio2)
car::Anova(lmm_Q1b, type=3)
lmm_Q1b <- update(lmm_Q1b, .~.- meanses:white)
car::Anova(lmm_Q1b, type=3)

summary(lmm_Q1b)


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
          age=age-6
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
   ylim(c(110, 150))
#' On voit que la relation entre l'âge des petites filles et leur grandeur est
#' très linéaire.


# Entraînement d'un modèle linéaire standard
lm_complet <- lm(height ~ group + age + age:group,
                 data=data)
ols_vif_tol(lm_complet)

Y_hat <- predict(lm_complet)
MSE_lm <- mean((height - Y_hat)^2)


#--- Analyse des résidus de régression ---
data_residuals <- data %>% 
   mutate(
      residus_lm = rstudent(lm_complet),
      id = row_number()
      ) %>%
   arrange(child, age)


plot_residuals_variable <- function(variable){
   #' Fonction qui affiche le graphique des résidus en fonction d'une variable
   #' explicative
   #' @param variable La variable explicative que l'on désire regarder 
   data_residuals %>% 
      ggplot() +
      geom_line(aes(
         x = age,
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
plot_residuals_variable(age)
#' On remarque que les courbes des résidus varient selon l'enfant.


# --- Entraînement des modèles linéaires mixtes (lmm) complets et sélection ---
# des matrices de variance appropriées
lm_complet

lmm_VC_UN <- lmer(# Modèle mixte avec structure de variance VC, UN
   height ~ group + age + group:age +
      (age | child),
   data = data, 
   REML = TRUE
   )
lmm_VC_UN1 <- lmer(# Modèle mixte avec structure de variance VC, UN(1)
   height ~ group + age + group:age +
      (age || child), 
   data = data, 
   REML = TRUE
   )
lmm_CS_UN <- lme(# Modèle mixte avec structure de variance CS, UN
   fixed = height ~ group + age + group:age,
   random = ~ age | child,
   data = data,
   method='REML',
   correlation = corCompSymm()
   )
lmm_AR1_UN <- lme(# Modèle mixte avec structure de variance AR(1), UN
   fixed = height ~ group + age + group:age,
   random = ~ 1 | child,
   data = data,
   method='REML',
   correlation = corAR1()
)


extractAIC_lme <- function(object) {
   #' Fonction qui prend un modèle de la classe lme et qui calcule l'AIC
   #' @param object Modèle entraîné de la classe lme.
   df = sum(object$dims$ncol, object$dims$qvec)
   aic = -2 * (object$logLik - df)
   return(c(df, aic))
}

# Comparaison des AIC
extractAIC(lmm_VC_UN)
extractAIC(lmm_VC_UN1)
extractAIC_lme(lmm_CS_UN)
extractAIC_lme(lmm_AR1_UN)
#' La structure de variance pour les résidus qui semble la plus adaptée selon 
#' le critère de l'AIC est la structure auto-régressive d'ordre 1 (AR1).
#' 
#' Comme cette méthode ne converge que si l'effet aléatoire ne concerne que
#' l'ordonnée à l'origine, on ne complexifiera pas ce dernier davantage.
#' Conséquemment, comme il n'y a qu'un seul effet aléatoire, la structure de
#' variance de ce dernier n'a aucune importance puique q=1 et que Di = sigma^2.
lmm_complet <- lmm_AR1_UN
summary(lmm_complet)


#--- Test des effets aléatoires ---
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


lmm_H1 <- lme(
   fixed = height ~ group + age + group:age,
   random = ~ 1 | child,
   data = data,
   method='ML',
   correlation = corAR1()
)

lmm_H0 <- lm(
   height ~ group + age + group:age,
   data = data
)
lrtest_lmm(lmm_H0, lmm_H1)
#' On rejette H0 et on conserve l'effet aléatoire (1 | child)
#' On en déduit que le modèle mixte est significativement plus efficace qu'un
#' modèle linéaire standard pour ce contexte.
lmm_Q2 <- lmm_complet


#--- Sélection des effets fixes
car::Anova(lmm_Q2, type=3)
# Toutes les variables sont significatives


Y_hat <- predict(lmm_Q2)
(MSE_lmm_Q2 <- mean((height - Y_hat)^2))
MSE_lm
#' L'ajout de l'effet aléatoire n'a apporté aucune amélioration de la prédiction
#' sur les données d'entraînement.
