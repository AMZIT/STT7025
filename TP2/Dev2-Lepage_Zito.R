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


# a) ---------------------------------------------------------------------------
lm_complet <- lm(math ~ (homework + white + ratio)^2, data=data)
ols_vif_tol(lm_complet)
#' Il y a présence de multicolinéarité


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
lm_complet2 <- update(lm_complet, .~.-white:ratio)
ols_vif_tol(lm_complet2) %>% top_n(1, VIF)

multicol_diagnosis(lm_complet2)
lm_complet2 <- update(lm_complet2, .~.-homework:ratio)
ols_vif_tol(lm_complet2) %>% top_n(1, VIF)

multicol_diagnosis(lm_complet2)
lm_complet2 <- update(lm_complet2, .~.-homework:white)
ols_vif_tol(lm_complet2) %>% top_n(1, VIF)

lm_complet <- lm_complet2


Y_hat <- predict(lm_complet)
MSE_lm <- mean((math - Y_hat)^2)

#--- Analyse des résidus de régression ---
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
#' Dépendamment de la grappe, il semblerait que les résidus aient une espérance
#' et une variance qui diffèrent.

plot_residuals_variable(homework)
plot_residuals_variable(white)
plot_residuals_variable(ratio)
#' Il est difficile de tirer des conclusions de ces graphiques.
#' Les tests du ratio des vraisemblances permettra de mettre les choses
#' au clair sur les effets aléatoires.


# --- Entraînement des modèles linéaires mixtes (lmm) complets et sélection ---
# des matrices de variance appropriées
lm_complet
set.seed(2020)
lmm_complet_UN <- lmer(math ~ homework + white + ratio +
                          (ratio|schid) + (white|schid) + (homework|schid), 
                       REML = TRUE)
lmm_complet_UN1 <- lmer(math ~ (homework + white + ratio)^2 +
                          (ratio||schid) + (white||schid) + (homework||schid), 
                        REML = TRUE)
extractAIC(lmm_complet_UN)
extractAIC(lmm_complet_UN1)

lmm_complet <- lmm_complet_UN
summary(lmm_complet)
#' On choisit la structure nonstructurée (UN) pour la matrice de variance D.
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


lmm_H0 <- update(lmm_complet, .~.- (homework | schid), REML=T)
lrtest_lmm(lmm_H0, lmm_complet)
#' On rejette H0 et on conserve l'effet aléatoire (homework | schid)

lmm_H0 <- update(lmm_complet, .~.- (white | schid), REML=F)
lrtest_lmm(lmm_H0, lmm_complet)
#' On rejete H0 au seuil de 5% et on conserve l'effet aléatoire (white | schid)

lmm_H0 <- update(lmm_complet, .~.- (ratio | schid), REML=F)
lrtest_lmm(lmm_H0, lmm_complet)
#' On rejete H0 au seuil de 5% et on conserve l'effet aléatoire (ratio | schid)
#' On conserve donc tous les effets aléatoires.

lmm_Q1a <- lmm_complet
summary(lmm_Q1a)


#--- Sélection des effets fixes
car::Anova(lmm_Q1a, type=3)
lmm_Q1a <- update(lmm_Q1a, .~.- ratio - (ratio | schid), REML=F)
# set.seed(2020)
# lmm_Q1a <- lmer(
#    math ~ homework + white + (white | schid) + (homework | schid),
#    REML = F
# )
car::Anova(lmm_Q1a, type=3)

summary(lmm_Q1a)


Y_hat <- predict(lmm_Q1a)
(MSE_lmm_Q1a <- mean((math - Y_hat)^2))
MSE_lm
#' Belle amélioration de l'erreur quadratique
 

# b) ---------------------------------------------------------------------------
#--- Analyse des résidus de régression ---
lm_complet <- lm(math ~ (meanses + homework + white + ratio)^2, data=data)
ols_vif_tol(lm_complet)
# Il y a présence de multicolinéarité

multicol_diagnosis(lm_complet)
lm_complet2 <- update(lm_complet, .~.-white:ratio)
ols_vif_tol(lm_complet2) %>% top_n(1, VIF)

multicol_diagnosis(lm_complet2)
lm_complet2 <- update(lm_complet2, .~.-homework:ratio)
ols_vif_tol(lm_complet2) %>% top_n(1, VIF)

multicol_diagnosis(lm_complet2)
lm_complet2 <- update(lm_complet2, .~.-homework:white)
ols_vif_tol(lm_complet2) %>% top_n(1, VIF)

multicol_diagnosis(lm_complet2)
lm_complet2 <- update(lm_complet2, .~.-meanses:ratio)
ols_vif_tol(lm_complet2) %>% top_n(1, VIF)

lm_complet <- lm_complet2


Y_hat <- predict(lm_complet)
MSE_lm2 <- mean((math - Y_hat)^2)


data_residuals <- data %>% mutate(
   residus_lm = rstudent(lm_complet),
   id = row_number())


plot_residuals_id()
#' Dépendamment de la grappe, il semblerait que les résidus aient une espérance
#' et une variance qui diffèrent.
plot_residuals_variable(meanses)
plot_residuals_variable(homework)
plot_residuals_variable(white)
plot_residuals_variable(ratio)
#' On voit que l'effet aléatoire (meanses|schid) pourrait expliquer la variation
#' des résidus d'une grappe à une autre dans le modèle linéaire.


# --- Entraînement des modèles linéaires mixtes (lmm) complets et sélection ---
# des matrices de variance appropriées
lm_complet

rm(lmm_complet_UN, lmm_complet_UN1)
set.seed(2020)
lmm_complet_UN <- lmer(
   math ~ meanses + homework + white + ratio + meanses:homework + 
      meanses:white + (meanses|schid) + (ratio|schid) + (white|schid) +
      (homework|schid),
   REML = TRUE
   )
lmm_complet_UN1 <- lmer(
   math ~ meanses + homework + white + ratio + meanses:homework + 
      meanses:white + (meanses||schid) + (ratio||schid) + (white||schid) +
      (homework||schid),
   REML = TRUE
)
extractAIC(lmm_complet_UN)
extractAIC(lmm_complet_UN1)
#' On a la même structure de variance qu'en a).

lmm_complet <- lmm_complet_UN
summary(lmm_complet)


#--- Test des effets aléatoires ---

set.seed(2020)
lmm_H0 <- update(lmm_complet, .~.- (homework | schid), REML=T)
lrtest_lmm(lmm_H0, lmm_complet)
# On rejette H0 et on conserve l'effet aléatoire (homework|schid)

set.seed(2020)
lmm_H0 <- update(lmm_complet, .~.- (white | schid), REML=T)
lrtest_lmm(lmm_H0, lmm_complet)
lmm_complet <- lmm_H0
# On ne peut rejeter H0 au seuil de 5% et on retire l'effet (white|schid)

set.seed(2020)
lmm_H0 <- update(lmm_complet, .~.- (ratio | schid), REML=T)
lrtest_lmm(lmm_H0, lmm_complet)
lmm_complet <- lmm_H0
# On ne peut rejeter H0 au seuil de 5% et on retire l'effet (ratio|schid)

set.seed(2020)
lmm_H0 <- update(lmm_complet, .~.- (meanses | schid), REML=T)
lrtest_lmm(lmm_H0, lmm_complet)
# On ne peut rejeter H0 au seuil de 5% et on retire l'effet (meanses|schid)

lmm_Q1b <- lmm_H0
summary(lmm_Q1b)


#--- Sélection des effets fixes
car::Anova(lmm_Q1b, type=3)
lmm_Q1b <- update(lmm_Q1b, .~.- meanses:homework)
car::Anova(lmm_Q1b, type=3)
lmm_Q1b <- update(lmm_Q1b, .~.- ratio)
car::Anova(lmm_Q1b, type=3)
lmm_Q1b <- update(lmm_Q1b, .~.- meanses:white)
car::Anova(lmm_Q1b, type=3)

summary(lmm_Q1b)


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
