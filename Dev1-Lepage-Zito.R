#-----------------------------------------------
# TRAVAIL PRATIQUE 1
# 
# Cours: STT-7120
# Professeur: M. Thierry Duchesne
# 
# Équipe 7
# Étudiants: - Alexandre Lepage 111 144 776
#            - Amedeo Zito
#-----------------------------------------------

library(tidyverse)
library(olsrr)
library(glmnet)
library(plotmo)
library(CASdatasets)

current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))

#============================= Question 1 ======================================
rm(list=ls())

data <- read.table("pollution_1.txt", row.names=1)
#      A1,  the average annual precipitation;
#      A2,  the average January temperature;
#      A3,  the average July temperature;
#      A4,  the size of the population older than 65;
#      A5,  the number of members per household;
#      A6,  the number of years of schooling for persons over 22;
#      A7,  the number of households with fully equipped kitchens;
#      A8,  the population per square mile; 
#      A9,  the size of the nonwhite population;
#      A10, the number of office workers;
#      A11, the number of families with an income less than $3000;
#      A12, the hydrocarbon pollution index;
#      A13, the nitric oxide pollution index;
#      A14, the sulfur dioxide pollution index;
#      A15, the degree of atmospheric moisture.
#      B,   the death rate.
colnames(data) <- colnames_ <- c(paste0('A', 1:15), 'B')
data %>% glimpse()
data %>% summary()


modele_complet <-  lm(B~., data=data, x=TRUE, y=TRUE)


# Analyse de la multicolinéarité -----------------------------------------------
ols_vif_tol(modele_complet)
#' Il y a une sérieuse présence de multicolinéarité


multicol_diagnosis <- function(modele, tol=0.6, verbose=TRUE){
   #' Fonction qui cible les variables redondantes dans un modèle linéaire
   #' en évaluant les proportions de variabilité (p_lj) en présence de 
   #' multicolinéarité..
   #' 
   #' Variables:
   #' @param modele Le modèle linéaire sur lequel on veut procéder au diagnostique;
   #' @param tol Seuil sur lequel on base le rejet d'une variable exogène;
   #' @param verbose Afficher les informations pertinentes.
   diagnosis <- ols_coll_diag(modele)$eig_cindex %>% as_tibble() %>%
      select('Condition Index', starts_with('A')) %>%
      arrange(desc(`Condition Index`)) %>%
      filter(`Condition Index`>30) %>% 
      pivot_longer(starts_with("A"), values_to = 'P_lj') %>% 
      filter(P_lj>=tol)
   
   if (verbose){
      print(diagnosis)
   }
   return(diagnosis[1,2])
}


retirer_variables_redondantes <- function(modele, tol=0.6, verbose=TRUE) {
   #' Fonction qui itère sur la fonction @multicol_diagnosis tant et aussi 
   #' longtemps que le problème de multicolinéarité n'est pas réglé (VIF>10).
   #' 
   #' Variables:
   #' @param modele Le modèle linéaire sur lequel on veut procéder au diagnostique;
   #' @param tol Seuil sur lequel on base le rejet d'une variable exogène;
   #' @param verbose Afficher les informations pertinentes.
   vif_max <-  ols_vif_tol(modele) %>% select(VIF) %>% max()
   modele_complet2 <- modele
   while (vif_max>10) 
      {
      if (verbose) print(paste("VIF max:", round(vif_max)), quote=F)
      
      variable_redondante <- multicol_diagnosis(modele_complet2, tol, verbose)
      
      if (is.null(variable_redondante))
      {
         message(
            "Le paramètre tol est trop élevé et ne permet pas de régler
              le problème de multicolinéarité.")
      }
      else
      {
         f <- as.formula(paste(".~. -", paste(variable_redondante)))
         modele_complet2 <-  update(modele_complet2, f)
      }
      
      vif_max <- ols_vif_tol(modele_complet2) %>% select(VIF) %>% max()
   }
   if (verbose) print(paste("VIF max:", round(vif_max,2)), quote=F)
   return(modele_complet2)
}


modele_complet2 <- retirer_variables_redondantes(modele_complet, 0.6)
#' Le problème de multicolinéarité est maintenant réglé.
modele_complet2 %>% summary()


# --- Autre approche: régression régularisée (LASSO) ---
modele_glmnet <- glmnet(modele_complet$x[,-1], modele_complet$y, 
                        family='gaussian', alpha=1)
plot_glmnet(modele_glmnet)
set.seed(2020)
cv_out <- cv.glmnet(modele_complet$x[,-1], modele_complet$y, 
                    family='gaussian', alpha=1)
plot(cv_out)

coefs <- coef(modele_glmnet, s = c(cv_out$lambda.min,
                                    cv_out$lambda.1se))
colnames(coefs) <- c("lambda.min", "lambda.1se")
coefs
predictors <- rownames(coefs)[which(coefs[, 1] != 0)][-1]
#' Variables à conserver : A1, A2, A6, A7, A8, A9 et A14

f <- as.formula(paste('B~', paste(predictors, collapse ='+')))
modele_complet2 <- lm(f, data=data)
ols_vif_tol(modele_complet2)
#' La multicolinéarité a été réglée.


#' L'approche itérative en utilisant les proportions de variation ne donne pas
#' les mêmes résultats que celle utilisant la régression régularisée.
#' Il reste à voir selon l'AIC laquelle des deux approches offre le meilleur 
#' pouvoir prédictif.


# Transformations ? ------------------------------------------------------------
car::boxCox(modele_complet2)
#' Aucune transformation nécessaire.


# Sélection de variables -------------------------------------------------------
# --- Approche de Hosmer et Lemeshow ---

# --- Test de tous les sous-modèles possibles ---
all_possible <- ols_step_all_possible(modele_complet2)
plot(all_possible)
(bests_models <- all_possible %>% as_tibble() %>% top_n(3, predrsq))

#' Deux modèles se démarquent du lot. Le plus simple des deux (Index 256) 
#' est celui qui minimise le R-carré de PRESS. Comme l'objectif de cette
#' question est de faire de la prédiction, on favorisera ce critère par rapport
#' aux autres.
predictors <- str_split(bests_models[1, 3], ' ') %>% unlist()
f <- as.formula(paste("B~", paste(predictors, collapse = '+')))
mortality_model <- lm(f, data=data)


# Ajout d'interactions ---------------------------------------------------------
add1(mortality_model, .~. +.^2 , test = "F")
#' Aucune interaction intéressante au seuil de 1%


# Appréciation du modèle -------------------------------------------------------
mortality_model %>% ols_regress()
#' Le R-carré de prédiction est de 0.64; le pouvoir prédictif de ce modèle est
#' donc très limité.
ame_model <- lm(B~A1 + A2 + A3 + A6 + A8 + A9 + A14, data=data)
anova(mortality_model, ame_model)

#' Comparaison de l'AIC selon la méthode de traitement de la multicolinéarité:
# AIC(mortality_model)
#' Méthode itérative: AIC: 603.0711,    Pred R-Squared 0.640
#' Régression LASSO : AIC: 598.7574,    Pred R-Squared 0.665
#' 
#' Verdict: Pour traiter la multicolinéarité, la régression LASSO permet de
#' conserver un meilleur pouvoir prédictif.


# Réponse à la question 1  ------------------------------------------------------------------
new_data <- read.table(text="40 30 80 9 3 10 77 4100 13 46 15 25 26 145 55") %>%
   as_tibble()
colnames(new_data) <- colnames_[-16]

predict.lm(mortality_model, newdata = new_data,
               type='response', interval = 'prediction')


#============================= Question 2 ======================================
rm(list=ls())

data <- read.table("processed.cleveland.data", sep=',') %>% as_tibble()
#' @param age : âge en années
#' @param sex : 1 = homme, 0 = femme
#' @param cp : nature des douleurs à la poitrine, variable qualitative à 4 modalités, où 1 dénote l’angine
#' typique, 2 l’angine atypique, 3 une douleur non anginienne et 4 une douleur asymptomatique
#' @param trestbps : tension artérielle au repos (en mm Hg) à l’admission à l’hôpital
#' @param chol : cholestérol sanguin en mg/dl
#' @param fbs : indicatrice qui vaut 1 si le taux de sucre sanguin à jeun > 120 mg/dl et qui vaut 0 sinon
#' @param restecg : résultat de l’électrocardiogramme au repos, variable qualitative à 3 modalités, où 0
#' signifie normal, 1 signifie anomalie des ondes ST-T et 2 signifie hypertrophie probable du ventricule
#' gauche
#' @param thalach : pouls maximum atteint
#' @param exang : indicatrice indiquant la présence d’angine induite par l’exercice (1 pour oui, 0 pour non)
#' @param oldpeak : baisse dans ST induite par l’exercise par rapport au repos
#' @param slope : pente du segment de ST lors de l’exercice maximal, variable qualitative à 3 modalités
#' soit 1 pour ascendante, 2 pour plate et 3 pour descendante
#' @param ca : nombre de vaissaux sanguins majeurs colorés par fluroscopie
#' @param thal : variable qualitative à 3 modalités où 3 = normal, 6 = défaut réparé, 7 = défaut réparable
#' @param num : la variable réponse que nous cherchons à prédire est Y = 1 si num> 0 et Y = 0 si num= 0
names(data) <- c("age", "sex", "cp", "trestbps", "chol", "fbs", "restecg",
                 "thalach", "exang", "oldpeak", "slope", "ca", "thal", "num")
data %>% glimpse()
data %>% summary()
data <- data %>% mutate(
   sex = as.factor(sex),
   cp = as.factor(cp),
   fbs = as.factor(fbs),
   restecg = as.factor(restecg),
   exang = as.factor(exang),
   slope = as.factor(slope),
   ca = as.integer(na_if(ca, '?')),
   thal = as.factor(as.integer(na_if(thal, '?'))),
   Y = as.factor(1*(num>0)),
   num = NULL
) %>% drop_na()


# Analyse de la multicolinéarité -----------------------------------------------

#' Comme on n'a pas besoin d'un modèle pour tester la multicolinéarité, on va
#' se créer un modèle linéaire bidon qui permettra d'utiliser les outils du
#' package olsrr.
bidon <- lm(rnorm(nrow(data))~., data=data)
ols_vif_tol(bidon)
#' Il n'y a aucun signe de multicolinéarité.


# Transformations ? ------------------------------------------------------------
modele.GAM <- gam::gam(Y ~ ., data = data, family = binomial())
plot(modele.GAM) 
#' Aucune transformation nécessaire.


# Sélection de variables -------------------------------------------------------
modele_nul <- glm(Y ~ 1, family = binomial("logit"), data = data)
modele_complet <-  glm(Y ~ ., data = data, family = binomial("logit"), x=T, y=T)

# --- Test de tous les sous-modèles ---
# f <- as.formula(paste("Y~", paste(colnames(data[,-14]), collapse = '+')))
# sick_tout <- glmbb::glmbb(f,
#                           family=binomial("logit"),
#                           criterion='AIC',
#                           cutoff = 30,
#                           data=data)
#' La fonction glmbb du package du même nom, n'est pas fait pour gérer les 
#' variables catégorielles. Il faut donc utiliser des méthodes alternatives.

# --- Méthode forward ---
sick_forward <- MASS::stepAIC(
   modele_nul,
   scope = list(upper = modele_complet, lower = modele_nul),
   trace = 0,
   direction = "forward",
   data = data,
   k = 2
)
sick_forward$anova
sick_forward %>% drop1(test="LRT")
sick_forward %>% rsq::rsq(adj=T)
#' Différentes valeurs de k ont été testées puisque, avec k=2, ce n'est pas
#' toutes les variables explicatives qui sont significatives au seuil de 5%.
#' Les valeurs de k testées sont k=2,3,4. Pour chacune d'elle, la statistique
#' du R-carré ajusté a été calculée.
#' 
#' k=2 -> Adj R-Squared = 0.5824615
#' k=3 -> Adj R-Squared = 0.5792386
#' k=4 -> Adj R-Squared = 0.5550493
#' 
#' Verdict: k=2 est l'hyperparamètre qui maximise l'explicabilité de la variable
#' endogène par les variables exogènes. Ainsi, dans un contexte d'explication,
#' il est préférable d'augmenter le seuil du test de Wald.


# --- Méthode backward ---
sick_backward <- MASS::stepAIC(
   modele_complet,
   trace = 0,
   direction = "backward",
   data = data,
   k = 2
)
sick_backward$anova
sick_backward %>% drop1(test="LRT")
sick_backward %>% rsq::rsq(adj=T)
#' Même constat
#' k=2 -> Adj R-Squared = 0.5824615
#' k=3 -> Adj R-Squared = 0.5792386
#' k=4 -> Adj R-Squared = 0.5732573


# --- Méthode stepwise ---
sick_stepwise <- MASS::stepAIC(
   modele_nul, 
   scope=list(upper=modele_complet, lower=modele_nul),
   trace = 0,
   direction = "both",
   data = data,
   k = 2
)
sick_stepwise$anova
sick_stepwise %>% drop1(test="LRT")
sick_stepwise %>% rsq::rsq(adj=T)
#' Idem
#' k=2 -> Adj R-Squared = 0.5824615
#' k=3 -> Adj R-Squared = 0.5792386
#' k=4 -> Adj R-Squared = 0.5550493


AIC(sick_forward, sick_backward, sick_stepwise)
rsq::rsq(sick_forward, adj=T)
rsq::rsq(sick_backward, adj=T)
rsq::rsq(sick_stepwise, adj=T)
#' On trouve donc que les trois algorithmes de sélection de variables 
#' aboutissent au même modèle.
sick_model <- sick_backward


# Ajout d'interactions ---------------------------------------------------------
#' Approche stepwise avec seuil d'inclusion à 5% et seuil d'exclusion à 1% pour 
#' les interactions et à 5% pour les variables.
sick_model %>% add1(.~. +.^2 , test = "LRT")
sick_model <- sick_model %>% update(.~. + ca:thal)

sick_model %>% drop1(test="LRT")
sick_model <- sick_model %>% update(.~. - exang)

sick_model %>% add1(.~. +.^2 , test = "LRT")
sick_model <- sick_model %>% update(.~. + cp:slope )

sick_model %>% drop1(test="LRT")
sick_model <- sick_model %>% update(.~. - thalach)

sick_model %>% drop1(test="LRT")
sick_model <- sick_model %>% update(.~. - cp:slope )

sick_model %>% add1(.~. +.^2 , test = "LRT")
sick_model %>% drop1(test="LRT")
#' Seule l'interaction  ca:thal est significative au seuil de 1%.
sick_model %>% rsq::rsq(adj=T) # 0.5875979
#' Il y a une légère amélioration du R-carré ajusté par rapport au modèle sans
#' interaction.


# Réponse à la question 2 ------------------------------------------------------
sick_model %>% anova()
sick_model %>% coef()



#============================= Question 3 ======================================
rm(list=ls())

data(ausprivauto0405)
data <- ausprivauto0405 %>% as_tibble()
rm(ausprivauto0405)

data %>% glimpse()
data <- data %>% select(-ClaimOcc, -ClaimAmount)
data %>% summary()

# Analyse de la multicolinéarité -----------------------------------------------
bidon <- lm(rnorm(nrow(data))~., data=data)
ols_vif_tol(bidon)
#' On voit que certaines valeurs de VIF sont plus grandes que 10. Cependant,
#' celles-ci sont calculées sur toutes les variables catégorielles unitaires.
#' On veut donc connaître le VIF agrégé. 
car::vif(bidon)


# Transformations ? ------------------------------------------------------------
modele.GAM <- gam::gam(ClaimNb ~ VehValue + offset(log(Exposure)),
                  data = data, family = poisson("log"))
plot(modele.GAM) 
#' Aucune transformation nécessaire.


# Analyse de la variance --------------------------------------------------
# Modèle de Poisson
modele_Poisson <-  glm(ClaimNb ~ .-Exposure, offset = log(Exposure),
                       data = data, family = poisson("log"))
modele_Poisson %>% summary()
(phi <- modele_Poisson$deviance / (modele_Poisson$df.residual))
#' Il semble y avoir sous-dispersion

# Modèle binomial négative
modele_BinomNeg <- MASS::glm.nb(ClaimNb ~ .-Exposure + offset(log(Exposure)), 
                          data = data, link="log")


AIC(modele_Poisson, modele_BinomNeg)
BIC(modele_Poisson, modele_BinomNeg)
0.5 *(1 - pchisq(modele_Poisson$deviance - modele_BinomNeg$deviance, 1))
#' Avec un test du ratio de vraisemblance, on trouve que le modèle binomial 
#' négative est significativement mieux ajusté aux données qu'un modèle 
#' poissonien. Les statistiques de l'AIC et du BIC confirment ce constat.

# Sélection de variables -------------------------------------------------------
modele_nul <- MASS::glm.nb(ClaimNb ~ 1 + offset(log(Exposure)), 
                           data = data, link="log")
modele_complet <- modele_BinomNeg

# --- Méthode forward ---
freq_forward <- MASS::stepAIC(
   modele_nul,
   scope = list(upper = modele_complet, lower = modele_nul),
   trace = 0,
   direction = "forward",
   data = data,
   k = 2
)
freq_forward$anova
freq_forward %>% drop1(test="LRT")
freq_forward %>% rsq::rsq(adj=T)
#' Toutes les variables sélectionnées sont significatives
#' Le R-carré ajusté est ridiculement petit...quelque chose cloche...
#' possiblement que la fonction rsq est non compatible avec la fonction glm.nb.


# --- Méthode backward ---
freq_backward <- MASS::stepAIC(
   modele_complet, 
   scope = list(upper = modele_complet, lower = modele_nul),
   trace = 0,
   direction = "backward",
   data = data,
   k = 2
)
freq_backward$anova
freq_backward %>% drop1(test="LRT")
freq_backward %>% rsq::rsq(adj=T)
#' Même constat


# --- Méthode stepwise ---
freq_stepwise <- MASS::stepAIC(
   modele_nul, 
   scope=list(upper=modele_complet, lower=modele_nul),
   trace = 0,
   direction = "both",
   data = data,
   k = 2
)
freq_stepwise$anova
freq_stepwise %>% drop1(test="LRT")
freq_stepwise %>% rsq::rsq(adj=T)
#' Idem


AIC(freq_forward, freq_backward, freq_stepwise)
#' On trouve que les trois algorithmes de sélection de variables 
#' aboutissent au même modèle.
freq_model <- freq_forward


# Ajout d'interactions ---------------------------------------------------------
freq_model %>% add1(.~. +.^2 , test = "LRT")
# Aucune interaction significative.


# Réponse à la question 3 ------------------------------------------------------
freq_model %>% anova()
freq_model %>% coef()
