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


# Analyse préliminaire
par(mfrow=c(3,3))
for (l in colnames(data[-16])) {
   plot(data[,l], data$B, xlab=l, ylab="B")
   plot(I(log(data[,l])), data$B, xlab=l, ylab="B")
   plot(I(sqrt(data[,l])), data$B, xlab=l, ylab="B")
}
par(mfrow=c(1,1))
#' Il semblerait que les variables A12, A13 et A14 gagneraient à recevoir une 
#' transformation logarithmique afin d'améliorer la relation linéaire les 
#' unissant à la variable B

f <- as.formula('B~. -A12 -A13 - A14 + I(log(A12)) + I(log(A13)) +
                I(sqrt(A14))')
modele_complet <-  lm(f, data=data, x=TRUE, y=TRUE)
AIC(modele_complet)

f <- as.formula('B~.-A9 -A12 -A13 - A14 + I(log(A12)) + I(log(A13)) +
                I(sqrt(A14)) + I(sqrt(A9))')
modele_complet <-  lm(f, data=data, x=TRUE, y=TRUE)
AIC(modele_complet)

f <- as.formula('B~.-A9 -A12 -A13 - A14 + I(log(A12)) + I(log(A13)) +
                I(log(A14)) + I(sqrt(A9))')
modele_complet <-  lm(f, data=data, x=TRUE, y=TRUE)
AIC(modele_complet)

f <- as.formula('B~.-A9 -A12 -A13 - A14 + I(log(A12)) + I(log(A13)) +
                I(sqrt(A14)) + I(log(A9))')
modele_complet <-  lm(f, data=data, x=TRUE, y=TRUE)
AIC(modele_complet)

f <- as.formula('B~.-A9 -A12 -A13 - A14 + I(log(A12)) + I(log(A13)) +
                I(log(A14)) + I(log(A9))')
modele_complet <-  lm(f, data=data, x=TRUE, y=TRUE)
AIC(modele_complet)

# Meilleur modèle:
f <- as.formula('B~.-A9 -A12 -A13 - A14 + I(log(A12)) + I(log(A13)) +
                I(log(A14)) + I(sqrt(A9))')
modele_complet <-  lm(f, data=data, x=TRUE, y=TRUE)


# Analyse de la multicolinéarité -----------------------------------------------
ols_vif_tol(modele_complet)
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


multicol_diagnosis(modele_complet)
modele_complet2 <- update(modele_complet, .~.-A5)
ols_vif_tol(modele_complet2) %>% top_n(1, VIF)

multicol_diagnosis(modele_complet2)
modele_complet2 <- update(modele_complet2, .~.-A3)
ols_vif_tol(modele_complet2) %>% top_n(1, VIF)

multicol_diagnosis(modele_complet2)
modele_complet2 <- update(modele_complet2, .~.- A7)
ols_vif_tol(modele_complet2) %>% top_n(1, VIF)

multicol_diagnosis(modele_complet2)
modele_complet2 <- update(modele_complet2, .~.- A6)
ols_vif_tol(modele_complet2) %>% top_n(1, VIF)

multicol_diagnosis(modele_complet2)
modele_complet2 <- 
   update(modele_complet2, .~.- A4 - I(sqrt(A9)) + I(A4 + sqrt(A9)))
ols_vif_tol(modele_complet2) %>% top_n(1, VIF)

multicol_diagnosis(modele_complet2)
modele_complet2 <- update(modele_complet2, .~.- I(A4 + sqrt(A9)))
ols_vif_tol(modele_complet2) %>% top_n(1, VIF)

multicol_diagnosis(modele_complet2)
modele_complet2 <- update(modele_complet2, .~.-A10 -A15 + I((A10 + A15)/2))
ols_vif_tol(modele_complet2) %>% top_n(1, VIF)

multicol_diagnosis(modele_complet2)
modele_complet2 <- update(modele_complet2, .~.- I((A10 + A15)/2))
ols_vif_tol(modele_complet2) %>% top_n(1, VIF)

multicol_diagnosis(modele_complet2)
modele_complet2 <- 
   update(modele_complet2, .~.-I(log(A12)) -I(log(A13)) + I(log(sqrt(A12*A13))))
ols_vif_tol(modele_complet2) %>% top_n(1, VIF)
#' Le problème de multicolinéarité est maintenant réglé.
drop1(modele_complet2)


# --- Autre approche: régression régularisée (LASSO) ---
set.seed(2020)
modele_glmnet_lasso <- glmnet(modele_complet$x[,-1], modele_complet$y, 
                              family='gaussian', alpha=1)
#' On va prendre la valeur de lambda qui conserve le plus de variables 
#' explicatives possibles et qui minimise la statistique de déviance.


tester_lambda <- function(modele, df, return_model=FALSE) {
   #' Fonction qui extrait la valeur de lambda minimisant la déviance pour un 
   #' nombre de degrés de liberté donné. Par la suite, elle présente le 
   #' VIF maximal afin de valider que la régression régularisée a permis de 
   #' traiter la multicollinéarité.
   #' @param modele Modèle glm.net sur lequel on cherche à régler un problème de
   #' multicollinéarité;
   #' @param df Nombre de degrés de libertés utilisé pour sélectionner le 
   #' meilleur lambda;
   #' @param return_model Variable booléenne. Si TRUE, retourne le modèle testé.
   lambda <- modele$lambda[which(modele$df==df)]
   dev_ratio <- modele$dev.ratio[which(modele$df==df)]
   
   meilleur_lambda <- lambda[which.min(dev_ratio)]
   coefs <- coef(modele, s = meilleur_lambda)
   
   predictors <- rownames(coefs)[which(coefs != 0)][-1]
   f <- as.formula(paste('B~', paste(predictors, collapse ='+')))
   modele_LASSO <- lm(f, data=data)
   
   print(paste("Valeur de lambda =", meilleur_lambda))
   print(ols_vif_tol(modele_LASSO) %>% top_n(1, VIF))
   if(return_model) return(modele_LASSO)
}


tester_lambda(modele_glmnet_lasso, df=15)
tester_lambda(modele_glmnet_lasso, df=14)
tester_lambda(modele_glmnet_lasso, df=13)
tester_lambda(modele_glmnet_lasso, df=11)
tester_lambda(modele_glmnet_lasso, df=9)
tester_lambda(modele_glmnet_lasso, df=8)
#' On règle le problème de multicolinéarité lorsque le nombre maximal de 
#' variables explicatives est de 8.
modele_filtre_Lasso <- tester_lambda(modele_glmnet_lasso, df=8, return_model=T)

AIC(modele_complet2, modele_filtre_Lasso)
ols_regress(modele_complet2)$prsq
ols_regress(modele_filtre_Lasso)$prsq
#' Le modèle utilisant la régression LASSO est plus performant selon l'AIC et le
#' R-carré de prédiction.


# Sélection de variables -------------------------------------------------------
# --- Test de tous les sous-modèles possibles ---
all_possible <- ols_step_all_possible(modele_filtre_Lasso)
plot(all_possible)
(best_models <- all_possible %>% as_tibble() %>% top_n(3, predrsq))
#' Deux modèles se démarquent du lot. Le plus simple des deux (Index 382) 
#' est celui qui minimise . Comme l'objectif de cette
#' question est de faire de la prédiction, on favorisera le R-carré de PRESS
#' par rapport aux autres critères.
best_models[1,3]
mortality_model <- lm(B ~ A1 +A2 +A6 +A8 +I(log(A13)) +I(sqrt(A9)), data = data)
ols_regress(mortality_model)


# Ajout d'interactions ---------------------------------------------------------
add1(mortality_model, .~. +.^2 , test = "F")
#' Aucune interaction intéressante au seuil de 1% 


# Appréciation du modèle -------------------------------------------------------
mortality_model %>% ols_regress()
#' Le R-carré de prédiction est de 0.715; le pouvoir prédictif de ce modèle est
#' tout à fait adéquat considérant les données disponibles.


# Réponse à la question 1  -----------------------------------------------------
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
cc <- c("age", "trestbps","chol","thalach" ,"oldpeak","ca"  )
# Analyse préliminaire
par(mfrow=c(3,3))
for (l in cc) {
   boxplot( as.formula(paste0(l," ~ Y")),data = data)
   boxplot( as.formula(paste0("I(log(",l,")) ~ Y")),data = data)
   boxplot( as.formula(paste0("I(sqrt(",l,")) ~ Y")),data = data)
}
par(mfrow=c(1,1))


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

# Analyse Preliminaire
cc <- c("Exposure" ,"VehValue" )
par(mfrow=c(2,3))
for (l in cc) {
   boxplot( as.formula(paste0(l," ~ ClaimNb")),data = data)
   boxplot( as.formula(paste0("I(log(",l,")) ~ ClaimNb")),data = data)
   boxplot( as.formula(paste0("I(sqrt(",l,")) ~ ClaimNb")),data = data)
}
par(mfrow=c(1,1))

# Transformations ? ------------------------------------------------------------
modele.GAM <- gam::gam(ClaimNb ~ VehValue + offset(log(Exposure)),
                  data = data, family = poisson("log"))
plot(modele.GAM) 
#' Aucune transformation nécessaire.


# Analyse de la surdispersion --------------------------------------------------
# Modèle de Poisson
modele_Poisson <-  glm(ClaimNb ~ .-Exposure, offset = log(Exposure),
                       data = data, family = poisson("log"))
modele_Poisson %>% summary()
(phi <- modele_Poisson$deviance / (modele_Poisson$df.residual))
#' Il semble y avoir sous-dispersion

# Modèle binomial négative
modele_BinomNeg <- MASS::glm.nb(ClaimNb ~ .-Exposure + offset(log(Exposure)), 
                          data = data, link="log")
lmtest::lrtest(modele_Poisson, modele_BinomNeg)


AIC(modele_Poisson, modele_BinomNeg)
BIC(modele_Poisson, modele_BinomNeg)

l0 <- logLik(modele_Poisson)
l1 <- logLik(modele_BinomNeg)
0.5 *(1 - pchisq(2*(l1 - l0), 1)) %>% as.numeric()
#' Avec un test du ratio de vraisemblance, on trouve que le modèle binomial 
#' négative est significativement mieux ajusté aux données qu'un modèle 
#' poissonien. Les statistiques de l'AIC et du BIC confirment ce constat.

# Sélection de variables -------------------------------------------------------
modele_nul <- MASS::glm.nb(ClaimNb ~ 1 + offset(log(Exposure)), 
                           data = data, link="log")
modele_complet <- modele_BinomNeg

# --- Tester tous les sous-modèles ---
# all_models <- glmulti::glmulti(modele_complet)
#' Nécessite l'installation de JavaScript pour être utilisé...

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

