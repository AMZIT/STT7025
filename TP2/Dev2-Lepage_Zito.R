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


current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))

#============================= Question 1 ======================================
rm(list=ls())
