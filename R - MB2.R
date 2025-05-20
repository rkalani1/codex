library('haven'); library('dplyr'); library('purrr'); library('car'); library('tidyverse'); library('broom'); library('officer'); library('flextable'); 
library('kableExtra'); library('knitr'); library('corrr'); library('corrplot'); library('ggcorrplot'); library('GGally'); library('sandwich')
library('glmnet'); library('caret'); library('tibble')

#Olink3K Data
MESAOlink3K<-read_csv("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/Proteomics_Olink/Intensity_Data/SMP_IntensityNormalized_03062024.csv")
MESAOlink3K_ProteinKeys<-read_csv("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/Proteomics_Olink/Mapping_Olink3K_OlinkID_Key_03062024.csv")
MESAOlink3K_LinkingData<-read_csv("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/Proteomics_Olink/Mapping_SMP_Plate_20240514.csv")
Olink3KFinal<-merge.data.frame(MESAOlink3K,MESAOlink3K_LinkingData,by="SampleID") #MESAOlink3K data merged by TOM_ID
Olink3Ke6<-Olink3KFinal[Olink3KFinal$exam==6,] #Olink3K data for Exam 6
Olink3Ke5<-Olink3KFinal[Olink3KFinal$exam==5,] #Olink3K data for Exam 5
Olink3Ke1<-Olink3KFinal[Olink3KFinal$exam==1,] #Olink3K data for Exam 1
#Rename "sidno" column to "subject_id" for merging with covariate data
Olink3Ke6<-Olink3Ke6%>%rename(subject_id=sidno)
Olink3Ke5<-Olink3Ke5%>%rename(subject_id=sidno)
Olink3Ke1<-Olink3Ke1%>%rename(subject_id=sidno)

#Covariate Data
CoV_e6<-read.delim("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/Data from MESA CC/SHARe_Exam6Main/SHARe_Exam6Main.txt")
CoV_e5<-read.delim("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/Data from MESA CC/SHARe_MesaExam5Main/SHARe_MesaExam5Main_DS.txt")
CoV_e1<-read.delim("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/Data from MESA CC/SHARe_MesaExam1Main/SHARe_MesaExam1Main_DS.txt")
Edu_e1<-CoV_e1%>%dplyr::select(subject_id, educ1) #Education data for Exam 1
Edu_e1<-Edu_e1%>%
  mutate(Edu = case_when(
    educ1 %in% 0:2 ~ 0,
    educ1 %in% 3:8 ~ 1
  )) # Education data recoded to 0=less than high school education, 1= completed high school
Edu_e1<-Edu_e1[!is.na(Edu_e1$Edu), ] #Remove NA education values
CoV_e6<-CoV_e6[!is.na(CoV_e6$sbp6c), ] #Remove NA SBP values
CoV_e6<-CoV_e6[!is.na(CoV_e6$htnmed6c), ] #Remove NA HTN med values
CoV_e6<-CoV_e6[!is.na(CoV_e6$cig6c), ] #Remove NA cig use values
CoV_e6<-CoV_e6[!is.na(CoV_e6$ldl6), ] #Remove NA LDLC values
CoV_e6<-CoV_e6[!is.na(CoV_e6$afib6), ] #Remove NA AF Prevalent values
ApoE<-read.delim("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/Data from MESA CC/SHARe_AncilMesaApoE/SHARe_AncilMesaApoE.txt")
ApoE<-ApoE%>%
  mutate(E4 = case_when(
    ApoE == 22 ~ 0,
    ApoE == 23 ~ 0,
    ApoE == 24 ~ 1,
    ApoE == 33 ~ 0,
    ApoE == 34 ~ 1,
    ApoE == 44 ~ 1
  )) # ApoE data recoded to E4 carrier status (0=non-carrier, 1=carrier)
ApoE<-ApoE[!is.na(ApoE$E4), ] #Remove NA ApoE values
PrevalentCVD<-read.delim("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/Data from MESA CC/SHARe_EventsThruYear2019/SHARe_EventsThruYear2019_DS.txt")

#MRI Outcome Data
MB<-read.delim("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/Data from MESA CC/SHARe_AncilMesaAF_BMRIMB/SHARe_AncilMesaAF_BMRIMB_DS.txt")
MB$mb_present<-ifelse(MB$mb_n_total>0, 1, 0)

#Merged Dataframes for Analyses
Prot_CoV<-merge.data.frame(Olink3Ke6,CoV_e6,by="subject_id") #Merge Olink3K and Covariate Data
Prot_CoV_ApoE<-merge.data.frame(Prot_CoV,ApoE,by="subject_id") #Merge Olink3K, Covariate, and ApoE Data
Prot_CoV_ApoE_Edu<-merge.data.frame(Prot_CoV_ApoE,Edu_e1,by="subject_id") #Merge Olink3K, Covariate, ApoE and Education Data
Prot_CoV_ApoE_Edu_CVD<-merge.data.frame(Prot_CoV_ApoE_Edu,PrevalentCVD,by="subject_id") #Merge Olink3K, Covariate, ApoE, Education and Prevalent CVD Data
Prot_CoV_ApoE_Edu_CVD_MB<-merge.data.frame(Prot_CoV_ApoE_Edu_CVD,MB,by="subject_id") #Merge Olink3K, Covariate, ApoE, Education, and MB Data
MBFinala<-Prot_CoV_ApoE_Edu_CVD_MB[Prot_CoV_ApoE_Edu_CVD_MB$qsm_swi_image_quality != 4,] #QC (qsm_swi_image_quality=4) exclusion
MBFinal<-MBFinala[!rownames(MBFinala) %in% c("NA", "NA.1"), ] #Remove NA rows

### MB Analyses ###

#MB Analyses - Baseline characteristics
mean(MBFinal$age6c); sd(MBFinal$age6c); summary(is.na(MBFinal$age6c)) # Age Mean and SD
table(MBFinal$gender1); summary(is.na(MBFinal$gender1)) # 0=female, 1=male
table(MBFinal$race1c); table(is.na(MBFinal$race1c)) # 1=white, 2=chinese, 3=black, #4=hispanic/latino
table(MBFinal$site6c); summary(is.na(MBFinal$site6c)) # 3=WFU, 4=COL, 5=JHU, 6=UMN, 7=NWU, 8=UCLA
table(MBFinal$Edu); summary(is.na(MBFinal$Edu)) # 0=less than high school education, 1= completed high school
mean(MBFinal$cepgfr6c); sd(MBFinal$cepgfr6c); summary(is.na(MBFinal$cepgfr6c)) # CKD-EPI eGFR Mean and SD (mL/min/1.73m2)
mean(MBFinal$bmi6c); sd(MBFinal$bmi6c); summary(is.na(MBFinal$bmi6c)) # BMI Mean and SD (kg/m2)
mean(MBFinal$sbp6c); sd(MBFinal$sbp6c); summary(is.na(MBFinal$sbp6c)) # SBP Mean and SD (mmHg)
table(MBFinal$htnmed6c); summary(is.na(MBFinal$htnmed6c)) # 1=Yes, 0=No
table(MBFinal$cig6c); summary(is.na(MBFinal$cig6c)) # 2=Current, 1=Former, 0=Never
mean(MBFinal$ldl6); sd(MBFinal$ldl6); summary(is.na(MBFinal$ldl6)) # LDL Mean and SD (mg/dL)
table(MBFinal$E4); summary(is.na(MBFinal$E4)) # 1=Presence of 1 or more E4 Allele, 0=Absence of E4 Allele
table(MBFinal$afib6); summary(is.na(MBFinal$afib6)); MBFinal<-MBFinal %>% mutate(AFprevalent = case_when( 
  afib6 == 0 ~ 0,
  afib6 == 1 ~ 1,
  afib6 == 9 ~ 0,
)); table(MBFinal$AFprevalent); summary(is.na(MBFinal$AFprevalent)) # 1=Prevalent AF, 0=No
table(MBFinal$dm036t); summary(is.na(MBFinal$dm036t)) # 1=Prevalent DM, 0=No
table(MBFinal$mi); summary(is.na(MBFinal$mi)) # 1=Prevalent MI, 0=No
table(MBFinal$chf); summary(is.na(MBFinal$chf)) # 1=Prevalent CHF, 0=No

#MB Analyses
table(MBFinal$mb_present); summary(is.na(MBFinal$mb_present)) # 1=MB Present, 0=No

#MB Analyses - Unadjusted
results<-data.frame(Protein = character(), RR = numeric(), CI_lower = numeric(), CI_upper = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)
for(proteins in colnames(MBFinal)[3:2943]) {
  formula <- as.formula(paste("mb_present ~", proteins))
  model <- glm(formula, family = poisson(link = "log"), data=MBFinal)
  robust_se<-sqrt(diag(vcovHC(model, type = "HC0")))
  RR<-exp(coef(model)[2])
  CI<-confint(model, level=0.95)[2,] 
  CI_adj <- exp(coef(model)[2] + qnorm(c(0.025, 0.975)) * robust_se[2])
  z_value<-coef(model)[2]/robust_se[2]
  p_value<-2*pnorm(abs(z_value), lower.tail = FALSE)
  results<-rbind(results, data.frame(Protein = proteins, RR = RR, CI_lower = CI[1], CI_upper = CI[2], P_Value = p_value))
}
results$P_Value_Adjusted<-p.adjust(results$P_Value, method = "BH")
significant_results<-results%>%filter(P_Value_Adjusted<0.05)
final_results<-merge(significant_results, MESAOlink3K_ProteinKeys, by.x = "Protein", by.y = "OlinkID")
final_table<-final_results%>%dplyr::select(ID = Assay, RR, P_Value_Adjusted)%>%arrange(P_Value_Adjusted); print(final_table)
ft<-flextable(final_table); ft<-set_table_properties(ft, width = 1, layout = "autofit"); ft<-padding(ft, padding.top = 0, padding.bottom = 0)
doc<-read_docx("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/results draft.docx"); doc<-body_add_flextable(doc, value=ft); print(doc, target = "/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/results draft.docx")

#Fully Adjusted LR Model
CovariatesE6<-c("age6c", "gender1", "race1c", "site6c", "Edu", "cepgfr6c", "bmi6c", "sbp6c", "htnmed6c", "cig6c", "ldl6", "E4", "AFprevalent", "dm036t", "mi", "chf")
results<-data.frame(Protein = character(), RR = numeric(), CI_lower = numeric(), CI_upper = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)
for(proteins in colnames(MBFinal)[3:2943]) {
  formula <- as.formula(paste("mb_present ~", proteins, "+", paste(CovariatesE6, collapse = "+")))
  model <- glm(formula, family = poisson(link = "log"), data=MBFinal)
  robust_se<-sqrt(diag(vcovHC(model, type = "HC0")))
  RR<-exp(coef(model)[2])
  CI<-confint(model, level=0.95)[2,] 
  CI_adj <- exp(coef(model)[2] + qnorm(c(0.025, 0.975)) * robust_se[2])
  z_value<-coef(model)[2]/robust_se[2]
  p_value<-2*pnorm(abs(z_value), lower.tail = FALSE)
  results<-rbind(results, data.frame(Protein = proteins, RR = RR, CI_lower = CI[1], CI_upper = CI[2], P_Value = p_value))
}
results$P_Value_Adjusted<-p.adjust(results$P_Value, method = "BH")
significant_results<-results%>%filter(P_Value_Adjusted<0.05)
final_results<-merge(significant_results, MESAOlink3K_ProteinKeys, by.x = "Protein", by.y = "OlinkID")
final_table<-final_results%>%dplyr::select(ID = Assay, RR, P_Value_Adjusted)%>%arrange(P_Value_Adjusted); print(final_table)
ft<-flextable(final_table); ft<-set_table_properties(ft, width = 1, layout = "autofit"); ft<-padding(ft, padding.top = 0, padding.bottom = 0)
doc<-read_docx("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/results draft.docx"); doc<-body_add_flextable(doc, value=ft); print(doc, target = "/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/results draft.docx")

###########################################################

#Adjusted for Age, eGFR, SBP, HTN med
CovariatesE6<-c("age6c", "cepgfr6c", "sbp6c", "htnmed6c")
results<-data.frame(Protein = character(), RR = numeric(), CI_lower = numeric(), CI_upper = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)
for(proteins in colnames(MBFinal)[4:2944]) {
  formula <- as.formula(paste("mb_present ~", proteins, "+", paste(CovariatesE6, collapse = "+")))
  model <- glm(formula, family = poisson(link = "log"), data=MBFinal)
  robust_se<-sqrt(diag(vcovHC(model, type = "HC0")))
  RR<-exp(coef(model)[2])
  CI<-confint(model, level=0.95)[2,] 
  CI_adj<-exp(coef(model)[2] + qnorm(c(0.025, 0.975))*robust_se[2])
  z_value<-coef(model)[2]/robust_se[2]
  p_value<-2*pnorm(abs(z_value), lower.tail = FALSE)
  results<-rbind(results, data.frame(Protein = proteins, RR = RR, CI_lower = CI[1], CI_upper = CI[2], P_Value = p_value))
}
results$P_Value_Adjusted<-p.adjust(results$P_Value, method = "BH")
significant_results<-results%>%filter(P_Value_Adjusted<0.05)
final_results<-merge(significant_results, MESAOlink3K_ProteinKeys, by.x = "Protein", by.y = "OlinkID")
final_table<-final_results%>%dplyr::select(ID = Assay, RR, P_Value_Adjusted)%>%arrange(P_Value_Adjusted); print(final_table)
ft<-flextable(final_table); ft<-set_table_properties(ft, width = 1, layout = "autofit"); ft<-padding(ft, padding.top = 0, padding.bottom = 0)
doc<-read_docx("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/results draft.docx"); doc<-body_add_flextable(doc, value=ft); print(doc, target = "/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/results draft.docx")
