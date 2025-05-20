library('haven'); library('dplyr'); library('purrr'); library('car'); library('tidyverse'); library('broom'); library('officer'); library('flextable'); library('EnhancedVolcano')

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
Edu_e1<-CoV_e1%>%select(subject_id, educ1) #Education data for Exam 1
# Education data recoded to 0=less than high school education, 1= completed high school
Edu_e1 <- CoV_e1 %>%
  select(subject_id, educ1) %>%
  mutate(Edu = case_when(
    educ1 %in% 0:2 ~ 0,
    educ1 %in% 3:8 ~ 1
  )) %>%
  filter(!is.na(Edu)) 
#Remove NA values for covariates (SBP, HTN med use, cig use, LDLC, AF Prevalent)
CoV_e6 <- CoV_e6 %>%
  filter(!is.na(sbp6c), !is.na(htnmed6c), !is.na(cig6c), !is.na(ldl6), !is.na(afib6))
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
MRIVol<-read.delim("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/Data from MESA CC/SHARe_AncilMesaAF_BMRIROIVol/SHARe_AncilMesaAF_BMRIROIVol_DS.txt")
ePVS<-read.delim("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/Data from MESA CC/MESAe6as253_BMRIPVS_SIDNO_20231121.txt")
ePVS2<-ePVS%>%rename(subject_id=sidno) #Rename subject_id column

#Merged Dataframes for Analyses
Prot_CoV<-merge.data.frame(Olink3Ke6,CoV_e6,by="subject_id") #Merge Olink3K and Covariate Data
Prot_CoV_Edu<-merge.data.frame(Prot_CoV,Edu_e1,by="subject_id") #Merge Olink3K, Covariate, Education Data
Prot_CoV_Edu_CVD<-merge.data.frame(Prot_CoV_Edu,PrevalentCVD,by="subject_id") #Merge Olink3K, Covariate, Education and Prevalent CVD Data
Prot_CoV_Edu_CVD_MRIVol<-merge.data.frame(Prot_CoV_Edu_CVD,MRIVol,by="subject_id") #Merge Olink3K, Covariate, Education, and MRI Volume Data
ePVSFinala<-merge.data.frame(Prot_CoV_Edu_CVD_MRIVol,ePVS2,by="subject_id") #Merge Olink3K, Covariate, MRI Volume and WMH Volume Data
ePVSFinal<-ePVSFinala[ePVSFinala$pvs_exclude !=1,] #QC exclusions (pvs_exclude=1)

### ePVS Analyses ###

#ePVS Analyses - Baseline characteristics
mean(ePVSFinal$age6c); sd(ePVSFinal$age6c); summary(is.na(ePVSFinal$age6c)) # Age Mean and SD
table(ePVSFinal$gender1); summary(is.na(ePVSFinal$gender1)) # 0=female, 1=male
table(ePVSFinal$race1c); table(is.na(ePVSFinal$race1c)) # 1=white, 2=chinese, 3=black, #4=hispanic/latino
table(ePVSFinal$site6c); summary(is.na(ePVSFinal$site6c)) # 3=WFU, 4=COL, 5=JHU, 6=UMN, 7=NWU, 8=UCLA
table(ePVSFinal$Edu); summary(is.na(ePVSFinal$Edu)) # 0=less than high school education, 1= completed high school
mean(ePVSFinal$cepgfr6c); sd(ePVSFinal$cepgfr6c); summary(is.na(ePVSFinal$cepgfr6c)) # CKD-EPI eGFR Mean and SD (mL/min/1.73m2)
mean(ePVSFinal$bmi6c); sd(ePVSFinal$bmi6c); summary(is.na(ePVSFinal$bmi6c)) # BMI Mean and SD (kg/m2)
mean(ePVSFinal$sbp6c); sd(ePVSFinal$sbp6c); summary(is.na(ePVSFinal$sbp6c)) # SBP Mean and SD (mmHg)
table(ePVSFinal$htnmed6c); summary(is.na(ePVSFinal$htnmed6c)) # 1=Yes, 0=No)
table(ePVSFinal$cig6c); summary(is.na(ePVSFinal$cig6c)) # 2=Current, 1=Former, 0=Never
mean(ePVSFinal$ldl6); sd(ePVSFinal$ldl6); summary(is.na(ePVSFinal$ldl6)) # LDL Mean and SD (mg/dL)
table(ePVSFinal$E4); summary(is.na(ePVSFinal$E4)) # 1=Presence of 1 or more E4 Allele, 0=Absence of E4 Allele
table(ePVSFinal$afib6); summary(is.na(ePVSFinal$afib6)); ePVSFinal<-ePVSFinal %>% mutate(AFprevalent = case_when( 
  afib6 == 0 ~ 0,
  afib6 == 1 ~ 1,
  afib6 == 9 ~ 0,
)); table(ePVSFinal$AFprevalent); summary(is.na(ePVSFinal$AFprevalent)) # 1=Prevalent AF, 0=No
table(ePVSFinal$dm036t); summary(is.na(ePVSFinal$dm036t)) # 1=Prevalent DM, 0=No
table(ePVSFinal$mi); summary(is.na(ePVSFinal$mi)) # 1=Prevalent MI, 0=No
table(ePVSFinal$chf); summary(is.na(ePVSFinal$chf)) # 1=Prevalent CHF, 0=No
mean(ePVSFinal$icv); sd(ePVSFinal$icv) # Intracranial Volume Mean and SD (cm3)

#ePVS Analyses - ePVS Volume
mean(ePVSFinal$epvs_wholebrain_vol); sd(ePVSFinal$epvs_wholebrain_vol); summary(is.na(ePVSFinal$epvs_wholebrain_vol)) # ePVS Volume Mean and SD (cm3)
summary(ePVSFinal$epvs_wholebrain_vol); hist(ePVSFinal$epvs_wholebrain_vol) # WMH Volume (cm3) distribution
summary(log(ePVSFinal$epvs_wholebrain_vol)); hist(log(ePVSFinal$epvs_wholebrain_vol)) # log-transformmed WMH Volume (cm3) distribution

#LR Model - adjusting for ICV
CovariatesE6<-c("icv")
results<-data.frame(Protein = character(), BCoefficient = numeric(), AdjPValue = numeric(), stringsAsFactors = FALSE)
for(proteins in colnames(ePVSFinal)[3:2943]) {
  formula <- as.formula(paste("(log(epvs_wholebrain_vol)) ~", proteins, "+", paste(CovariatesE6, collapse = "+")))
  model <- lm(formula, data = ePVSFinal)
  tidy_model<-tidy(model)
  results<-rbind(results, data.frame(Protein = proteins, Beta = tidy_model$estimate[2], P_Value = tidy_model$p.value[2]))
}
results$P_Value_Adjusted<-p.adjust(results$P_Value, method = "BH")
significant_results<-results%>%filter(P_Value_Adjusted<0.05)
final_results<-merge(significant_results, MESAOlink3K_ProteinKeys, by.x = "Protein", by.y = "OlinkID")
final_table<-final_results%>%select(ID = Assay, Beta, P_Value_Adjusted)%>%arrange(P_Value_Adjusted); print(final_table)
ft<-flextable(final_table); ft<-set_table_properties(ft, width = 1, layout = "autofit"); ft<-padding(ft, padding.top = 0, padding.bottom = 0)
doc<-read_docx("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/results draft.docx"); doc<-body_add_flextable(doc, value=ft); print(doc, target = "/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/results draft.docx")

#Fully Adjusted LR Model
CovariatesE6<-c("icv", "age6c", "gender1", "race1c", "site6c", "Edu", "cepgfr6c", "bmi6c", "sbp6c", "htnmed6c", "cig6c", "ldl6", "AFprevalent", "dm036t", "mi", "chf")
results<-data.frame(Protein = character(), BCoefficient = numeric(), AdjPValue = numeric(), stringsAsFactors = FALSE)
for(proteins in colnames(ePVSFinal)[3:2943]) {
  formula <- as.formula(paste("(log(epvs_wholebrain_vol)) ~", proteins, "+", paste(CovariatesE6, collapse = "+")))
  model <- lm(formula, data = ePVSFinal)
  tidy_model<-tidy(model)
  results<-rbind(results, data.frame(Protein = proteins, Beta = tidy_model$estimate[2], P_Value = tidy_model$p.value[2]))}
results$P_Value_Adjusted<-p.adjust(results$P_Value, method = "BH")
significant_results<-results%>%filter(P_Value_Adjusted<0.05)
final_results<-merge(significant_results, MESAOlink3K_ProteinKeys, by.x = "Protein", by.y = "OlinkID")
final_table<-final_results%>%select(ID = Assay, Beta, P_Value_Adjusted)%>%arrange(P_Value_Adjusted); print(final_table)
ft<-flextable(final_table); ft<-set_table_properties(ft, width = 1, layout = "autofit"); ft<-padding(ft, padding.top = 0, padding.bottom = 0)
doc<-read_docx("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/results draft.docx"); doc<-body_add_flextable(doc, value=ft); print(doc, target = "/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/results draft.docx")

################################

#Model Diagnostics
for (proteins in significant_results$Protein) {
  model<-lm(reformulate(proteins, response = "log(wm.y)"), data = ePVSFinal)
  plot(model, which = 1)
  plot(model, which = 2)
  plot(model, which = 3)
  plot(model, which = 5)
}

#LR Model - unadjusted (not adjusted for ICV)
results<-data.frame(Protein = character(), BCoefficient = numeric(), AdjPValue = numeric(), stringsAsFactors = FALSE)
for(proteins in colnames(ePVSFinal)[4:2944]) {
  formula <- as.formula(paste("(log(epvs_wholebrain_vol)) ~", proteins))
  model <- lm(formula, data = ePVSFinal)
  tidy_model<-tidy(model)
  results<-rbind(results, data.frame(Protein = proteins, Beta = tidy_model$estimate[2], P_Value = tidy_model$p.value[2]))
}
results$P_Value_Adjusted<-p.adjust(results$P_Value, method = "BH")
significant_results<-results%>%filter(P_Value_Adjusted<0.05)
final_results<-merge(significant_results, MESAOlink3K_ProteinKeys, by.x = "Protein", by.y = "OlinkID")
final_table<-final_results%>%select(ID = Assay, Beta, P_Value_Adjusted)%>%arrange(P_Value_Adjusted); print(final_table)
ft<-flextable(final_table); ft<-set_table_properties(ft, width = 1, layout = "autofit"); ft<-padding(ft, padding.top = 0, padding.bottom = 0)
doc<-read_docx("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/results draft.docx"); doc<-body_add_flextable(doc, value=ft); print(doc, target = "/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/results draft.docx")