library('haven'); library('dplyr'); library('purrr'); library('car'); library('tidyverse'); library('broom'); library('officer'); library('flextable'); 
library('kableExtra'); library('knitr'); library('corrr'); library('corrplot'); library('ggcorrplot'); library('GGally'); 
library('glmnet'); library('caret'); library('tibble')
library(EnhancedVolcano); library('OlinkAnalyze'); library('ggrepel'); library('grid')

#Olink3K Data
MESAOlink3K<-read_csv("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/Proteomics_Olink/Intensity_Data/SMP_IntensityNormalized_03062024.csv")
MESAOlink3K_ProteinKeys<-read_csv("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/Proteomics_Olink/Mapping_Olink3K_OlinkID_Key_03062024.csv")
MESAOlink3K_LinkingData<-read_csv("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/Proteomics_Olink/Mapping_SMP_Plate_20240514.csv")
Olink3KFinal<-merge.data.frame(MESAOlink3K,MESAOlink3K_LinkingData,by="SampleID") #MESAOlink3K data merged by SampleID
Olink3Ke6<-Olink3KFinal[Olink3KFinal$exam==6,] #Olink3K data for Exam 6
#Rename "sidno" column to "subject_id" for merging with covariate data
Olink3Ke6<-Olink3Ke6%>%rename(subject_id=sidno)

#Covariate Data
CoV_e6<-read.delim("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/Data from MESA CC/SHARe_Exam6Main/SHARe_Exam6Main.txt")
CoV_e1<-read.delim("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/Data from MESA CC/SHARe_MesaExam1Main/SHARe_MesaExam1Main_DS.txt")
Edu_e1<-CoV_e1%>%select(subject_id, educ1) #Education data for Exam 1
Edu_e1<-Edu_e1%>%
  mutate(Edu = case_when(
    educ1 %in% 0:2 ~ 0,
    educ1 %in% 3:8 ~ 1
  )) # Education data recoded to 0=less than high school education, 1= completed high school
Edu_e1<-Edu_e1[!is.na(Edu_e1$Edu), ] #Remove NA education values
CoV_e6<-CoV_e6[!is.na(CoV_e6$age6c), ] #Remove NA age values
CoV_e6<-CoV_e6[!is.na(CoV_e6$site6c), ] #Remove NA site values
CoV_e6<-CoV_e6[!is.na(CoV_e6$cepgfr6c), ] #Remove NA eGFR values
CoV_e6<-CoV_e6[!is.na(CoV_e6$bmi6c), ] #Remove NA BMI values
CoV_e6<-CoV_e6[!is.na(CoV_e6$sbp6c), ] #Remove NA SBP values
CoV_e6<-CoV_e6[!is.na(CoV_e6$htnmed6c), ] #Remove NA HTN med values
CoV_e6<-CoV_e6[!is.na(CoV_e6$cig6c), ] #Remove NA cig use values
CoV_e6<-CoV_e6[!is.na(CoV_e6$ldl6), ] #Remove NA LDLC values
CoV_e6<-CoV_e6[!is.na(CoV_e6$afib6), ] #Remove NA AF Prevalent values
CoV_e6<-CoV_e6[!is.na(CoV_e6$dm036t), ] #Remove NA DM Prevalent values
PrevalentCVD<-read.delim("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/Data from MESA CC/SHARe_EventsThruYear2019/SHARe_EventsThruYear2019_DS.txt")
PrevalentCVD<-PrevalentCVD[!is.na(PrevalentCVD$mi), ] #Remove NA MI Prevalent values
PrevalentCVD<-PrevalentCVD[!is.na(PrevalentCVD$chf), ] #Remove NA CHF Prevalent values

#MRI Outcome Data
WMFA<-read.delim("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/Data from MESA CC/SHARe_AncilMesaAF_BMRIFAMuse/SHARe_AncilMesaAF_BMRIFAMuse_DS.txt")

#Merged Dataframes for Analyses
Prot_CoV<-merge.data.frame(Olink3Ke6,CoV_e6,by="subject_id") #Merge Olink3K and Covariate Data
Prot_CoV_Edu<-merge.data.frame(Prot_CoV,Edu_e1,by="subject_id") #Merge Olink3K, Covariate, and Education Data
Prot_CoV_Edu_CVD<-merge.data.frame(Prot_CoV_Edu,PrevalentCVD,by="subject_id") #Merge Olink3K, Covariate, Education and Prevalent CVD Data
WMHFAFinala<-merge.data.frame(Prot_CoV_Edu_CVD,WMFA,by="subject_id") #Merge Olink3K, Covariate, Education, Prevalent CVD and WM FA Data
WMFAFinal<-WMHFAFinala[is.na(WMHFAFinala$qc_code),] #QC participants removed (no missing FA data)

### WMFA Analyses ###

#WMFA Analyses - Baseline characteristics
mean(WMFAFinal$age6c); sd(WMFAFinal$age6c); summary(is.na(WMFAFinal$age6c)) # Age Mean and SD
table(WMFAFinal$gender1); summary(is.na(WMFAFinal$gender1)) # 0=female, 1=male
table(WMFAFinal$race1c); table(is.na(WMFAFinal$race1c)) # 1=white, 2=chinese, 3=black, #4=hispanic/latino
table(WMFAFinal$site6c); summary(is.na(WMFAFinal$site6c)) # 3=WFU, 4=COL, 5=JHU, 6=UMN, 7=NWU, 8=UCLA
table(WMFAFinal$Edu); summary(is.na(WMFAFinal$Edu)) # 0=less than high school education, 1= completed high school
mean(WMFAFinal$cepgfr6c); sd(WMFAFinal$cepgfr6c); summary(is.na(WMFAFinal$cepgfr6c)) # CKD-EPI eGFR Mean and SD (mL/min/1.73m2)
mean(WMFAFinal$bmi6c); sd(WMFAFinal$bmi6c); summary(is.na(WMFAFinal$bmi6c)) # BMI Mean and SD (kg/m2)
mean(WMFAFinal$sbp6c); sd(WMFAFinal$sbp6c); summary(is.na(WMFAFinal$sbp6c)) # SBP Mean and SD (mmHg)
table(WMFAFinal$htnmed6c); summary(is.na(WMFAFinal$htnmed6c)) # 1=Yes, 0=No
table(WMFAFinal$cig6c); summary(is.na(WMFAFinal$cig6c)) # 2=Current, 1=Former, 0=Never
mean(WMFAFinal$ldl6); sd(WMFAFinal$ldl6); summary(is.na(WMFAFinal$ldl6)) # LDL Mean and SD (mg/dL)
table(WMFAFinal$afib6); summary(is.na(WMFAFinal$afib6)); WMFAFinal<-WMFAFinal %>% mutate(AFprevalent = case_when( 
  afib6 == 0 ~ 0,
  afib6 == 1 ~ 1,
  afib6 == 9 ~ 0,
)); table(WMFAFinal$AFprevalent); summary(is.na(WMFAFinal$AFprevalent)) # 1=Prevalent AF, 0=No
table(WMFAFinal$dm036t); summary(is.na(WMFAFinal$dm036t)) # 1=Prevalent DM, 0=No
table(WMFAFinal$mi); summary(is.na(WMFAFinal$mi)) # 1=Prevalent MI, 0=No
table(WMFAFinal$chf); summary(is.na(WMFAFinal$chf)) # 1=Prevalent CHF, 0=No

#WMFA Analyses - WMFA
mean(WMFAFinal$wm); sd(WMFAFinal$wm); sum(is.na(WMFAFinal$wm))
summary(WMFAFinal$wm); hist(WMFAFinal$wm) # WMFA distribution

#Fully Adjusted LR Model (FDR)
CovariatesE6<-c("age6c", "gender1", "race1c", "site6c", "Edu", "cepgfr6c", "bmi6c", "sbp6c", "htnmed6c", "cig6c", "ldl6", "AFprevalent", "dm036t", "mi", "chf")
results <- data.frame(Protein = character(), Beta = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)
for(proteins in colnames(WMFAFinal)[3:2943]) {
  formula <- as.formula(paste("(wm) ~", proteins, "+", paste(CovariatesE6, collapse = "+")))
  model <- lm(formula, data = WMFAFinal)
  tidy_model<-tidy(model)
  results<-rbind(results, data.frame(Protein = proteins, Beta = tidy_model$estimate[2], P_Value = tidy_model$p.value[2]))}
results$P_Value_Adjusted<-p.adjust(results$P_Value, method = "BH")
significant_results<-results%>%filter(P_Value_Adjusted<0.05)
results<-merge(results, MESAOlink3K_ProteinKeys, by.x = "Protein", by.y = "OlinkID")
final_results<-merge(significant_results, MESAOlink3K_ProteinKeys, by.x = "Protein", by.y = "OlinkID")
final_table<-final_results%>%select(ID = Assay, Beta, P_Value_Adjusted)%>%arrange(P_Value_Adjusted); print(final_table)
ft<-flextable(final_table); ft<-set_table_properties(ft, width = 1, layout = "autofit"); ft<-padding(ft, padding.top = 0, padding.bottom = 0)
FT<-final_table %>% filter(Beta>0.006 | Beta<(-0.006))
doc<-read_docx("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/results draft.docx"); doc<-body_add_flextable(doc, value=ft); print(doc, target = "/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/results draft.docx")
duplicates<-final_table%>%group_by(ID)%>%filter(n()>1); print(duplicates); duplicate_count<-duplicates%>%count(ID); print(duplicate_count)

##############
# Table w/ all proteins significantly associated with WMFA
# First, let's create our filtered and formatted table data
supp_table1 <- results %>%
  # Select and rename columns for the table
  select(
    Protein = Assay,
    `Beta-coefficient` = Beta,
    `P-value` = P_Value,
    `Adjusted P-value` = P_Value_Adjusted
  ) %>%
  # Filter for significant results first (adjusted p-value < 0.05)
  filter(`Adjusted P-value` < 0.05) %>%
  # Remove duplicates, keeping the entry with the lowest adjusted p-value for each protein
  group_by(Protein) %>%
  slice_min(order_by = `Adjusted P-value`, n = 1) %>%
  ungroup() %>%
  # Format numeric columns with appropriate precision
  mutate(
    `Beta-coefficient` = format(round(`Beta-coefficient`, 5), nsmall = 5, scientific = FALSE),
    `P-value` = format.pval(`P-value`, digits = 3),
    `Adjusted P-value` = format.pval(`Adjusted P-value`, digits = 3)
  ) %>%
  # Sort by adjusted p-value and then by p-value for ties
  arrange(`Adjusted P-value`, `P-value`)

# Let's verify our filtering worked correctly
print(paste("Number of significant unique proteins:", nrow(supp_table1)))

# Create a flextable object with the specified formatting
ft_supp1 <- flextable(supp_table1) %>%
  # Set theme for consistent formatting
  theme_box() %>%
  # Set font to Arial 12
  font(fontname = "Arial", part = "all") %>%
  fontsize(size = 12, part = "all") %>%
  # Align columns appropriately - left align protein names, right align numeric values
  align(align = "left", part = "all") %>%
  align(j = 2:4, align = "right", part = "body") %>%
  # Add borders
  border_outer(border = fp_border(width = 1)) %>%
  border_inner_h(border = fp_border(width = 0.5)) %>%
  # Set column headers to bold
  bold(part = "header") %>%
  # Set line spacing for double spacing
  height_all(height = 24) %>%
  # Add padding for better readability
  padding(padding.top = 6, padding.bottom = 6, padding.left = 4, padding.right = 4) %>%
  # Set table properties
  set_table_properties(layout = "autofit", width = 1)

# Create a new Word document and add content
doc <- read_docx() %>%
  # Add title
  body_add_par(value = "Supplementary Table 1. Plasma Protein Associations with White Matter Fractional Anisotropy.", 
               style = "Normal") %>%
  # Add spacing after title
  body_add_par(value = "", style = "Normal") %>%
  # Add the table
  body_add_flextable(value = ft_supp1) %>%
  # Add spacing before footnote
  body_add_par(value = "", style = "Normal") %>%
  # Add footnote
  body_add_par(value = "Table shows proteins with statistically significant associations (adjusted p-value < 0.05). Beta-coefficients and p-values are from multivariable linear regression models adjusted for age, sex, race/ethnicity, education, estimated glomerular filtration rate, body-mass index, systolic blood pressure, antihypertensive drug therapy, tobacco use, low-density lipoprotein cholesterol, and prevalent atrial fibrillation, diabetes, myocardial infarction, and heart failure.", 
               style = "Normal")

# Save the document
print(doc, target = "/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/Supplementary_Table_1.docx")
########################################

# VP
buffer_fc <- 0.001
min_fc <- min(results$Beta) - buffer_fc
max_fc <- max(results$Beta) + buffer_fc
min_p <- -log10(max(results$P_Value_Adjusted))
max_p <- -log10(min(results$P_Value_Adjusted, na.rm = TRUE)) + 0.3
EnhancedVolcano(results,
                lab = results$Assay,
                x = 'Beta',
                y = 'P_Value_Adjusted',
                xlab = 'Beta-coefficient',
                title = 'Volcano Plot',
                subtitle = '',
                legendPosition = 'none',
                xlim = c(min_fc, max_fc),
                ylim = c(0, max_p),
                pCutoff = 0.05,
                FCcutoff = 0.006,
                pointSize = 1.5,
                labSize = 3,
                col = c('grey30', 'grey30', 'blue', 'blue'),
                selectLab = FT$ID,
                drawConnectors = TRUE,
                widthConnectors = 1,
                colConnectors = 'black',
                boxedLabels = TRUE,
                max.overlaps = Inf)

#LASSO regression
X<-as.matrix(WMFAFinal[,significant_results$Protein]) # protein matrix
Y<-log(WMFAFinal$wm) # log-transformed WMFA
set.seed(123) # set seed for reproducibility
trainIndex<-createDataPartition(Y, p=0.5, list=FALSE) # 50/50 partition of training and test data
X_train<-X[trainIndex,]; X_test<-X[-trainIndex,]; Y_train<-Y[trainIndex]; Y_test<-Y[-trainIndex]
set.seed(123) # set seed
cv.lasso<-cv.glmnet(X_train, Y_train, alpha = 1, family = "gaussian", type.measure = "mse", nfolds = 10) # LASSO regression
best_lambda<-cv.lasso$lambda.min # best lambda
print(paste("Best lambda:", best_lambda))
final_model<-glmnet(X_train, Y_train, alpha = 1, lambda = best_lambda) # final model
lasso_coef<-coef(final_model) # coefficients
nonzero_coef<-lasso_coef[lasso_coef[,1]!=0,] ; print(nonzero_coef) # non-zero coefficients
nonzero_coef <- round(nonzero_coef, 5) # round coefficients to 4 decimal places
predictions<-predict(final_model, s = best_lambda, newx = as.matrix(X_test)) # predictions
test_mse<-mean((predictions-Y_test)^2); print(test_mse) # test MSE
plot(cv.lasso) # plot cross-validation results
nonzero_coef<-rownames_to_column(as.data.frame(nonzero_coef), var = "Protein")
final_lasso<-merge(nonzero_coef, MESAOlink3K_ProteinKeys, by.x = "Protein", by.y = "OlinkID", all.x = TRUE)
final_lasso<-select(final_lasso, -Protein, -UniProt, -Panel); view(final_lasso)
fl<-flextable(final_lasso); fl<-set_table_properties(fl, width = 1, layout = "autofit"); fl<-padding(fl, padding.top = 0, padding.bottom = 0)
fl<-font(fl, fontname = "Arial", part = "all"); fl<-fontsize(fl, size = 12, part = "all"); fl<-line_spacing(fl, space = 2)
doc<-read_docx("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/results draft.docx"); doc<-body_add_flextable(doc, value=fl); print(doc, target = "/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/results draft.docx")

#################################################################################

#LR Model - unadjusted
results<-data.frame(Protein = character(), BCoefficient = numeric(), AdjPValue = numeric(), stringsAsFactors = FALSE)
for(proteins in colnames(WMFAFinal)[3:2943]) {
  formula <- as.formula(paste("(wm) ~", proteins))
  model <- lm(formula, data = WMFAFinal)
  tidy_model<-tidy(model)
  results<-rbind(results, data.frame(Protein = proteins, Beta = tidy_model$estimate[2], P_Value = tidy_model$p.value[2]))
}
results$P_Value_Adjusted<-p.adjust(results$P_Value, method = "BH")
significant_results<-results%>%filter(P_Value_Adjusted<0.05)
final_results<-merge(significant_results, MESAOlink3K_ProteinKeys, by.x = "Protein", by.y = "OlinkID")
final_table<-final_results%>%select(ID = Assay, Beta, P_Value_Adjusted)%>%arrange(P_Value_Adjusted); print(final_table)
ft<-flextable(final_table); ft<-set_table_properties(ft, width = 1, layout = "autofit"); ft<-padding(ft, padding.top = 0, padding.bottom = 0)
doc<-read_docx("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/results draft.docx"); doc<-body_add_flextable(doc, value=ft); print(doc, target = "/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/results draft.docx")

#Model Diagnostics
for (proteins in significant_results$Protein) {
  model<-lm(reformulate(proteins, response = "wm"), data = WMFAFinal)
  plot(model, which = 1)
  plot(model, which = 2)
  plot(model, which = 3)
  plot(model, which = 5)
}

#Fully Adjusted LR Model with bonferroni correction
CovariatesE6<-c("age6c", "race1c", "site6c", "Edu", "cepgfr6c", "bmi6c", "sbp6c", "htnmed6c", "cig6c", "ldl6", "AFprevalent", "dm036t", "mi", "chf")
results<-data.frame(Protein = character(), BCoefficient = numeric(), AdjPValue = numeric(), stringsAsFactors = FALSE)
for(proteins in colnames(WMFAFinal)[4:2944]) {
  formula <- as.formula(paste("(wm) ~", proteins, "+", paste(CovariatesE6, collapse = "+")))
  model <- lm(formula, data = WMFAFinal)
  tidy_model<-tidy(model)
  results<-rbind(results, data.frame(Protein = proteins, Beta = tidy_model$estimate[2], P_Value = tidy_model$p.value[2]))}
results$P_Value_Adjusted<-p.adjust(results$P_Value, method = "bonferroni")
significant_results<-results%>%filter(P_Value_Adjusted<0.05)
final_results<-merge(significant_results, MESAOlink3K_ProteinKeys, by.x = "Protein", by.y = "OlinkID")
final_table<-final_results%>%select(ID = Assay, Beta, P_Value_Adjusted)%>%arrange(P_Value_Adjusted); print(final_table)
ft<-flextable(final_table); ft<-set_table_properties(ft, width = 1, layout = "autofit"); ft<-padding(ft, padding.top = 0, padding.bottom = 0)
doc<-read_docx("/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/results draft.docx"); doc<-body_add_flextable(doc, value=ft); print(doc, target = "/Users/rizwankalani/Library/CloudStorage/OneDrive-UW/MESA Olink Proteomics and cSVD-Cog/results draft.docx")