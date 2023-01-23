# Ýlayda Yýlmaz
# Dilay Gümüþ

# Libraries
library("dplyr")
library("ggpubr")
library(MVN)
library(AID)
library(robustbase)
library(tidyverse)
library(ggplot2)
library(MASS)
library(tigerstats)
library(e1071)
library(mvnormalTest)
library(psych)
library(readr)
library(heplots)
library("ICSNP")
library(rstatix)
library(gridExtra)
library(mice)
library(ggplot2)
library(klaR)
library(psych)
library(factoextra)

# DATA
# all data
patientDataset <- read_csv("archive.zip")

################################################################################################################################
# Ýlayda Yýlmaz

##### Cleaning 

# select columns - 91713 obs
selected_data=patientDataset[,c(1,2,3,4,5,8,13,24,27,30,34,35,55,76,85)]

## complete NA's - 90973
categoric_omitted_data=selected_data %>% drop_na(c(encounter_id,patient_id,hospital_id,gender,icu_type,diabetes_mellitus,hospital_death))
imputed_Data <- mice(categoric_omitted_data, m=5, maxit = 10, method = 'pmm', seed = 500)
completeData <- complete(imputed_Data,2)
# this process takes time.

################################################################################################################################
# Dilay Gümüþ

#### Visualizing

my_hist=hist(completeData$age , breaks=100  , plot=F)
# Color vector
my_color= ifelse(my_hist$breaks < 40, rgb(0.2,0.8,0.5,0.5) , ifelse (my_hist$breaks >= 75, "purple", rgb(0.2,0.2,0.2,0.2) ))
# Final plot
plot(my_hist, col=my_color , border=F , main="" , xlab="value of the variable", xlim=c(0,100) )

#sample
ggplot(completeData, aes(x=gender, y=heart_rate_apache)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) + 
  xlab("gender")+
  geom_boxplot(fill="orange")


ggplot(completeData, aes(x=as.factor(gender), y=resprate_apache)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) + 
  xlab("gender")+
  geom_boxplot(fill="yellow")

################################################################################################################################
# Ýlayda Yýlmaz

#### Factor Analysis

selected_data_fa=patientDataset[,c("encounter_id","patient_id","hospital_id","gender",
                                   "icu_type","diabetes_mellitus","hospital_death",
                                   "heart_rate_apache", "resprate_apache",
                                   "d1_diasbp_min","d1_mbp_min",
                                   "d1_spo2_min","d1_sysbp_min",
                                   "d1_temp_min","h1_diasbp_min" ,
                                   "h1_mbp_min","h1_spo2_min",
                                   "h1_sysbp_min","d1_glucose_min","d1_potassium_min")]

categoric_omitted_data_fa=selected_data_fa %>% drop_na(c(encounter_id,patient_id,hospital_id,gender,icu_type,diabetes_mellitus,hospital_death))
imputed_Data_fa <- mice(categoric_omitted_data_fa, m=5, maxit = 10, method = 'pmm', seed = 500)
completeData_fa <- complete(imputed_Data_fa,2)
faData<-completeData_fa[,-c(1,2,3,4,5,6,7,8,9)]

cm <- cor(faData, method="pearson")
corrplot::corrplot(cm, method= "number", order = "hclust")
#we can observe that there are some correlated variables
KMO(r=cm)
#Since MSA > 0.5, we can run Factor Analysis on this data. Besides, Bartletts test of sphericity should be significant.
print(cortest.bartlett(cm,nrow(faData)))
#The Kaiser-Meyer Olkin (KMO) and Bartletts Test measure of sampling adequacy were used to examine the appropriateness of Factor Analysis. 
#The approximate of Chi-square is 1078296 with 36 degrees of freedom, which is significant at 0.05 Level of significance. 
#The KMO statistic of 0.77 is also large (greater than 0.50). Hence Factor Analysis is considered as an appropriate technique for further analysis of the data.
parallel <- fa.parallel(faData, fm = "minres", fa = "fa")
parallel

# Lets see whether 5 factor is enough to group the variables. 
factanal(faData, factors = 5, method ="mle", lower = 0.01,nstart=5)$PVAL
#Since p value is less than alpha, we reject H0. 6 factor solution is not adequate.
factanal(faData, factors = 6, method ="mle", lower = 0.01,nstart=5)$PVAL
factanal(faData, factors = 7, method ="mle", lower = 0.01,nstart=11)$PVAL
# We tried with more factors but still they are not adequate. So we used all variables.

################################################################################################################################
# Ýlayda Yýlmaz

####  Multivariate Multiple Linear Regression:

mlm1 <- lm(cbind(heart_rate_apache, resprate_apache) ~ d1_diasbp_min + d1_mbp_min + d1_spo2_min  +
             d1_sysbp_min + d1_spo2_min+d1_sysbp_min + d1_temp_min + h1_diasbp_min + h1_mbp_min + h1_spo2_min +
             h1_sysbp_min + d1_glucose_min + d1_potassium_min,
           data = completeData_fa)
summary(mlm1)


################################################################################################################################
# Ýlayda Yýlmaz


#### MANOVA :

# icu_admit_source : The location of the patient prior to being admitted to the unit
# heart_rate_apache   : The heart rate measured during the first 24 hours which results in the highest APACHE III score
# resprate_apache : The respiratory rate measured during the first 24 hours which results in the highest APACHE III score

# In order to run some tests, we created a small sample that represents our data well.
set.seed(2049)
sample200=popsamp(200,completeData)
sample200$icu_type<-as.factor(sample200$icu_type)
# the descriptive statistics
sample200 %>% group_by(icu_type) %>%  summarise(n = n(), mean = mean(resprate_apache,na.rm=T ), sd = sd(resprate_apache,na.rm=T))
# visualize the data
p1 <- ggplot(sample200, aes(x = icu_type, y = heart_rate_apache, fill = icu_type)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) + theme(legend.position="top")+
  labs(title = "The Box Plot of heart_rate_apache by icu_admit_source" )
p2 <- ggplot(sample200, aes(x = icu_type, y = resprate_apache, fill = icu_type)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) + theme(legend.position="top")+
  labs(title = "The Box Plot of resprate_apache by icu_admit_source" )
grid.arrange(p1, p2, ncol=2)

### Assumptions
## For normality :
mvn(sample200[, c("heart_rate_apache", "resprate_apache")],mvnTest = "mardia")
y<-sample200%>%select(heart_rate_apache,resprate_apache)
error.bars (y, ylab="Group Means", xlab=" Dependent Variables")

## Homogeneity of the variance-covariance matrices:
# We will use Box’s M test to assess the homogeneity of the variance-covariance matrices.
# Null hypothesis: variance-covariance matrices are equal for each combination formed by each group in the independent variable
boxM(Y = sample200[, c("heart_rate_apache", "resprate_apache")], group = factor(sample200$icu_type))
# As the p-value is non-significant (p > 0.05) for Box’s M test, 
# we fail to reject the null hypothesis,
# and conclude that variance-covariance matrices are equal for each combination of the dependent variable formed by each group in the independent variable.


### 
m1 <- manova(cbind(heart_rate_apache,resprate_apache) ~ icu_type, data = sample200)
summary(m1)
# We’ll interpret this table as usual ANOVA table. The value which should take into account is Pr(>F) value which indicates p value.
# Therefore, we are 95% confident that at least one drug type is significantly different than others since p<alpha.
# If we reject the null hypothesis in MANOVA, we need to hold a post-hoc analysis to see which one causes the difference. To do so,
summary.aov(m1)
# From the output above, it can be seen that the heart_rate_apache variable highly significantly different among icu type.
# However, we cannot make a similar conclusion for resprate_apache.



################################################################################################################################
# Ýlayda Yýlmaz



#### Hotelling's T^2 :

### Assumptions

## |cov|>0

## cov1=cov2=cov
boxM(Y = sample200[, c("d1_heartrate_min", "h1_heartrate_min")], group = factor(sample200$gender))
# As the p-value is non-significant (p = 0.6235 > 0.05) for Box’s M test, we fail to reject the null hypothesis and 
# conclude that variance-covariance matrices are equal for each combination of the dependent variable formed by each group in the independent variable.

## Both samples follow normal distributions.
table(sample200$gender)
sample200 %>% group_by(gender) %>%  shapiro_test(d1_heartrate_min,h1_heartrate_min)
# As the p-value is non-significant (p > 0.05) for each combination of independent and dependent variables, 
# we fail to reject the null hypothesis and conclude that data follows univariate normality.

HotellingsT2(cbind(sample200$d1_heartrate_min,sample200$h1_heartrate_min) ~ sample200$gender)

# Since p<alpha we fail to reject H0. Therefore, we don’t have enough evidence to prove that the mean of responses change with respect to gender.
y<-sample200%>%select(d1_heartrate_min,h1_heartrate_min)

error.bars (y, ylab="Group Means", xlab=" Dependent Variables")





################################################################################################################################
# Dilay Gümüþ

#### Principal Component Analysis

PCAData <- patientDataset[,c("d1_diasbp_min","d1_mbp_min","h1_sysbp_min",
                             "h1_mbp_min","intubated_apache","diabetes_mellitus")]

categoric_omitted_data_pca=PCAData %>% drop_na(c(diabetes_mellitus))
imputed_Data_pca <- mice(categoric_omitted_data_pca, m=5, maxit = 10, method = 'pmm', seed = 500)
completeData_pca <- complete(imputed_Data_pca,2)

str(completeData_pca)
summary(completeData_pca)

pca_new <- completeData_pca[,-c(5,6)]



res <- cor(pca_new, method="pearson")
corrplot::corrplot(res, method= "color", order = "hclust")

cor(pca_new)

pca_new1 <- pca_new[,-4]
pca1 <- prcomp(pca_new1)
summary(pca1)

pca1$rotation
pca1$x


pca1$sdev

fviz_eig(pca1,addlabels=TRUE) #represent the proportion values

pca<-pca1$x[,1:2]
head(pca)

res1 <- cor(pca, method="pearson")
corrplot::corrplot(res1, method= "color", order = "hclust")

cor(pca_new1,pca)

#Loading Plots
fviz_pca_var(pca1, col.var = "contrib")+ scale_color_gradient2( low="red", mid="green",
                                                                high="blue", midpoint=96, space = "Lab")

fviz_pca_var(pca1,axes = c(2, 3))

fviz_pca_var(pca1, select.var = list(contrib = 2))

fviz_contrib(pca1, choice = "ind", axes = 1:2) + coord_flip()

fviz_pca_ind(pca1, label="none", habillage = PCAData$intubated_apache,
             addEllipses=TRUE, ellipse.level=0.95)

################################################################################################################################
# Dilay Gümüþ


####  Discr & Class:
# diabetes_mellitus -- d1_heartrate_max - d1_heartrate_min


selected_data_discr=patientDataset[,c("encounter_id","patient_id","hospital_id","gender",
                                    "icu_type","diabetes_mellitus","hospital_death",
                                    "h1_diasbp_min","h1_heartrate_min",
                                    "h1_resprate_min","h1_spo2_min","h1_sysbp_min")]

categoric_omitted_data_try=selected_data_discr %>% drop_na(c(encounter_id,patient_id,hospital_id,gender,icu_type,diabetes_mellitus,hospital_death))
imputed_Data_try <- mice(categoric_omitted_data_try, m=5, maxit = 10, method = 'pmm', seed = 500)
completeData_try <- complete(imputed_Data_try,2)

smalldata<-completeData_try[,c("diabetes_mellitus",
                               "h1_diasbp_min","h1_heartrate_min",
                               "h1_resprate_min","h1_spo2_min","h1_sysbp_min")]

str(smalldata)
summary(smalldata)

set.seed(2049)
sample5000=popsamp(2000,smalldata)
sample <- sample(c(TRUE, FALSE), nrow(sample5000), replace=TRUE, prob=c(0.7,0.3))
train <- sample5000[sample, ]
test <- sample5000[!sample, ]

model <- lda(diabetes_mellitus~.,data = train)
model

plot(model)
#partition plot
train$diabetes_mellitus<-as.factor(train$diabetes_mellitus)

# Train performance
train_predict<- predict(model,train)$class
table_train <- table(Predicted =train_predict, Actual = train$diabetes_mellitus)
table_train


sum(diag(table_train))/sum(table_train)






