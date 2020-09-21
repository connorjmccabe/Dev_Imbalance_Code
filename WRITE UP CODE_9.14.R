rm(list=ls())
require(dplyr)
require(tidyr)
require(FactoMineR)
require(missMDA)
require(RMediation)
require(semTools)
require(lavaan)
options(scipen=20)
options(digits=4) #show 4 digits when displaying results (e.g., 0.000 or 00.00)
options(warn=-1) #supress warnings for readability

setwd("~/Dropbox/Postdoc/NCANDA Project/Analyses")
# 
#Minor data curation
ncanda_w1_w4<-read.csv("df_final_wide_mi.csv",header=T) %>% na_if(-9999) %>% rowwise() %>% 
  mutate(upps_ss_16r=mean(upps23_ss_16r,upps31_ss_16r,upps36_ss_16r,upps46_ss_16r,na.rm=T),
         upps_ss_17r=mean(upps23_ss_17r,upps31_ss_17r,upps36_ss_17r,upps46_ss_17r,na.rm=T),
         upps_ss_18r=mean(upps23_ss_18r,upps31_ss_18r,upps36_ss_18r,upps46_ss_18r,na.rm=T),
         upps_ss_19r=mean(upps23_ss_19r,upps31_ss_19r,upps36_ss_19r,upps46_ss_19r,na.rm=T),
         upps_ss_20r=mean(upps23_ss_20r,upps31_ss_20r,upps36_ss_20r,upps46_ss_20r,na.rm=T),
         upps_ss_21r=mean(upps23_ss_21r,upps31_ss_21r,upps36_ss_21r,upps46_ss_21r,na.rm=T),
         upps_pmt_16r=mean(upps5_pmt_16r,upps16_pmt_16r,upps28_pmt_16r,upps48_pmt_16r,na.rm=T),
         upps_pmt_17r=mean(upps5_pmt_17r,upps16_pmt_17r,upps28_pmt_17r,upps48_pmt_17r,na.rm=T),
         upps_pmt_18r=mean(upps5_pmt_18r,upps16_pmt_18r,upps28_pmt_18r,upps48_pmt_18r,na.rm=T),
         upps_pmt_19r=mean(upps5_pmt_19r,upps16_pmt_19r,upps28_pmt_19r,upps48_pmt_19r,na.rm=T),
         upps_pmt_20r=mean(upps5_pmt_20r,upps16_pmt_20r,upps28_pmt_20r,upps48_pmt_20r,na.rm=T),
         upps_pmt_21r=mean(upps5_pmt_21r,upps16_pmt_21r,upps28_pmt_21r,upps48_pmt_21r,na.rm=T))

#Computing alphas for each age
require(psych)
psych::alpha(ncanda_w1_w4%>%dplyr::select(upps23_ss_16,upps31_ss_16,upps36_ss_16,upps46_ss_16))
psych::alpha(ncanda_w1_w4%>%dplyr::select(upps23_ss_17,upps31_ss_17,upps36_ss_17,upps46_ss_17))
psych::alpha(ncanda_w1_w4%>%dplyr::select(upps23_ss_18,upps31_ss_18,upps36_ss_18,upps46_ss_18))
psych::alpha(ncanda_w1_w4%>%dplyr::select(upps23_ss_19,upps31_ss_19,upps36_ss_19,upps46_ss_19))
psych::alpha(ncanda_w1_w4%>%dplyr::select(upps23_ss_20,upps31_ss_20,upps36_ss_20,upps46_ss_20))
psych::alpha(ncanda_w1_w4%>%dplyr::select(upps23_ss_21,upps31_ss_21,upps36_ss_21,upps46_ss_21))
psych::alpha(ncanda_w1_w4%>%dplyr::select(upps16_pmt_16,upps5_pmt_16,upps28_pmt_16,upps48_pmt_16))
psych::alpha(ncanda_w1_w4%>%dplyr::select(upps16_pmt_17,upps5_pmt_17,upps28_pmt_17,upps48_pmt_17))
psych::alpha(ncanda_w1_w4%>%dplyr::select(upps16_pmt_18,upps5_pmt_18,upps28_pmt_18,upps48_pmt_18))
psych::alpha(ncanda_w1_w4%>%dplyr::select(upps16_pmt_19,upps5_pmt_19,upps28_pmt_19,upps48_pmt_19))
psych::alpha(ncanda_w1_w4%>%dplyr::select(upps16_pmt_20,upps5_pmt_20,upps28_pmt_20,upps48_pmt_20))
psych::alpha(ncanda_w1_w4%>%dplyr::select(upps16_pmt_21,upps5_pmt_21,upps28_pmt_21,upps48_pmt_21))

########SEMS#####

####MEASUREMENT INVARIANCE
#Filter data to include only age 16 to 20
df.long<-read.csv("df_final_long.csv",header=T) %>% filter(visit_age_integer>=16 & visit_age_integer<= 20)

#SENSATION SEEKING
sscfa<-"
ss =~ upps23_ss + upps31_ss + upps36_ss + upps46_ss
"
require(lavaan)

ssmi<-cfa(sscfa, data=df.long)

ssmi.config<-cfa(sscfa, data=df.long,  group="visit_age_integer", estimator="ML")

#Chi-square test comparing age-agnostic CFA to age-stratified CFA. No significant difference in fit.
pchisq(abs(fitmeasures(ssmi)[c("chisq")]-fitmeasures(ssmi.config)[c("chisq")]),abs(fitmeasures(ssmi)[c("df")]-fitmeasures(ssmi.config)[c("df")]))

#Metric and scalar models for measurement invaraince testing
ssmi.metric<-cfa(sscfa, data=df.long,  group="visit_age_integer", group.equal="loadings")
ssmi.scalar<-cfa(sscfa, data=df.long,  group="visit_age_integer", group.equal=c("loadings","intercepts"))

#Measurement invariance test
anova(ssmi.config, ssmi.metric,ssmi.scalar)

#PREMEDITATION
pmtcfa<-"
pmt =~ upps16_pmt + upps5_pmt + upps28_pmt + upps48_pmt
"
pmtmi<-cfa(pmtcfa, data=df.long)
pmtmi.config<-cfa(pmtcfa, data=df.long,  group="visit_age_integer")

pchisq(abs(fitmeasures(pmtmi)[c("chisq")]-fitmeasures(pmtmi.config)[c("chisq")]),abs(fitmeasures(pmtmi)[c("df")]-fitmeasures(pmtmi.config)[c("df")]))

pmtmi.metric<-cfa(pmtcfa, data=df.long,  group="visit_age_integer", group.equal="loadings")
pmtmi.scalar<-cfa(pmtcfa, data=df.long,  group="visit_age_integer", group.equal=c("loadings","intercepts"))
#Measurement invariance test
anova(pmtmi.config, pmtmi.metric,pmtmi.scalar)

########
#Unconditional Growth Models
########

m.uncon<-"

ss.l=~ 0*upps_ss_16r + 1*upps_ss_17r + 2*upps_ss_18r + 3*upps_ss_19r + 4*upps_ss_20r
#ss.q=~ 0*upps_ss.s_16 + 1*upps_ss.s_17 + 4*upps_ss.s_18 + 9*upps_ss.s_19 + 16*upps_ss.s_20 + 5*upps_ss.s_21
ss.i=~ 1*upps_ss_16r + 1*upps_ss_17r + 1*upps_ss_18r + 1*upps_ss_19r + 1*upps_ss_20r

pmt.l=~ 0*upps_pmt_16r + 1*upps_pmt_17r + 2*upps_pmt_18r + 3*upps_pmt_19r + 4*upps_pmt_20r
#pmt.q=~ 0*upps_pmt_16r + 1*upps_pmt_17r + 4*upps_pmt_18r + 9*upps_pmt_19r + 16*upps_pmt_20r
pmt.i=~ 1*upps_pmt_16r + 1*upps_pmt_17r + 1*upps_pmt_18r + 1*upps_pmt_19r + 1*upps_pmt_20r

pmt.i~~pmt.l
ss.i~~ss.l
ss.i~~pmt.i
pmt.l~~ss.i

"

#Model using the full sample (i.e. not stratified by sex)
uncon<-growth(model=m.uncon,data=ncanda_w1_w4,mimic="MPlus")
fitmeasures(uncon)[c("df","chisq","rmsea","cfi","tli","srmr")]

#Model stratefied by sex
uncon.sex<-growth(model=m.uncon,data=ncanda_w1_w4,mimic="MPlus",group="male")
fitmeasures(uncon.sex)[c("df","chisq","rmsea","cfi","tli","srmr")]

#SEX MODERATION

#Fixing Age 16 premeditation intercept to equivalence
uncon.sex.pmti<-growth(model=paste0(m.uncon,"pmt.i~c(g1,g1)*1"),data=ncanda_w1_w4,mimic="MPlus",group="male")
#Fixing premeditation growth to equivalence
uncon.sex.pmts<-growth(model=paste0(m.uncon,"pmt.l~c(g1,g1)*1"),data=ncanda_w1_w4,mimic="MPlus",group="male")
#Fixing Age 16 sensation seeking intercept to equivalence
uncon.sex.ssi<-growth(model=paste0(m.uncon,"ss.i~c(g1,g1)*1"),data=ncanda_w1_w4,mimic="MPlus",group="male")
#Fixing sensation seeking growth to equivalence
uncon.sex.sss<-growth(model=paste0(m.uncon,"ss.l~c(g1,g1)*1"),data=ncanda_w1_w4,mimic="MPlus",group="male")

#Chi-square difference tests assessing sex moderation of regression effects

lavTestLRT(uncon.sex,uncon.sex.pmti)
lavTestLRT(uncon.sex,uncon.sex.pmts)
lavTestLRT(uncon.sex,uncon.sex.ssi)
lavTestLRT(uncon.sex,uncon.sex.sss)

#View model results
View(standardizedsolution(uncon))

#########
#Parallel Process Growth Curve Model
#########

m.sspmthed<-'
ss.l=~ 0*upps_ss_16r + 1*upps_ss_17r + 2*upps_ss_18r + 3*upps_ss_19r + 4*upps_ss_20r
ss.i=~ 1*upps_ss_16r + 1*upps_ss_17r + 1*upps_ss_18r + 1*upps_ss_19r + 1*upps_ss_20r

pmt.l=~ 0*upps_pmt_16r + 1*upps_pmt_17r + 2*upps_pmt_18r + 3*upps_pmt_19r + 4*upps_pmt_20r
pmt.i=~ 1*upps_pmt_16r + 1*upps_pmt_17r + 1*upps_pmt_18r + 1*upps_pmt_19r + 1*upps_pmt_20r

hed.l=~ -4*alc_bingeyr_ord_16 + -3*alc_bingeyr_ord_17 + -2*alc_bingeyr_ord_18 + -1*alc_bingeyr_ord_19 + 0*alc_bingeyr_ord_20
hed.i=~ 1*alc_bingeyr_ord_16 + 1*alc_bingeyr_ord_17 + 1*alc_bingeyr_ord_18 + 1*alc_bingeyr_ord_19 + 1*alc_bingeyr_ord_20

hed.i~pmt.l + pmt.i + ss.l + ss.i + white + ses + male
hed.l~pmt.l + pmt.i + ss.l + ss.i + white + ses + male

hed.i~~hed.l

ss.i~~pmt.i
ss.i~~ss.l
pmt.i~~pmt.l
ss.l~~pmt.l

'
sspmthed<-growth(m.sspmthed,data=ncanda_w1_w4, estimator = "MLR", mimic="MPlus")
fitmeasures(sspmthed)[c("chisq","rmsea","cfi","tli", "srmr")]
View(standardizedsolution(sspmthed,ci=T))

sspmthed.sex<-growth(gsub("[+] male","",m.sspmthed), data=ncanda_w1_w4, check.gradient=F,estimator="MLR",mimic="MPlus",group="male")
fitmeasures(sspmthed.sex)[c("chisq","rmsea","cfi","tli", "srmr")]

#Moderation Tests
sspmthed.sex.hedi_ssi<-growth(gsub("hed.i~pmt.l [+] pmt.i [+] ss.l [+] ss.i [+] white [+] ses [+] male\nhed.l~pmt.l [+] pmt.i [+] ss.l [+] ss.i [+] white [+] ses [+] male",
                                     "hed.i~pmt.l + pmt.i + ss.l + c(g1,g1)*ss.i + white + ses\nhed.l~pmt.l + pmt.i + ss.l + ss.i + white + ses",m.sspmthed), data=ncanda_w1_w4, check.gradient=F,estimator = "MLR", mimic="MPlus", group="male")

sspmthed.sex.heds_ssi<-growth(gsub("hed.i~pmt.l [+] pmt.i [+] ss.l [+] ss.i [+] white [+] ses [+] male\nhed.l~pmt.l [+] pmt.i [+] ss.l [+] ss.i [+] white [+] ses [+] male",
                                     "hed.i~pmt.l + pmt.i + ss.l + ss.i + white + ses\n
                                     hed.l~pmt.l + pmt.i + ss.l + c(g1,g1)*ss.i + white + ses",m.sspmthed), data=ncanda_w1_w4, check.gradient=F,estimator="mlr",mimic="mplus",group="male")


sspmthed.sex.hedi_sss<-growth(gsub("hed.i~pmt.l [+] pmt.i [+] ss.l [+] ss.i [+] white [+] ses [+] male\nhed.l~pmt.l [+] pmt.i [+] ss.l [+] ss.i [+] white [+] ses [+] male",
                                     "hed.i~pmt.l + pmt.i + c(g1,g1)*ss.l + ss.i + white + ses\n
                                     hed.l~pmt.l + pmt.i + ss.l + ss.i + white + ses",m.sspmthed), data=ncanda_w1_w4, check.gradient=F,estimator="mlr",mimic="mplus",group="male")

sspmthed.sex.heds_sss<-growth(gsub("hed.i~pmt.l [+] pmt.i [+] ss.l [+] ss.i [+] white [+] ses [+] male\nhed.l~pmt.l [+] pmt.i [+] ss.l [+] ss.i [+] white [+] ses [+] male",
                                     "hed.i~pmt.l + pmt.i + ss.l + ss.i + white + ses\n
                                     hed.l~pmt.l + pmt.i + c(g1,g1)*ss.l + ss.i + white + ses",m.sspmthed), data=ncanda_w1_w4, check.gradient=F,estimator="mlr",mimic="mplus",group="male")

sspmthed.sex.hedi_prei<-growth(gsub("hed.i~pmt.l [+] pmt.i [+] ss.l [+] ss.i [+] white [+] ses [+] male\nhed.l~pmt.l [+] pmt.i [+] ss.l [+] ss.i [+] white [+] ses [+] male",
                                     "hed.i~pmt.l + c(g1,g1)*pmt.i + ss.l + ss.i + white + ses\n
                                     hed.l~pmt.l + pmt.i + ss.l + ss.i + white + ses",m.sspmthed), data=ncanda_w1_w4, check.gradient=F,estimator="mlr",mimic="mplus",group="male")

sspmthed.sex.heds_prei<-growth(gsub("hed.i~pmt.l [+] pmt.i [+] ss.l [+] ss.i [+] white [+] ses [+] male\nhed.l~pmt.l [+] pmt.i [+] ss.l [+] ss.i [+] white [+] ses [+] male",
                                     "hed.i~pmt.l + pmt.i + ss.l + ss.i + white + ses\n
                                     hed.l~pmt.l + c(g1,g1)*pmt.i + ss.l + ss.i + white + ses",m.sspmthed), data=ncanda_w1_w4, check.gradient=F,estimator="mlr",mimic="mplus",group="male")


sspmthed.sex.hedi_pres<-growth(gsub("hed.i~pmt.l [+] pmt.i [+] ss.l [+] ss.i [+] white [+] ses [+] male\nhed.l~pmt.l [+] pmt.i [+] ss.l [+] ss.i [+] white [+] ses [+] male",
                                     "hed.i~c(g1,g1)*pmt.l + pmt.i + ss.l + ss.i + white + ses\n
                                     hed.l~pmt.l + pmt.i + ss.l + ss.i + white + ses",m.sspmthed), data=ncanda_w1_w4, check.gradient=F,estimator="mlr",mimic="mplus",group="male")

sspmthed.sex.heds_pres<-growth(gsub("hed.i~pmt.l [+] pmt.i [+] ss.l [+] ss.i [+] white [+] ses [+] male\nhed.l~pmt.l [+] pmt.i [+] ss.l [+] ss.i [+] white [+] ses [+] male",
                                     "hed.i~pmt.l + pmt.i + ss.l + ss.i + white + ses\n
                                     hed.l~c(g1,g1)*pmt.l + pmt.i + ss.l + ss.i + white + ses",m.sspmthed), data=ncanda_w1_w4, check.gradient=F,estimator="mlr",mimic="mplus",group="male")

fitmeasures(sspmthed.sex)["chisq"]-fitmeasures(sspmthed.sex.hedi_ssi)["chisq"]

sspmthed.sex@test[[2]]$test
sspmthed.sex.hedi_ssi@test[[2]]$test

unlist(sspmthed.sex)
lavTestLRT(sspmthed.sex,sspmthed.sex.hedi_ssi)
lavTestLRT(sspmthed.sex,sspmthed.sex.heds_ssi,method="satorra.bentler.2001")
lavTestLRT(sspmthed.sex,sspmthed.sex.hedi_sss,method="satorra.bentler.2001")
lavTestLRT(sspmthed.sex,sspmthed.sex.heds_sss,method="satorra.bentler.2001")

lavTestLRT(sspmthed.sex,sspmthed.sex.hedi_prei,method="satorra.bentler.2001")
lavTestLRT(sspmthed.sex,sspmthed.sex.heds_prei,method="satorra.bentler.2001")
lavTestLRT(sspmthed.sex,sspmthed.sex.hedi_pres,method="satorra.bentler.2001")
lavTestLRT(sspmthed.sex,sspmthed.sex.heds_pres,method="satorra.bentler.2001")

######
#Latent Difference Score Growth Model
######

m.lds.full<-"
ss17=~ 1*upps23_ss_17r + upps31_ss_17r + upps36_ss_17r + upps46_ss_17r
ss18=~ 1*upps23_ss_18r + upps31_ss_18r + upps36_ss_18r + upps46_ss_18r
ss19=~ 1*upps23_ss_19r + upps31_ss_19r + upps36_ss_19r + upps46_ss_19r
ss20=~ 1*upps23_ss_20r + upps31_ss_20r + upps36_ss_20r + upps46_ss_20r

pre17=~ 1*upps16_pmt_17r + upps5_pmt_17r + upps28_pmt_17r + upps48_pmt_17r
pre18=~ 1*upps16_pmt_18r + upps5_pmt_18r + upps28_pmt_18r + upps48_pmt_18r
pre19=~ 1*upps16_pmt_19r + upps5_pmt_19r + upps28_pmt_19r + upps48_pmt_19r
pre20=~ 1*upps16_pmt_20r + upps5_pmt_20r + upps28_pmt_20r + upps48_pmt_20r

# upps23_ss_17r	~~	upps23_ss_17r
# upps31_ss_17r	~~	upps31_ss_17r
# upps36_ss_17r	~~	upps36_ss_17r
# upps46_ss_17r	~~	upps46_ss_17r
# upps23_ss_18r	~~	upps23_ss_18r
# upps31_ss_18r	~~	upps31_ss_18r
# upps36_ss_18r	~~  upps36_ss_18r
# upps46_ss_18r	~~	upps46_ss_18r
# upps23_ss_19r	~~	upps23_ss_19r
# upps31_ss_19r	~~	upps31_ss_19r
# upps36_ss_19r	~~	upps36_ss_19r
# upps46_ss_19r	~~	upps46_ss_19r
# upps23_ss_20r	~~	upps23_ss_20r
# upps31_ss_20r	~~	upps31_ss_20r
# upps36_ss_20r	~~	upps36_ss_20r
# upps46_ss_20r	~~	upps46_ss_20r
# upps23_ss_21	~~	upps23_ss_21
# upps31_ss_21	~~	upps31_ss_21
# upps36_ss_21	~~	upps36_ss_21
# upps46_ss_21	~~	upps46_ss_21
# upps16_pmt_17r	~~	upps16_pmt_17r
# upps5_pmt_17r	~~	upps5_pmt_17r
# upps28_pmt_17r	~~	upps28_pmt_17r
# upps48_pmt_17r	~~	upps48_pmt_17r
# 
# upps16_pmt_18r	~~	upps16_pmt_18r
# upps5_pmt_18r	~~	upps5_pmt_18r
# upps28_pmt_18r	~~	upps28_pmt_18r
# upps48_pmt_18r	~~	upps48_pmt_18r
# upps16_pmt_19r	~~	upps16_pmt_19r
# upps5_pmt_19r	~~	upps5_pmt_19r
# upps28_pmt_19r	~~	upps28_pmt_19r
# upps48_pmt_19r	~~	upps48_pmt_19r
# upps16_pmt_20r	~~	upps16_pmt_20r
# upps5_pmt_20r	~~	upps5_pmt_20r
# upps28_pmt_20r	~~	upps28_pmt_20r
# upps48_pmt_20r	~~	upps48_pmt_20r

#Item intercepts fixed to equivalence across ages

#Sensation Seeking
upps23_ss_17r	~	ssm*1
upps31_ss_17r	~	ssm*1
upps36_ss_17r	~	ssm*1
upps46_ss_17r	~	ssm*1
upps23_ss_18r	~	ssm*1
upps31_ss_18r	~	ssm*1
upps36_ss_18r	~  ssm*1
upps46_ss_18r	~	ssm*1
upps23_ss_19r	~	ssm*1
upps31_ss_19r	~	ssm*1
upps36_ss_19r	~	ssm*1
upps46_ss_19r	~	ssm*1
upps23_ss_20r	~	ssm*1
upps31_ss_20r	~	ssm*1
upps36_ss_20r	~	ssm*1
upps46_ss_20r	~	ssm*1

#Premeditation
upps16_pmt_17r	~	prem*1
upps5_pmt_17r	~	prem*1
upps28_pmt_17r	~	prem*1
upps48_pmt_17r	~	prem*1
upps16_pmt_18r	~	prem*1
upps5_pmt_18r	~	prem*1
upps28_pmt_18r	~	prem*1
upps48_pmt_18r	~	prem*1
upps16_pmt_19r	~	prem*1
upps5_pmt_19r	~	prem*1
upps28_pmt_19r	~	prem*1
upps48_pmt_19r	~	prem*1
upps16_pmt_20r	~	prem*1
upps5_pmt_20r	~	prem*1
upps28_pmt_20r	~	prem*1
upps48_pmt_20r	~	prem*1

#Item variances fixed to equivalence across ages

#Sensation Seeking
upps23_ss_17r	~~	i23*upps23_ss_17r
upps31_ss_17r	~~	i31*upps31_ss_17r
upps36_ss_17r	~~	i36*upps36_ss_17r
upps46_ss_17r	~~	i46*upps46_ss_17r
upps23_ss_18r	~~	i23*upps23_ss_18r
upps31_ss_18r	~~	i31*upps31_ss_18r
upps36_ss_18r	~~  i36*upps36_ss_18r
upps46_ss_18r	~~	i46*upps46_ss_18r
upps23_ss_19r	~~	i23*upps23_ss_19r
upps31_ss_19r	~~	i31*upps31_ss_19r
upps36_ss_19r	~~	i36*upps36_ss_19r
upps46_ss_19r	~~	i46*upps46_ss_19r
upps23_ss_20r	~~	i23*upps23_ss_20r
upps31_ss_20r	~~	i31*upps31_ss_20r
upps36_ss_20r	~~	i36*upps36_ss_20r
upps46_ss_20r	~~	i46*upps46_ss_20r

#Premeditation
upps16_pmt_17r	~~	i16*upps16_pmt_17r
upps5_pmt_17r	~~	i5*upps5_pmt_17r
upps28_pmt_17r	~~	i28*upps28_pmt_17r
upps48_pmt_17r	~~	i48*upps48_pmt_17r
upps16_pmt_18r	~~	i16*upps16_pmt_18r
upps5_pmt_18r	~~	i5*upps5_pmt_18r
upps28_pmt_18r	~~	i28*upps28_pmt_18r
upps48_pmt_18r	~~	i48*upps48_pmt_18r
upps16_pmt_19r	~~	i16*upps16_pmt_19r
upps5_pmt_19r	~~	i5*upps5_pmt_19r
upps28_pmt_19r	~~	i28*upps28_pmt_19r
upps48_pmt_19r	~~	i48*upps48_pmt_19r
upps16_pmt_20r	~~	i16*upps16_pmt_20r
upps5_pmt_20r	~~	i5*upps5_pmt_20r
upps28_pmt_20r	~~	i28*upps28_pmt_20r
upps48_pmt_20r	~~	i48*upps48_pmt_20r

ss17	~~	0*pre17
ss17	~~	0*ss18
ss17	~~	0*pre18
ss17	~~	0*ss19
ss17	~~	0*pre19
ss17	~~	0*ss20
ss17	~~	0*pre20
ss18	~~	0*pre17
pre17	~~	0*pre18
ss19	~~	0*pre17
pre17	~~	0*pre19
ss20	~~	0*pre17
pre17	~~	0*pre20
ss18	~~	0*pre18
ss18	~~	0*ss19
ss18	~~	0*pre19
ss18	~~	0*ss20
ss18	~~	0*pre20
ss19	~~	0*pre18
pre18	~~	0*pre19
ss20	~~	0*pre18
pre18	~~	0*pre20

ss19	~~	0*pre19
ss19	~~	0*ss20
ss19	~~	0*pre20

ss20	~~	0*pre19
pre19	~~	0*pre20

ss20	~~	0*pre20


####LDS MODEL; modelled off of Kievit et al (2018)


ss17 ~ 1*pre17 #regresses SS perfectly on PRE
da17 =~ 1*ss17 #Defines latent difference score factored as measured perfectly by SS
da17 ~ 1*0 #Estimates conditinoal mean of DA. Fix these and others to zero
pre17 ~ 1 #Estimates the mean of PRE
da17~~da17 #Variance of developmental imbalance
pre17~~pre17 #Variance of PRE
da17~~pre17 #self-feedback parameter

ss18 ~ 1*pre18
da18 =~ 1*ss18
da18 ~ 1*0
pre18 ~ 1
da18~~da18
pre18~~pre18
da18~~pre18

ss19 ~ 1*pre19
da19 =~ 1*ss19
da19 ~ 1*0
pre19 ~ 1
da19~~da19
pre19~~pre19
da19~~pre19

ss20 ~ 1*pre20
da20 =~ 1*ss20
da20 ~ 1*0
pre20 ~ 1
da20~~da20
pre20~~pre20
da20~~pre20


ss17 ~ 0*1 #Constrains intercept of SS to 0
ss17 ~~0*ss17 #fix variance of SS to 0
ss18 ~ 0*1
ss18 ~~0*ss18
ss19 ~ 0*1
ss19 ~~0*ss19
ss20 ~ 0*1
ss20 ~~0*ss20


da_i =~ 1*da17 + 1*da18 + 1*da19 + 1*da20
da_l =~ 0*da17 + 1*da18 + 2*da19 + 3*da20

da_i~~da_i
da_l~~da_l
da_i~~da_l

hed_i =~1*alc_bingeyr_ord_17 + 1*alc_bingeyr_ord_18 + 1*alc_bingeyr_ord_19 + 1*alc_bingeyr_ord_20
hed_s =~ -4*alc_bingeyr_ord_17 + -3*alc_bingeyr_ord_18 + -2*alc_bingeyr_ord_19 + -1*alc_bingeyr_ord_20

hed_i~1
hed_s~1
hed_i~~hed_i
hed_s~~hed_s
hed_i~~hed_s


alc_bingeyr_ord_17~0*1
alc_bingeyr_ord_18~0*1
alc_bingeyr_ord_19~0*1
alc_bingeyr_ord_20~0*1


alc_bingeyr_ord_17~~alc_bingeyr_ord_17
alc_bingeyr_ord_18~~alc_bingeyr_ord_18
alc_bingeyr_ord_19~~alc_bingeyr_ord_19
alc_bingeyr_ord_20~~alc_bingeyr_ord_20

hed_i~ da_i + da_l + white + ses + male
hed_s~ da_i + da_l + white + ses + male
da_i~1
da_l~1

"

lds<-lavaan(m.lds.full, data=ncanda_w1_w4, check.gradient=F,estimator="mlr", mimic="MPlus")

fitmeasures(lds)[c("df","chisq","rmsea","cfi","tli","srmr")]

#Stratified by sex
lds.sex<-lavaan(gsub("hed_i~ da_i [+] da_l [+] white [+] ses [+] male\nhed_s~ da_i [+] da_l [+] white [+] ses [+] male",
                     "hed_i~ da_i + da_l + white + ses\nhed_s~ da_i + da_l + white + ses",m.lds.full), data=ncanda_w1_w4, check.gradient=F,estimator="mlr",mimic="mplus",group="male")

fitmeasures(lds.sex)[c("df","chisq","rmsea","cfi","tli","srmr")]

#Sex moderation models
lds.sex.fix.hedi_dai<-lavaan(gsub("hed_i~ da_i [+] da_l [+] white [+] ses [+] male\nhed_s~ da_i [+] da_l [+] white [+] ses [+] male",
                                  "hed_i~ c(g1,g1)*da_i + da_l + white + ses\nhed_s~ da_i + da_l + white + ses",m.lds.full), 
                             data=ncanda_w1_w4, check.gradient=F,estimator="mlr",mimic="mplus",group="male")

lds.sex.fix.heds_dai<-lavaan(gsub("hed_i~ da_i [+] da_l [+] white [+] ses [+] male\nhed_s~ da_i [+] da_l [+] white [+] ses [+] male",
                                  "hed_i~ da_i + da_l + white + ses\nhed_s~ c(g1,g1)*da_i + da_l + white + ses",m.lds.full), 
                             data=ncanda_w1_w4, check.gradient=F,estimator="mlr",mimic="mplus",group="male")
lds.sex.fix.hedi_das<-lavaan(gsub("hed_i~ da_i [+] da_l [+] white [+] ses [+] male\nhed_s~ da_i [+] da_l [+] white [+] ses [+] male",
                                  "hed_i~ da_i + c(g1,g1)*da_l + white + ses\nhed_s~ da_i + da_l + white + ses",m.lds.full), 
                             data=ncanda_w1_w4, check.gradient=F,estimator="mlr",mimic="mplus",group="male")
lds.sex.fix.heds_das<-lavaan(gsub("hed_i~ da_i [+] da_l [+] white [+] ses [+] male\nhed_s~ da_i [+] da_l [+] white [+] ses [+] male",
                                  "hed_i~ da_i + da_l + white + ses\nhed_s~ da_i + c(g1,g1)*da_l + white + ses",m.lds.full), 
                             data=ncanda_w1_w4, check.gradient=F,estimator="mlr",mimic="mplus",group="male")

#Chi-square difference tests
lavTestLRT(lds.sex,lds.sex.fix.hedi_dai,method = "satorra.bentler.2001")
lavTestLRT(lds.sex,lds.sex.fix.heds_dai,method = "satorra.bentler.2001")
lavTestLRT(lds.sex,lds.sex.fix.hedi_das,method = "satorra.bentler.2001")
lavTestLRT(lds.sex,lds.sex.fix.heds_das,method = "satorra.bentler.2001")