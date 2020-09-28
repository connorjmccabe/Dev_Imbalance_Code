#"Developmental imbalance between sensation seeking and premeditation in adolescence predicts greater heavy episodic drinking in young adulthood"

###Analyses Code####

#####NOTES

#This code provides comprehensive code for all SEM models described in the manuscript.
#These are re-constructed using the covariance matrix from the sample data given that we were unable to share the data directly.
#As such, fit measures and standard errors may not be directly comparable to the mmodel results estimated on the observed sample. Nonetheless, estimates provided by the code below will reflect the associations reported in the manuscript.


###Minor setup code

packages.req <- c("lavaan","RCurl")
new.packages <- packages.req[!(packages.req %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

require(lavaan)
require(RCurl)
########SEMS#####

#NOTE:

#upps_ss = UPPS Sensation Seeking
#upps_pmt = UPPS Premeditatino
#alc_bingeyr_ord = Past Year Heavy Episodic Drinking Quantity
########
#Unconditional Growth Model of Sensation Seeking and Premeditation
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

#Read in covariance matrices and measure means:
uncon.full<-readRDS(url("https://github.com/connorjmccabe/Dev_Imbalance_Code/blob/master/uncon_full.RDS?raw=true"))

#Fit using the full sample (i.e. not stratified by sex)
uncon.matinput<-growth(model=m.uncon,
                       sample.cov=uncon.full$cov,
                       sample.mean=uncon.full$mean,
                       sample.nobs=758,
                       mimic="Mplus")

#Fit Measures
#NOTE: these will be biased toward exceptionally-fitting models due to the use of covariances rather than observed data.

fitmeasures(uncon.matinput)[c("df","chisq","rmsea","cfi","tli","srmr")]

summary(uncon.matinput,ci=T,standardized=T)

#Model stratefied by sex
m.uncon.sex<-"

ss.l=~ 0*upps_ss_16r + 1*upps_ss_17r + 2*upps_ss_18r + 3*upps_ss_19r + 4*upps_ss_20r
#ss.q=~ 0*upps_ss_16 + 1*upps_ss_17 + 4*upps_ss_18 + 9*upps_ss_19 + 16*upps_ss.s_20 + 5*upps_ss_21
ss.i=~ 1*upps_ss_16r + 1*upps_ss_17r + 1*upps_ss_18r + 1*upps_ss_19r + 1*upps_ss_20r

pmt.l=~ 0*upps_pmt_16r + 1*upps_pmt_17r + 2*upps_pmt_18r + 3*upps_pmt_19r + 4*upps_pmt_20r
#pmt.q=~ 0*upps_pmt_16r + 1*upps_pmt_17r + 4*upps_pmt_18r + 9*upps_pmt_19r + 16*upps_pmt_20r
pmt.i=~ 1*upps_pmt_16r + 1*upps_pmt_17r + 1*upps_pmt_18r + 1*upps_pmt_19r + 1*upps_pmt_20r

pmt.i~~pmt.l
ss.i~~ss.l
ss.i~~pmt.i
pmt.l~~ss.i

pmt.i~c(pmtiF,pmtiM)*1
pmt.l~c(pmtlF,pmtlM)*1
ss.i~c(ssiF,ssiM)*1
ss.l~c(sslF,sslM)*1

#Sex moderation parameters

sexmod_pmti := pmtiM - pmtiF
sexmod_pmtl := pmtlM - pmtlF
sexmod_ssi := ssiM - ssiF
sexmod_ssl := sslM - sslF
"

uncon.sex<-readRDS(url("https://github.com/connorjmccabe/Dev_Imbalance_Code/blob/master/uncon_sex.RDS?raw=true"))
uncon.sex.matinput<-growth(model=m.uncon.sex,
                           sample.cov=list(uncon.sex$`0`$cov,uncon.sex$`1`$cov),
                           sample.mean=list(uncon.sex$`0`$mean,uncon.sex$`1`$mean),
                           sample.nobs=list(387, 371),
                           mimic="Mplus")

fitmeasures(uncon.sex.matinput)[c("df","chisq","rmsea","cfi","tli","srmr")]

summary(uncon.sex.matinput,ci=T,standardized=T)

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

sspmthed<-readRDS(url("https://github.com/connorjmccabe/Dev_Imbalance_Code/blob/master/sspmthed_full.RDS?raw=true"))

sspmthed.matinput<-growth(m.sspmthed,
                          sample.cov=sspmthed$cov,
                          sample.mean=sspmthed$mean,
                          sample.nobs=774,
                          mimic="MPlus")

fitmeasures(sspmthed.matinput)[c("chisq","rmsea","cfi","tli", "srmr")]
summary(sspmthed.matinput,standardized=T,ci=T)

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

#Item intercepts fixed to equivalence across ages

#Sensation Seeking
upps23_ss_17r	~	ssm*1
upps31_ss_17r	~	ssm*1
upps36_ss_17r	~	ssm*1
upps46_ss_17r	~	ssm*1
upps23_ss_18r	~	ssm*1
upps31_ss_18r	~	ssm*1
upps36_ss_18r	~ ssm*1
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
hed_s =~ -3*alc_bingeyr_ord_17 + -2*alc_bingeyr_ord_18 + -1*alc_bingeyr_ord_19 + 0*alc_bingeyr_ord_20

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

lds.full<-readRDS(url("https://github.com/connorjmccabe/Dev_Imbalance_Code/blob/master/lds_full.RDS?raw=true"))



lds.matinput<-lavaan(model=m.lds.full,
       sample.cov=lds.full$cov,
       sample.mean=lds.full$mean,
       sample.nobs=774)

fitmeasures(lds.matinput)[c("df","chisq","rmsea","cfi","tli","srmr")]

summary(lds.matinput,standardized=T,ci=T)
