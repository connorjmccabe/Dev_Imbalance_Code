# Minor setup
packages.req <- c("semTools","dplyr")
new.packages <- packages.req[!(packages.req %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
require(dplyr)
require(semTools)
#Below read in the RDS file that contains all available fits from the lavaan object.
#These fits are then coverted to a list for easy reading.
fits.lds<-readRDS(url("https://github.com/connorjmccabe/Dev_Imbalance_Code/blob/master/DI_LDS_fits.rds?raw=true"))%>% as.list()


##the object fits.lds provides 72 fit measures that readers are welcomed/encouraged
# to examine. We report several in our manuscript, but describe below a consideration
#for interpreting relative fit measures (e.g., CFIs).
#Run the line of code below to view all available fit measures:

fits.lds

#Below, we compute scaled RMSEA value manually using the scaled variant given the use of
#MLR in analyses. We compare a manual computation to the estimated one as a "check"
#before computing this value for the baseline model.

sqrt(fits.lds$chisq.scaled-fits.lds$df.scaled)/sqrt(fits.lds$df.scaled*(fits.lds$ntotal-1))
fits.lds$rmsea.scaled#We can see these are essentially identical.

#Using the same formula, we can now compute this value for the baseline (null) model:
(fits.lds$rmsea.baseline.scaled<-sqrt(fits.lds$baseline.chisq.scaled-fits.lds$baseline.df.scaled)/sqrt(fits.lds$baseline.df.scaled*(fits.lds$ntotal-1)))

#The baseline RMSEA value was 0.1078, which is below the recommended RMSEA value
#of 0.15.

#As a more formal final check, we used the nullRMSEA function from semTools
#to compute this value. When supplied with our estimated model, results
#were nearly identical, as follows:

#> nullRMSEA(lds.final,scaled=T)
#The baseline model's RMSEA = 0.1112

#CFI, TLI, and other incremental fit indices may not be very informative because the baseline model's RMSEA < 0.158 (Kenny, Kaniskan, & McCoach, 2015).

