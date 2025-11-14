rm(list=ls())

vector<-c("J01CE","J01CF","J01CA","J01A","J01MA","J01FF")
vector2<-c("scapis","simpler","mos")
#vector2<-c("simpler")
vector3<-c("richness","shannon","invsimpson")
#vector3<-c("invsimpson")

for (type in vector3){
  for (cohort in vector2){
    for (j in 1:250){
    #type="richness"
    #cohort="scapis"
    if (type=="richness"){
      low_lim=-50
      upp_lim=50
    }
    
    if (type=="shannon"){
      low_lim=-0.2
      upp_lim=0.2
    }
    
    setwd("/home/ulfha881/Desktop/Ulf/Antibiotics/Results")
    #file=list.files()
    #i=paste0("FDA_",cohort,"_",anti_type,".dta")
    #i="FDA_J01CE.dta"
    
    #for (i in file){
    library(data.table)
    library(fda)
    library(fda.usc)
    library(haven)
    library(fastDummies)
    library(car)
    data=read_dta(paste0("/proj/sens2019512/boot/FDA_",cohort,j,".dta"))
    #data$sex=factor(data$sex)
    data$education=factor(data$education)
    data$placebirth=factor(data$placebirth)
    data$site_plate=factor(data$site_plate)
    data$smokestatus=factor(data$smokestatus)
    data$polypharmacy12m_cat=factor(data$polypharmacy12m_cat)
    #data$polypharmacy12m_cat=factor(data$polypharmacy12m_cat)
    n=colnames(data)
    cds=c("J01CE","J01CF","J01CA","J01A","J01FF","J01MA","J01ZZ")
    for (vtype in cds){
      nm=grep(vtype,n,value=T)
      x=as.matrix(data[,nm])
      assign(paste0("x",vtype),fdata(x))
    }
    y=data.frame(data[[type]])
    
    data <- as.data.frame(data)
    if (cohort=="scapis") data=dummy_cols(data,select_columns=c("education","placebirth","site_plate","smokestatus","polypharmacy12m_cat"))
    if (cohort=="mos" | cohort=="simpler") data=dummy_cols(data,select_columns=c("education","site_plate","smokestatus","polypharmacy12m_cat"))
    
    #also "statins", "polypharmacy_12m","sex"
    #colnames(y)="richness"
    tt <- xJ01CA[["argvals"]]
    nbasis.x <- 15
    nbasis.b <- 5
    #,breaks=c(2,10,80)
    #,breaks=c(0,40,99)
    basis1 <- create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x,norder=4)
    basis2 <- create.bspline.basis(rangeval=range(tt),nbasis=nbasis.b,norder=4)
    basis.x <- list("xJ01CE"=basis1,"xJ01CF"=basis1,"xJ01CA"=basis1,"xJ01A"=basis1,
                    "xJ01FF"=basis1,"xJ01MA"=basis1,"xJ01ZZ"=basis1)
    basis.b <- list("xJ01CE"=basis2,"xJ01CF"=basis2,"xJ01CA"=basis2,"xJ01A"=basis2,
                    "xJ01FF"=basis2,"xJ01MA"=basis2,"xJ01ZZ"=basis2)
    #    cds=c("J01CE","J01CF","J01CA","J01A","J01FF","J01ZZ")
    #+ smokestatus + sex + education + placebirth + site_plate +diabd+bmi + 
    q=colnames(data)
    nm2=grep("^site_plate",q,value=T)
    nm2=nm2[-(1:2)]
    if (cohort=="scapis"){ 
      other_cov=c("age","bmi","sex","education_university","education_upper_secondary",
                  "smokestatus_former","smokestatus_current",
                  "placebirth_Other","placebirth_Asia","placebirth_Europe",
                  "polypharmacy12m_cat_1_4","polypharmacy12m_cat_5","cciw","ppi","statins","metformin",
                  "betablock","ssri","psych",
                  "xJ01CE","xJ01CF","xJ01CA","xJ01A","xJ01FF","xJ01MA","xJ01ZZ")
    }
    if (cohort=="mos" | cohort=="simpler"){
      other_cov=c("age","bmi","sex","education_university","education_upper_secondary",
                  "smokestatus_former","smokestatus_current",
                  "polypharmacy12m_cat_1_4","polypharmacy12m_cat_5","cciw","ppi","statins","metformin",
                  "betablock","ssri","psych",
                  "xJ01CE","xJ01CF","xJ01CA","xJ01A","xJ01FF","xJ01MA","xJ01ZZ")
    }
    
    nm2=c(nm2,other_cov)
    #nm2=tolower(nm2)
    #nm2=gsub(nm2,"ö","o")
    nm3=paste(nm2,collapse="+")
    #nm2=collapse(nm2,sep="+")
    #f <- richness ~ x + age + bmi + sex_male + education_University + education_Upper_secondary + nm2
    #f <- richness ~ age + bmi + sex_male + education_University + education_Upper_secondary +
    #  placebirth_Asia + placebirth_Europe + 
    #  placebirth_Other + x
    library(dplyr)
    data=rename(data,education_university=education_University)
    data=rename(data,education_upper_secondary=`education_Upper secondary`)
    
    #g <- richness ~ age + bmi + sex_male + education_university + education_upper_secondary
    f = as.formula(paste0(type,"~",nm3))
    
    #+ 
    #site_plate_malmoplate2_1+site_plate_malmoplate2_10+site_plate_malmoplate2_11+site_plate_malmoplate2_12+site_plate_malmoplate2_13+site_plate_malmoplate2_14+site_plate_malmoplate2_15+site_plate_malmoplate2_16+site_plate_malmoplate2_17+site_plate_malmoplate2_18+site_plate_malmoplate2_19+site_plate_malmoplate2_2+site_plate_malmoplate2_20+site_plate_malmoplate2_21+site_plate_malmoplate2_22+site_plate_malmoplate2_23+site_plate_malmoplate2_24+site_plate_malmoplate2_25+site_plate_malmoplate2_26+site_plate_malmoplate2_27+site_plate_malmoplate2_28+site_plate_malmoplate2_29+site_plate_malmoplate2_3+site_plate_malmoplate2_30+site_plate_malmoplate2_31+site_plate_malmoplate2_32+site_plate_malmoplate2_33+site_plate_malmoplate2_34+site_plate_malmoplate2_35+site_plate_malmoplate2_36+site_plate_malmoplate2_37+site_plate_malmoplate2_38+site_plate_malmoplate2_39+site_plate_malmoplate2_4+site_plate_malmoplate2_40+site_plate_malmoplate2_41+site_plate_malmoplate2_42+site_plate_malmoplate2_43+site_plate_malmoplate2_5+site_plate_malmoplate2_6+site_plate_malmoplate2_67+site_plate_malmoplate2_68+site_plate_malmoplate2_69+site_plate_malmoplate2_7+site_plate_malmoplate2_70+site_plate_malmoplate2_71+site_plate_malmoplate2_72+site_plate_malmoplate2_73+site_plate_malmoplate2_74+site_plate_malmoplate2_75+site_plate_malmoplate2_76+site_plate_malmoplate2_77+site_plate_malmoplate2_78+site_plate_malmoplate2_79+site_plate_malmoplate2_8+site_plate_malmoplate2_81+site_plate_malmoplate2_9+site_plate_uppsalaplate5_1+site_plate_uppsalaplate5_10+site_plate_uppsalaplate5_11+site_plate_uppsalaplate5_12+site_plate_uppsalaplate5_13+site_plate_uppsalaplate5_14+site_plate_uppsalaplate5_15+site_plate_uppsalaplate5_16+site_plate_uppsalaplate5_17+site_plate_uppsalaplate5_18+site_plate_uppsalaplate5_19+site_plate_uppsalaplate5_2+site_plate_uppsalaplate5_20+site_plate_uppsalaplate5_21+site_plate_uppsalaplate5_22+site_plate_uppsalaplate5_23+site_plate_uppsalaplate5_24+site_plate_uppsalaplate5_25+site_plate_uppsalaplate5_26+site_plate_uppsalaplate5_27+site_plate_uppsalaplate5_28+site_plate_uppsalaplate5_29+site_plate_uppsalaplate5_3+site_plate_uppsalaplate5_30+site_plate_uppsalaplate5_31+site_plate_uppsalaplate5_32+site_plate_uppsalaplate5_33+site_plate_uppsalaplate5_34+site_plate_uppsalaplate5_35+site_plate_uppsalaplate5_36+site_plate_uppsalaplate5_37+site_plate_uppsalaplate5_38+site_plate_uppsalaplate5_39+site_plate_uppsalaplate5_4+site_plate_uppsalaplate5_40+site_plate_uppsalaplate5_41+site_plate_uppsalaplate5_42+site_plate_uppsalaplate5_43+site_plate_uppsalaplate5_44+site_plate_uppsalaplate5_45+site_plate_uppsalaplate5_46+site_plate_uppsalaplate5_47+site_plate_uppsalaplate5_48+site_plate_uppsalaplate5_49+site_plate_uppsalaplate5_5+site_plate_uppsalaplate5_50+site_plate_uppsalaplate5_51+site_plate_uppsalaplate5_52+site_plate_uppsalaplate5_53+site_plate_uppsalaplate5_54+site_plate_uppsalaplate5_55+site_plate_uppsalaplate5_6+site_plate_uppsalaplate5_7+site_plate_uppsalaplate5_8+site_plate_uppsalaplate5_9+x
    ldat <- ldata("df"=data,"xJ01CE"=xJ01CE,"xJ01CF"=xJ01CF,"xJ01CA"=xJ01CA,
                  "xJ01A"=xJ01A,"xJ01FF"=xJ01FF,"xJ01MA"=xJ01MA,"xJ01ZZ"=xJ01ZZ)
    res <-  fregre.lm(f,data=ldat,basis.b=basis.b)
        nobs(res)

        summary(res)
        vector<-c("xJ01CE","xJ01CF","xJ01CA","xJ01A","xJ01FF","xJ01MA","xJ01ZZ")
        for (anti_type in vector){
          dt=numeric(0)
          summary(res)
          beta.l=res$beta.l[[anti_type]]
          str(beta.l)
          y_vals <- eval.fd(1:98, beta.l)
          d=cbind(y_vals,j)
          dt=rbind(d,dt)
          dt=data.frame(dt)
          write_dta(dt,paste0("/home/ulfha881/Desktop/Ulf/Antibiotics/Results/boot_data_",type,"_",cohort,"_",anti_type,j,".dta"))
      }
    }
  }
}
folder_path <- "/proj/sens2019512/boot/"
files <- list.files(folder_path, full.names = TRUE)
#file.remove(files)

#}