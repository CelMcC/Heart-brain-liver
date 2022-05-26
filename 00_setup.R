

myPackages<- c("tidyverse",
               "stringr",
               "biglm",
               "data.table",
               "readxl",
               "writexl",
               "janitor",
               "skimr",
               "viridis",
               "car",
               "lubridate",
               "e1071",
               "broom",
               "ggh4x",
               "UKBTools")

pacman::p_load(myPackages, character.only = TRUE)


expectedRef<- data.frame(stringsAsFactors = FALSE,
                           ExpEF = c("healthy","healthy","adverse","adverse",
                                     "healthy","healthy","adverse","adverse"),
                           OutEF = c("healthy","adverse","healthy","adverse",
                                     "healthy","adverse","healthy","adverse"),
                       Direction = c("Pos","Pos","Pos","Pos","Neg","Neg",
                                     "Neg","Neg"),
                        Expected = c("Expected","Unexpected","Unexpected",
                                     "Expected","Unexpected","Expected","Expected",
                                     "Unexpected"))

get_checksum<- function(df){
  # Quick way to check if data.frame has returned the same data as last time
  df<-         data.frame(df)
  N<-          nrow(df)
  check_idx<-  seq(1, N, length.out= 100) %>% floor() %>% unique()
  is_numeric<- which(sapply(df, class) %like% "integer|numeric")
  res<-        sum(df[check_idx, is_numeric], na.rm=TRUE)
  return(res)
}



mod_format<-         "%2.3f"
p_format<-           "%2.4f"
exp_format<-         "%1.2E"
skewness_threshhold<- 1

## Some functions

compile_formulas<- function(t0_start, imaging_confounders= FALSE) {
  
  if (imaging_confounders== TRUE) {
    cfs<- paste0(" + ", confound_string)
  } else {
    cfs<- ""
  }
  t0<- t0_start %>%
    mutate(covString = 
             case_when(Type %like% "1" ~  paste0(str_c(c("age_2", "age_sq", "sex"), collapse = " + "), cfs), 
                       Type %like% "2" ~  paste0(str_c(covariates1, collapse= " + "), cfs),
                       Type %like% "3" ~  paste0(str_c(covariates2, collapse= " + "), cfs)),
           form= paste(Outcomes, "~", Exposures, "+", covString),
           ID= as.character(1:nrow(t0_start)),
           formulas= map(form, as.formula))
  print(t0$formulas[[nrow(t0)]])
  return(t0)
}


apply_fdr_to_set<- function(tidy_object){
    x<- tidy_object %>% mutate(fdr= p.adjust(p.value, method= "BH"),
                               sig_by_fdr= ifelse(fdr < 0.05, "*", ""))
    return(x)
}

formatP<- . %>% mutate(PValue = ifelse(p.value < 0.001, 
                                       sprintf(exp_format, p.value), 
                                       sprintf(p_format, p.value)), 
                       PValue = str_replace(PValue, "E-0", "E-"), 
                       PValue = str_replace(PValue, "E-", "x10-"))


expand_lin<- . %>% mutate(Sig= ifelse(sig_by_fdr=="*", 1, 0.25),
                          Direction= ifelse(estimate > 0, "Pos", "Neg"),
                          Comment1= paste0(sprintf(mod_format, estimate), sig_by_fdr),
                          Comment2= paste0("[", sprintf(mod_format, conf.low), ", ",
                                           sprintf(mod_format, conf.high), "]"))

report_cutoff<- function(res) {
  themax<- max(res$p.value[res$sig_by_fdr=="*"], na.rm=TRUE) %>% round(3)
  cat("Sig threshhold sits at",themax)
  return(res)
}

standard_linear<- function(tidy_object){
  res<- tidy_object %>% left_join(t0, by= "ID") %>%
    apply_fdr_to_set() %>%
    expand_lin() %>% formatP() %>%
    mutate(Analysis= titleText,
           is_exposure= (term== Exposures),
           Comment3= PValue, Comment4= N) %>%
    left_join(proc %>% select(term= Var, ExpName= CovName, ExpEF= Effect), by= "term") %>%
    left_join(proc %>% select(Outcomes= Var, OutName= CovName, OutEF= Effect), by= "Outcomes") %>%
    left_join(expectedRef, by= c("Direction", "ExpEF", "OutEF")) %>% report_cutoff()
  return(res)
}


safelm<- function(i, dataset){
  targetvar<-      t0$Exposures[i]
  dataset$target<- dataset[, targetvar]
  tempset<-        dataset %>% filter(!is.na(target))
  fmla<-           t0$formulas[[i]]
  mdl<-            biglm(fmla, tempset) 
  N<-              mdl$n
  tid<-            mdl %>% tidy(conf.int= TRUE) %>% mutate(N= N)
  return(tid)
  
}


selective_scale<- function(sett){
  ## Scale only the numeric fields (leave out f.eid, survival and binary fields)
  
  ## Work out which fields to scale
  feid_idx<-       which(names(sett)=="f.eid")
  binaries<-       which(sapply(sett, function(x) howmany(na.omit(x)))==2)
  non_numeric<-    which(sapply(sett, function(x) !is.numeric(x)))
  surv_cols<-      which(names(sett) %like% "_surv")
  non_scale_cols<- unique(c(feid_idx, binaries, non_numeric, surv_cols)) %>% sort()
  scale_cols<-     setdiff(1:ncol(sett), non_scale_cols)
  
  # Create scaled sett
  scaled_sett<-     sett
  scaled_sett[, scale_cols]<- scale(scaled_sett[, scale_cols]) %>% data.frame()
  return(scaled_sett)
}

putHere<- function(the_object){
  nm <-deparse(substitute(the_object))
  saveRDS(the_object, paste0(thisFileFolder, nm, ".rds"))
}

getHere<- function(object_name){
  the_object<- readRDS(paste0(thisFileFolder, object_name, ".rds"))
  return(the_object)
}

mystrips1 <- strip_nested(
  # Horizontal strips
  background_x = elem_list_rect(fill = c("white", "white", "grey90", "grey90", "grey90", "grey90", "grey90", "grey90"),
                                color= rep("grey10", 8),
                                size= rep(0.1, 8)),
  text_x = elem_list_text(colour = "black",
                          family= "Avenir Heavy"),
  by_layer_x = FALSE,
  # Vertical strips
  background_y = element_blank(),
  text_y = elem_list_text(angle = c(0, 90)),
  by_layer_y = FALSE
)

mystrips2 <- strip_nested(
  # Horizontal strips
  background_x = elem_list_rect(fill = c("white",  "white", "grey90", "grey90", "grey90", "grey90"),
                                color= c("grey10", rep("grey10", 6)),
                                size= c(0.1, rep(0.1, 5))),
  text_x = elem_list_text(colour = c("black", rep("black", 6)),
                          family= rep("Avenir Heavy", 6)),
  by_layer_x = FALSE,
  # Vertical strips
  background_y = element_blank(),
  text_y = elem_list_text(angle = c(0, 90)),
  by_layer_y = FALSE
)

fmtc<- function(x) {
  format(x, big.mark= ",", digits= 0, format= "f")
}

table1_fun<- function(tt){
  
  # Extract N's by exposure and format as a range
  n0<- tt %>% group_by(ExpName) %>% summarise(minN= min(N), maxN= max(N)) %>%
    mutate(Nbracket= paste0("N = ", fmtc(minN), " - ", fmtc(maxN))) %>%
    select(ExpName, Nbracket)
  
  tx<- tt %>%
    select(ExpName, Outcomes, Comment1, Comment2, PValue) %>%
    left_join(proc %>% select(Outcomes= Var, OutName= CovName)) %>%
    select(-Outcomes) %>%
    gather(key, value, -c(ExpName, OutName)) %>%
    spread(OutName, value) %>% left_join(n0)
  
  return(tx)
}

myTheme<- theme_light() + theme(
  panel.grid = element_blank(),
  axis.ticks.y = element_blank(),
  text = element_text(family= "serif"),
  strip.background = element_rect(fill= "grey20"),
  strip.text= element_text(face= "bold", color= "white"))
theme_set(myTheme)

asN<- function(number){
  fn<- format(number, big.mark = ",") %>% trimws()
  out<- paste0("\n(n = ",fn,")")
  return(out)
}

numeric_summary<- function(df, std_format, med_format) {
  res<- df %>% 
    summarise(Mean = mean(number, na.rm = TRUE),
              SD = sd(number, na.rm = TRUE),
              Median = median(number, na.rm =TRUE),
              skew = abs(e1071::skewness(number, na.rm= TRUE)),
              Q25 = quantile(number, p = 0.25, na.rm = TRUE),
              Q75 = quantile(number, p = 0.75, na.rm = TRUE),
              Missing = sum(is.na(number)),
              N = n()) %>%
    mutate(Comment1 = ifelse(
      summary_use== "median",
      paste0( sprintf(med_format, Median),
              " [", sprintf(med_format, Q25),
              ", ",  sprintf(med_format, Q75),  "]" ),
      paste0(sprintf(std_format, Mean), " (±", sprintf(std_format, SD) , ")")))
  return(res)
}


categorical_summary<- function(df, std_format){
res<- df %>% summarise(Count= n()) %>% 
  add_count(Set, CovName, descTable, wt = Count, name= "Varcount") %>%
  mutate(Percent= Count/Varcount,
         value= ifelse(is.na(value), "(Missing)", value)) %>% 
  mutate(Comment1= paste0(format(Count, big.mark = ","), " (",sprintf(std_format, Percent*100),"%)"),
         CovName_temp= paste(CovName, value, sep= ": ")) %>%
  filter(is.na(Collapse) | (Collapse=="yes" & value== 1))
  return(res)
}


flat_frame<- function(named_vec){
  out<- data.frame(Name= names(named_vec), Value= named_vec)
  return(out)
}

