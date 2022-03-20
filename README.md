1. Overview
The code in this replication packages constructs our analysis using R Studio. There are eight figures and two tables in this project. The following packages is needed before running the code.
library(readr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(kableExtra)
library(car)
library(patchwork)
library(gcookbook)  

2. Data Availability and Provenance Statements
2.1 Statement about Rights We certify that the authors of the manuscript have legitimate access to and permission to use the data used in this manuscript.

2.2 Summary of Availability All data in the replication package is available here, though some births data relies on aggregations from restricted-use microdata

datamain <- read_csv("AA6vcBgY.csv")
datatit<- read_file("AAOnAZHJ.txt")
dict <- read_lines("gss_dict.txt", skip = 18) 

3.Here are the instructions and codes to replicate the graph. 
3.1 By the following data, mutate a more cleaning data called gss.
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(readr)
library(tidyverse)
library(janitor)
```
```{r}
datamain <- read_csv("AA6vcBgY.csv")
datatit<- read_file("AAOnAZHJ.txt")
dict <- read_lines("gss_dict.txt", skip = 18) 
```
```{r}
variable_descriptions <- as_tibble(dict) %>% 
  filter(value!="}") %>% 
  mutate(value = str_replace(value, ".+%[0-9].*f[ ]{2,}", "")) %>% 
  mutate(value = str_remove_all(value, "\"")) %>% 
  rename(variable_description = value) %>% 
  bind_cols(tibble(variable_name = colnames(datamain)[-1]))
variable_descriptions
```

```{r,warning=FALSE}
# Now we want a variable name and the possible values
labels_raw_tibble <- as_tibble(str_split(datatit, ";")[[1]]) %>% 
  filter(row_number()!=1) %>% 
  mutate(value = str_remove(value, "\nlabel define ")) %>% 
  mutate(value = str_replace(value, "[ ]{2,}", "XXX")) %>% 
  mutate(splits = str_split(value, "XXX")) %>% 
  rowwise() %>% 
  mutate(variable_name = splits[1], cases = splits[2]) %>% 
  mutate(cases = str_replace_all(cases, "\n [ ]{2,}", "")) %>%
  select(variable_name, cases) %>% 
  drop_na()
labels_raw_tibble <- labels_raw_tibble %>% 
  mutate(splits = str_split(cases, "[ ]{0,}\"[ ]{0,}"))
add_cw_text <- function(x, y){
  if(!is.na(as.numeric(x))){
    x_new <- paste0(y, "==", x,"~")
  }
  else{
    x_new <- paste0("\"",x,"\",")
  }
  return(x_new)
}
cw_statements <- labels_raw_tibble %>% 
  rowwise() %>% 
  mutate(splits_with_cw_text = list(modify(splits, add_cw_text, y = variable_name))) %>% 
  mutate(cw_statement = paste(splits_with_cw_text, collapse = "")) %>% 
  mutate(cw_statement = paste0("case_when(", cw_statement,"TRUE~\"NA\")")) %>% 
  mutate(cw_statement = str_replace(cw_statement, ",\"\",",",")) %>% 
  select(variable_name, cw_statement)
```
```{r}
cw_statements <- 
  cw_statements %>% 
  mutate(variable_name = str_remove_all(variable_name, "\\r")) %>% 
  mutate(cw_statement = str_remove_all(cw_statement, "\\r"))
```
```{r}
gss <- datamain %>% 
  select(CASEID, 
         agedc, 
         achd_1c, 
         achdmpl, 
         totchdc, 
         acu0c,
         agema1c,
         achb1c,
         rsh_131a,
         arretwk,
         slm_01, 
         sex, 
         brthcan, 
         brthfcan,
         brthmcan,
         brthmacr,
         brthprvc,
         yrarri,
         prv, 
         region, 
         luc_rst, 
         marstat, 
         amb_01, 
         vismin, 
         alndimmg,
         bpr_16, 
         bpr_19,
         ehg3_01b, 
         odr_10, 
         livarr12, 
         dwelc, 
         hsdsizec,
         brthpcan,
         brtpprvc, 
         visminpr,
         rsh_125a, 
         eop_200,
         uhw_16gr,
         lmam_01, 
         acmpryr,
         srh_110,
         srh_115,
         religflg, 
         rlr_110,
         lanhome, 
         lan_01,
         famincg2, 
         ttlincg2, 
         noc1610, 
         cc_20_1,
         cc_30_1,
         ccmoc1c,
         cor_031,
         cor_041,
         cu0rnkc,
         pr_cl,
         chh0014c,
         nochricc,
         grndpa,
         gparliv,
         evermar,
         ma0_220,
         nmarevrc,
         ree_02,
         rsh_131b,
         rto_101,
         rto_110,
         rto_120,
         rtw_300,
         sts_410,
         csp_105,
         csp_110a,
         csp_110b,
         csp_110c,
         csp_110d,
         csp_160,
         fi_110) %>% 
  mutate_at(vars(agedc:fi_110), .funs = funs(ifelse(.>=96, NA, .))) %>% 
  mutate_at(.vars = vars(sex:fi_110),
            .funs = funs(eval(parse(text = cw_statements %>%
                                      filter(variable_name==deparse(substitute(.))) %>%
                                      select(cw_statement) %>%
                                      pull()))))

# Fix the names
gss <- gss %>% 
  clean_names() %>% 
  rename(age = agedc,
         age_first_child = achd_1c,
         age_youngest_child_under_6 = achdmpl,
         total_children = totchdc,
         age_start_relationship = acu0c,
         age_at_first_marriage = agema1c,
         age_at_first_birth = achb1c,
         distance_between_houses = rsh_131a,
         age_youngest_child_returned_work = arretwk,
         feelings_life = slm_01,
         sex = sex,
         place_birth_canada = brthcan,
         place_birth_father = brthfcan,
         place_birth_mother = brthmcan,
         place_birth_macro_region = brthmacr,
         place_birth_province = brthprvc,
         year_arrived_canada = yrarri,
         province = prv,
         region = region,
         pop_center = luc_rst,
         marital_status = marstat,
         aboriginal = amb_01,
         vis_minority = vismin,
         age_immigration = alndimmg,
         landed_immigrant = bpr_16,
         citizenship_status = bpr_19,
         education = ehg3_01b,
         own_rent = odr_10,
         living_arrangement = livarr12,
         hh_type = dwelc,
         hh_size = hsdsizec,
         partner_birth_country = brthpcan,
         partner_birth_province = brtpprvc,
         partner_vis_minority = visminpr,
         partner_sex = rsh_125a,
         partner_education = eop_200,
         average_hours_worked = uhw_16gr,
         worked_last_week = lmam_01,
         partner_main_activity = acmpryr,
         self_rated_health = srh_110,
         self_rated_mental_health = srh_115,
         religion_has_affiliation = religflg,
         regilion_importance = rlr_110,
         language_home = lanhome,
         language_knowledge = lan_01,
         income_family = famincg2,
         income_respondent = ttlincg2,
         occupation = noc1610,
         childcare_regular = cc_20_1,
         childcare_type = cc_30_1,
         childcare_monthly_cost = ccmoc1c,
         ever_fathered_child = cor_031,
         ever_given_birth = cor_041,
         number_of_current_union = cu0rnkc,
         lives_with_partner = pr_cl,
         children_in_household = chh0014c,
         number_total_children_intention = nochricc,
         has_grandchildren = grndpa,
         grandparents_still_living = gparliv,
         ever_married = evermar,
         current_marriage_is_first = ma0_220,
         number_marriages = nmarevrc,
         religion_participation = ree_02,
         partner_location_residence = rsh_131b,
         full_part_time_work = rto_101,
         time_off_work_birth = rto_110,
         reason_no_time_off_birth = rto_120,
         returned_same_job = rtw_300,
         satisfied_time_children = sts_410,
         provide_or_receive_fin_supp = csp_105,
         fin_supp_child_supp = csp_110a,
         fin_supp_child_exp = csp_110b,
         fin_supp_lump = csp_110c,
         fin_supp_other = csp_110d,
         fin_supp_agreement = csp_160,
         future_children_intention = fi_110) 

#### Clean up ####
gss <- gss %>% 
  mutate_at(vars(age:future_children_intention), 
            .funs = funs(ifelse(.=="Valid skip"|.=="Refusal"|.=="Not stated", "NA", .))) 

gss <- gss %>% 
  mutate(is_male = ifelse(sex=="Male", 1, 0)) 

gss <- gss %>% 
  mutate_at(vars(fin_supp_child_supp:fin_supp_other), .funs = funs(case_when(
    .=="Yes"~1,
    .=="No"~0,
    .=="NA"~as.numeric(NA)
  )))

main_act <- datamain %>% 
  mutate(main_activity = case_when(
    mpl_105a=="Yes"~ "Working at a paid job/business",
    mpl_105b=="Yes" ~ "Looking for paid work",
    mpl_105c=="Yes" ~ "Going to school",
    mpl_105d=="Yes" ~ "Caring for children",
    mpl_105e=="Yes" ~ "Household work", 
    mpl_105i=="Yes" ~ "Other", 
    TRUE~ "NA")) %>% 
  select(main_activity) %>% 
  pull()

age_diff <- datamain %>% 
  select(marstat, aprcu0c, adfgrma0) %>% 
  mutate_at(.vars = vars(aprcu0c:adfgrma0),
            .funs = funs(eval(parse(text = cw_statements %>%
                                      filter(variable_name==deparse(substitute(.))) %>%
                                      select(cw_statement) %>%
                                      pull())))) %>% 
  mutate(age_diff = ifelse(marstat=="Living common-law", aprcu0c, adfgrma0)) %>% 
  mutate_at(vars(age_diff), .funs = funs(ifelse(.=="Valid skip"|.=="Refusal"|.=="Not stated", "NA", .))) %>% 
  select(age_diff) %>% 
  pull()

gss <- gss %>% mutate(main_activity = main_act, age_diff = age_diff)

# Change some from strings into numbers
gss <- gss %>% 
  rowwise() %>% 
  mutate(hh_size = str_remove(string = hh_size, pattern = "\\ .*")) %>% 
  mutate(hh_size = case_when(
    hh_size=="One" ~ 1,
    hh_size=="Two" ~ 2,
    hh_size=="Three" ~ 3,
    hh_size=="Four" ~ 4,
    hh_size=="Five" ~ 5,
    hh_size=="Six" ~ 6
  )) 



gss <- gss %>% 
  rowwise() %>% 
  mutate(number_marriages = str_remove(string = number_marriages, pattern = "\\ .*")) %>% 
  mutate(number_marriages = case_when(
    number_marriages=="No" ~ 0,
    number_marriages=="One" ~ 1,
    number_marriages=="Two" ~ 2,
    number_marriages=="Three" ~ 3,
    number_marriages=="Four" ~ 4
  )) 

gss <- gss %>% 
  rowwise() %>% 
  mutate(number_total_children_known = ifelse(number_total_children_intention=="Don't know"|number_total_children_intention=="NA", 0, 1)) %>% 
  mutate(number_total_children_intention = str_remove(string = number_total_children_intention, pattern = "\\ .*")) %>% 
  mutate(number_total_children_intention = case_when(
    number_total_children_intention=="None" ~ 0,
    number_total_children_intention=="One" ~ 1,
    number_total_children_intention=="Two" ~ 2,
    number_total_children_intention=="Three" ~ 3,
    number_total_children_intention=="Four" ~ 4,
    number_total_children_intention=="Don't" ~ as.numeric(NA)
  )) 

write_csv(gss, "gss.csv")
```
```{r}
dataplot<-age_mom
i=1
tmp <- matrix(0, nr=61, nc=5)
agerange <- range(dataplot$age_first_child)
for(age in agerange[1]:agerange[2] ){
  # tmp[i] <- c(age, )
  tmpdata<- dataplot$par_age_range[dataplot$age_first_child==age]
  tmp[i,] <- c(age, 
              sum(tmpdata=="20"),
              sum(tmpdata=="20-30"), 
              sum(tmpdata=="30-40"), 
              sum(tmpdata=="40-50"))
  i<-i+1
}

colnames(tmp)<- c("age_first_child","<20","20-30","30-40", "40-50")
tmp <- as.data.frame(tmp)
write.csv(tmp, "tmp.csv")
class(tmp)
```
3.2
By using gss we can plot the graph and make a table 

 Figure 1:
 different feelings_life group
```{r,echo=FALSE,error=FALSE,warning=FALSE}
 ggplot(gss, aes(x= feelings_life)) +
  geom_bar(fill="#19b3a2") +
  theme_classic() +
  labs(y= "number of response", title="Fig.1 different feelings_life group")
feel <- gss %>% select(feelings_life)%>% filter(!is.na(feelings_life))
```
 Figure 2:
 six interested relevent groups
 ```{r,echo=FALSE, warning=FALSE, fig.height = 12, fig.width = 12}
gss_data <- gss %>% filter(!is.na(income_respondent)&!is.na(education)&!is.na(self_rated_health)&!is.na(self_rated_mental_health)&!is.na(average_hours_worked)&!is.na(average_hours_worked)&!is.na(marital_status)& self_rated_health != "Don't know"& self_rated_mental_health != "Don't know"&average_hours_worked != "Don't know"& marital_status != "Don't know")
g1<-ggplot(gss_data, aes(y= income_respondent)) +
  geom_bar(fill="#69b3a2") +
  theme_classic() +
  labs(x= "number of response", title="different income_respondent group")
g2<-ggplot(gss_data, aes(y= education)) +
  geom_bar(fill="#69b3a2") +
  theme_classic() +
  labs(x= "number of response", title="different education group")
g3 <-ggplot(gss_data, aes(y= self_rated_health)) +
  geom_bar(fill="#69b3a2") +
  theme_classic() +
  labs(x= "number of response", title="different self_rated_health group")
g4 <-ggplot(gss_data, aes(y= self_rated_mental_health)) +
  geom_bar(fill="#69b3a2") +
  theme_classic() +
  labs(x= "number of response", title="different self_rated_mental_health group")
g5 <-ggplot(gss_data, aes(y= average_hours_worked)) +
  geom_bar(fill="#69b3a2") +
  theme_classic() +
  labs(x= "number of response", title="different average_hours_worked group")
g6 <-ggplot(gss_data, aes(y= marital_status)) +
  geom_bar(fill="#69b3a2") +
  theme_classic() +
  labs(x= "number of response", title="different marital_status group")
g_all <- (g1 + g2) /
(g3 + g4) /
(g5 + g6)  +    # Create grid of plots with title
  plot_annotation(title = "fig.2 six interested relevent groups") & 
  theme(plot.title = element_text(hjust = 0.5))
g_all    
```

Figure.3:
relationship between age and feelings_life
```{r}
feeling_data <- gss %>% select(feelings_life, age) %>% filter(!is.na(feelings_life))
feeling_data$age <- round(feeling_data$age) 
feelingdata <- feeling_data %>% group_by(age) %>% summarise(mean_feelings_life = mean(feelings_life))
ggplot(feelingdata, aes(x=age, y = mean_feelings_life)) +
geom_line()+ 
  theme(panel.grid = element_line(color = "#8ccde3",
                                  size = 0.75,
                                  linetype = 2))+
  labs(x = "age", y = "average feelings of life from 0 to 10", title = "fig.3:relationship between age and feelings_life") 
```
Figure.4:
six interested relevent groups
```{r,echo=FALSE,error=FALSE,warning=FALSE}
a=c(1,2,3,4,5,6,7,8,9,10)
income_data <- gss %>% select(income_respondent, feelings_life) %>% filter(!is.na(feelings_life) & !is.na(income_respondent))
incomedata <- income_data %>%group_by(feelings_life) %>% summarise("<$25,000" = sum(income_respondent == "Less than $25,000"),"$25,000 to $49,999" = sum(income_respondent == "$25,000 to $49,999"), "$50,000 to $74,999" = sum(income_respondent == "$50,000 to $74,999"),"$75,000 to $99,999" = sum(income_respondent == "$75,000 to $99,999"), "$100,000 to $124,999" = sum(income_respondent == "$100,000 to $ 124,999"),"$125,000 and more" = sum(income_respondent == "$125,000 and more"))
incomedata$`<$25,000` = incomedata$`<$25,000`/sum(incomedata$`<$25,000`)
incomedata$`$25,000 to $49,999` = incomedata$`$25,000 to $49,999`/sum(incomedata$`$25,000 to $49,999`)
incomedata$`$50,000 to $74,999` = incomedata$`$50,000 to $74,999`/sum(incomedata$`$50,000 to $74,999`)
incomedata$`$75,000 to $99,999` = incomedata$`$75,000 to $99,999`/sum(incomedata$`$75,000 to $99,999`)
incomedata$`$100,000 to $124,999` = incomedata$`$100,000 to $124,999`/sum(incomedata$`$100,000 to $124,999`)
incomedata$`$125,000 and more` = incomedata$`$125,000 and more`/sum(incomedata$`$125,000 and more`)
dfb<-gather(incomedata,key=income_group,value=percent_of_Respondent,c("<$25,000","$25,000 to $49,999","$50,000 to $74,999","$75,000 to $99,999","$100,000 to $124,999","$125,000 and more"))
gfa<- ggplot(dfb, aes(x=feelings_life, y = percent_of_Respondent, group = income_group, colour = income_group)) +
  geom_line() +
  theme_classic() + scale_x_continuous(breaks = a)+
  labs( title="different income group relevant to feelings_life in percent")
```

```{r,echo=FALSE,error=FALSE,warning=FALSE}
edu_data <- gss %>% select(education, feelings_life) %>% filter(!is.na(feelings_life) & !is.na(education))
group <- edu_data %>%group_by(education) %>% summarise(n=n())
edudata <- edu_data %>%group_by(feelings_life) %>% summarise("high school" = sum(education == "High school diploma or a high school equivalency certificate"),"below high school" = sum(education == "Less than high school diploma or its equivalent"), "college" = sum(education == "College, CEGEP or other non-university certificate or di..."),"Bachelor" = sum(education == "Bachelor's degree (e.g. B.A., B.Sc., LL.B.)"), "trade certificate" = sum(education == "Trade certificate or diploma"),"under Bachelor" = sum(education == "University certificate or diploma below the bachelor's level"), "above Bachelor" = sum(education == "University certificate, diploma or degree above the bach...")) 
edudata$`high school` = edudata$`high school`/sum(edudata$`high school`)
edudata$`below high school` = edudata$`below high school`/sum(edudata$`below high school`)
edudata$college = edudata$college /sum(edudata$college)
edudata$Bachelor =edudata$Bachelor /sum(edudata$Bachelor)
edudata$`trade certificate` = edudata$`trade certificate`/ sum(edudata$`trade certificate`)
edudata$`under Bachelor` = edudata$`under Bachelor` /sum(edudata$`under Bachelor`)
edudata$`above Bachelor` =edudata$`above Bachelor` /sum(edudata$`above Bachelor`)
dfc <- gather(edudata, key = edu_group, value = percent_of_Respondent, c("below high school", "high school","trade certificate","college", "under Bachelor", "Bachelor","above Bachelor"))
gfb <- ggplot(dfc, aes(x=feelings_life, y = percent_of_Respondent, group = edu_group, colour = edu_group)) +
  geom_line() +
  theme_classic() + scale_x_continuous(breaks = a)+
  labs( title="different education group relevant to feelings_life")
```

```{r,echo=FALSE,error=FALSE,warning=FALSE}
health_data <- gss %>% select(self_rated_health, feelings_life) %>% filter(!is.na(feelings_life) & !is.na(self_rated_health)& self_rated_health != "Don't know")
group <- health_data %>%group_by(self_rated_health) %>% summarise(n=n())
healthdata <- health_data %>%group_by(feelings_life) %>% summarise("Excellent" = (sum(self_rated_health == "Excellent")),"Fair" = (sum(self_rated_health == "Fair")),"Good" = (sum(self_rated_health == "Good")),"Poor" = (sum(self_rated_health == "Poor")),"Very good" = (sum(self_rated_health == "Very good")))
healthdata$Excellent = healthdata$Excellent /sum(healthdata$Excellent)
healthdata$Fair = healthdata$Fair/sum(healthdata$Fair)
healthdata$Good = healthdata$Good/sum(healthdata$Good)
healthdata$Poor = healthdata$Poor/sum(healthdata$Poor)
healthdata$`Very good` = healthdata$`Very good`/ sum(healthdata$`Very good`)
dff <- gather(healthdata, key = health_group, value = percent_of_Respondent, c("Excellent", "Fair","Good", "Poor", "Very good"))
 gfc <- ggplot(dff, aes(x=feelings_life, y = percent_of_Respondent, group = health_group, colour = health_group)) +
  geom_line() +
  theme_classic() + scale_x_continuous(breaks = a)+
  labs( title="different health rated group relevant to feelings_life")
```

```{r,echo=FALSE,error=FALSE,warning=FALSE}
mental_health_data <- gss %>% select(self_rated_mental_health, feelings_life) %>% filter(!is.na(feelings_life) & !is.na(self_rated_mental_health)& self_rated_mental_health != "Don't know")
group <- mental_health_data %>%group_by(self_rated_mental_health) %>% summarise(n=n())
mentalhealth_data <- mental_health_data %>%group_by(feelings_life) %>% summarise("Excellent" = (sum(self_rated_mental_health == "Excellent")),"Fair" = (sum(self_rated_mental_health == "Fair")),"Good" = (sum(self_rated_mental_health == "Good")),"Poor" = (sum(self_rated_mental_health == "Poor")), "Very good"= (sum(self_rated_mental_health == "Very good")))
mentalhealth_data$Excellent = mentalhealth_data$Excellent/sum(mentalhealth_data$Excellent)
mentalhealth_data$Fair = mentalhealth_data$Fair /sum(mentalhealth_data$Fair)
mentalhealth_data$Good = mentalhealth_data$Good /sum(mentalhealth_data$Good)
mentalhealth_data$Poor = mentalhealth_data$Poor/sum(mentalhealth_data$Poor)
mentalhealth_data$`Very good` =mentalhealth_data$`Very good`/sum(mentalhealth_data$`Very good`)
dfh <- gather(mentalhealth_data, key = metal_group, value = percent_of_Respondent, c("Excellent", "Fair","Good", "Poor", "Very good"))
gfd <- ggplot(dfh, aes(x=feelings_life, y = percent_of_Respondent, group = metal_group, colour = metal_group)) +
  geom_line() +
  theme_classic() + scale_x_continuous(breaks = a)+
  labs( title="different mental_health group relevant to feelings_life")
```

```{r,echo=FALSE,error=FALSE,warning=FALSE}
work_data <- gss %>% select(average_hours_worked, feelings_life) %>% filter(!is.na(feelings_life) & !is.na(average_hours_worked)& average_hours_worked != "Don't know")
group <- work_data %>%group_by(average_hours_worked) %>% summarise(n=n())
workdata <- work_data %>%group_by(feelings_life) %>% summarise("0 hour" = (sum(average_hours_worked == "0 hour")),"0.1 to 29.9 hours" = (sum(average_hours_worked == "0.1 to 29.9 hours")),"30.0 to 40.0 hours" = (sum(average_hours_worked == "30.0 to 40.0 hours")),"40.1 to 50.0 hours" = (sum(average_hours_worked == "40.1 to 50.0 hours")),"50.1 hours and more" = (sum(average_hours_worked == "50.1 hours and more")))
workdata$`0 hour` = workdata$`0 hour`/sum(workdata$`0 hour`)
workdata$`0.1 to 29.9 hours` = workdata$`0.1 to 29.9 hours`/sum(workdata$`0.1 to 29.9 hours`)
workdata$`30.0 to 40.0 hours` = workdata$`30.0 to 40.0 hours`/ sum(workdata$`30.0 to 40.0 hours`)
workdata$`40.1 to 50.0 hours` = workdata$`40.1 to 50.0 hours`/sum(workdata$`40.1 to 50.0 hours`)
workdata$`50.1 hours and more` = workdata$`50.1 hours and more`/sum(workdata$`50.1 hours and more`)
 dff <- gather(workdata, key = work_group, value = percent_of_Respondent, c("0 hour", "0.1 to 29.9 hours","30.0 to 40.0 hours", "40.1 to 50.0 hours", "50.1 hours and more"))
gfe <- ggplot(dff, aes(x=feelings_life, y = percent_of_Respondent, group = work_group, colour = work_group)) +
  geom_line() +
  theme_classic() + scale_x_continuous(breaks = a)+
  labs( title="different number of hours worked per week group relevant to feelings_life")
```

```{r,echo=FALSE,error=FALSE,warning=FALSE}
marital_data <- gss %>% select(marital_status, feelings_life) %>% filter(!is.na(feelings_life) & !is.na(marital_status)& marital_status != "Don't know")
group <- marital_data %>%group_by(marital_status) %>% summarise(n=n())
maritaldata <- marital_data %>%group_by(feelings_life) %>% summarise("Divorced" = (sum(marital_status == "Divorced")),"Living common-law" = (sum(marital_status == "Living common-law")),"Married" = (sum(marital_status == "Married")),"Separated" = (sum(marital_status == "Separated")),"Single, never married" = (sum(marital_status == "Single, never married")), "Widowed" = (sum(marital_status == "Widowed")))
maritaldata$Divorced = maritaldata$Divorced/ sum(maritaldata$Divorced)
maritaldata$`Living common-law` = maritaldata$`Living common-law`/ sum(maritaldata$`Living common-law`)
maritaldata$Married = maritaldata$Married/ sum(maritaldata$Married)
maritaldata$Separated =maritaldata$Separated/sum(maritaldata$Separated)
maritaldata$`Single, never married` = maritaldata$`Single, never married` /sum(maritaldata$`Single, never married`)
maritaldata$Widowed = maritaldata$Widowed /sum(maritaldata$Widowed)
dfg <- gather(maritaldata, key = marital_group, value = percent_of_Respondent, c("Divorced", "Living common-law","Married", "Separated","Single, never married","Widowed"))
gff <- ggplot(dfg, aes(x=feelings_life, y = percent_of_Respondent, group = marital_group, colour = marital_group)) +
  geom_line() +
  theme_classic() + scale_x_continuous(breaks = a)+
  labs( title="different types marital status per week group relevant to feelings_life")
```
```{r,echo=FALSE, warning=FALSE, fig.height = 12, fig.width = 12}
gg_all <- (gfa + gfb) /
(gfc + gfd) /
(gfe + gff)  +    # Create grid of plots with title
  plot_annotation(title = "fig.4 six interested relevent groups") & 
  theme(plot.title = element_text(hjust = 0.5))
gg_all 
```
Figure.5:
different groups of working hours relevant to feelings_life
```{r,echo=FALSE,error=FALSE,warning=FALSE}
tmp2 <- gss %>% select(average_hours_worked, feelings_life, sex) %>% filter(!is.na(feelings_life) & !is.na(average_hours_worked) &!is.na(sex) & average_hours_worked != "Don't know")
tmp21<- tmp2 %>%filter(sex =="Male")
tmp22<- tmp2 %>%filter(sex =="Female")
tmp211 <- tmp21 %>%group_by(average_hours_worked) %>% summarise(feelings_life_male = mean(feelings_life))
tmp212 <-tmp211 %>% mutate(average_hour = case_when(average_hours_worked == "0 hour" ~ 0,
                                                    average_hours_worked == "0.1 to 29.9 hours" ~ 30,
                                                    average_hours_worked == "30.0 to 40.0 hours" ~ 40,
                                                    average_hours_worked == "40.1 to 50.0 hours" ~ 50,
                                                    average_hours_worked == "50.1 hours and more" ~ 60))
tmp221 <- tmp22 %>%group_by(average_hours_worked) %>% summarise(feelings_life_female = mean(feelings_life))
tmp222 <-tmp221 %>% mutate(average_hour = case_when(average_hours_worked == "0 hour" ~ 0,
                                                    average_hours_worked == "0.1 to 29.9 hours" ~ 30,
                                                    average_hours_worked == "30.0 to 40.0 hours" ~ 40,
                                                    average_hours_worked == "40.1 to 50.0 hours" ~ 50,
                                                    average_hours_worked == "50.1 hours and more" ~ 60))
tmp3 <- merge(tmp212,tmp222,by=c("average_hour"))
tmp4 <- gather(tmp3, key = feeling_of_different_sex, value = feeling_of_life, c("feelings_life_male","feelings_life_female"))
ggplot(tmp4, aes(x=average_hour, y = feeling_of_life, group = feeling_of_different_sex, colour = feeling_of_different_sex)) +
  geom_line() +
  theme_classic()+
  labs(x = "average_hour of work per week" ,title="fig.5: different groups of working hours relevant to feelings_life")
```
Figure.6
different education group relevant to feelings_life
```{r,echo=FALSE,error=FALSE,warning=FALSE}
t_mp2 <- gss %>% select(total_children, feelings_life, sex) %>% filter(!is.na(feelings_life) & !is.na(total_children) &!is.na(sex))
t_mp21<- t_mp2 %>%filter(sex =="Male")
t_mp22<- t_mp2 %>%filter(sex =="Female")
t_mp211 <- t_mp21 %>%group_by(total_children) %>% summarise(feelings_life_male = mean(feelings_life))

t_mp221 <- t_mp22 %>%group_by(total_children) %>% summarise(feelings_life_female = mean(feelings_life))
t_mp3 <- merge(t_mp211,t_mp221,by=c("total_children"))
t_mp4 <- gather(t_mp3, key = feeling_of_different_sex, value = feeling_of_life, c("feelings_life_male","feelings_life_female"))
ggplot(t_mp4, aes(x=total_children, y = feeling_of_life, group = feeling_of_different_sex, colour = feeling_of_different_sex)) +
  geom_line() +
  theme_classic() + scale_x_continuous(breaks = a)+ 
  labs( title="fig.6: different education group relevant to feelings_life")
```

Figure.7:
different groups of working hours relevant to feelings_life
```{r,echo=FALSE,error=FALSE,warning=FALSE}
age_diff_data <- gss %>% select(age_diff, feelings_life,sex) %>% filter(!is.na(feelings_life) & !is.na(age_diff)& age_diff != "Don't know" & sex == "Male")
tmp1 <- age_diff_data %>%group_by(age_diff) %>% summarise(feelings_life = mean(feelings_life))
tmp11 <-tmp1 %>% mutate(agediff = case_when(age_diff == "Respondent and spouse/partner are same age" ~ 0,
                                                    age_diff == "Respondent is 1 year older" ~ 1,
                                                    age_diff == "Respondent is 2 years older" ~ 2,
                                                    age_diff == "Respondent is 3 years older" ~ 3,
                                                    age_diff == "Respondent is 4 years older" ~ 4,
                                                    age_diff == "Respondent is 5 years older"~ 5,
                                                    age_diff == "Respondent is 6 to 10 years older" ~ 6,
                                                    age_diff == "Respondent is 11+ years older" ~ 7)) %>% select(agediff,feelings_life)%>% filter(!is.na(agediff))
tmp12 <-tmp1 %>% mutate(agediff = case_when(age_diff == "Respondent and spouse/partner are same age" ~ 0,
                                                    age_diff == "Respondent is 1 year younger" ~ 1,
                                                    age_diff == "Respondent is 2 years younger"~ 2,
                                                    age_diff == "Respondent is 3 years younger"~ 3,
                                                    age_diff == "Respondent is 4 years younger"~ 4,
                                                    age_diff == "Respondent is 5 years younger"~ 5,
                                                     age_diff == "Respondent is 6 to 10 years younger"~ 6,
                                                    age_diff == "Respondent is 11 + years younger"~ 7)) %>% select(agediff,feelings_life)%>% filter(!is.na(agediff))
tmp3 <- merge(tmp11,tmp12,by=c("agediff"))
colnames(tmp3)<- c("agediff","feelings_life older", "feelings_life younger")
tmp4 <- gather(tmp3, key = feeling_of_agediff, value = feeling_of_life, c("feelings_life older", "feelings_life younger"))
gfg <- ggplot(tmp4, aes(x=agediff, y = feeling_of_life, group = feeling_of_agediff, colour = feeling_of_agediff)) +
  geom_line() +
  theme_classic()+
  labs(x = "difference of age for male", title="different groups of working hours relevant to feelings_life")
```

```{r,echo=FALSE,error=FALSE,warning=FALSE,fig.height=6 , fig.width = 12}
age_diff_data <- gss %>% select(age_diff, feelings_life,sex) %>% filter(!is.na(feelings_life) & !is.na(age_diff)& age_diff != "Don't know" & sex == "Female")
tmp1 <- age_diff_data %>%group_by(age_diff) %>% summarise(feelings_life = mean(feelings_life))
tmp11 <-tmp1 %>% mutate(agediff = case_when(age_diff == "Respondent and spouse/partner are same age" ~ 0,
                                                    age_diff == "Respondent is 1 year older" ~ 1,
                                                    age_diff == "Respondent is 2 years older" ~ 2,
                                                    age_diff == "Respondent is 3 years older" ~ 3,
                                                    age_diff == "Respondent is 4 years older" ~ 4,
                                                    age_diff == "Respondent is 5 years older"~ 5,
                                                    age_diff == "Respondent is 6 to 10 years older" ~ 6,
                                                    age_diff == "Respondent is 11+ years older" ~ 7)) %>% select(agediff,feelings_life)%>% filter(!is.na(agediff))
tmp12 <-tmp1 %>% mutate(agediff = case_when(age_diff == "Respondent and spouse/partner are same age" ~ 0,
                                                    age_diff == "Respondent is 1 year younger" ~ 1,
                                                    age_diff == "Respondent is 2 years younger"~ 2,
                                                    age_diff == "Respondent is 3 years younger"~ 3,
                                                    age_diff == "Respondent is 4 years younger"~ 4,
                                                    age_diff == "Respondent is 5 years younger"~ 5,
                                                     age_diff == "Respondent is 6 to 10 years younger"~ 6,
                                                    age_diff == "Respondent is 11 + years younger"~ 7)) %>% select(agediff,feelings_life)%>% filter(!is.na(agediff))
tmp3 <- merge(tmp11,tmp12,by=c("agediff"))
colnames(tmp3)<- c("agediff","feelings_life older", "feelings_life younger")
tmp4 <- gather(tmp3, key = feeling_of_agediff, value = feeling_of_life, c("feelings_life older", "feelings_life younger"))
gfh <-ggplot(tmp4, aes(x=agediff, y = feeling_of_life, group = feeling_of_agediff, colour = feeling_of_agediff)) +
  geom_line() +
  theme_classic()+
  labs(x = "difference of age for female", title="fig.7: different groups of working hours relevant to feelings_life")
(gfg | gfh)
```
Figure.8:
different feelings_life group
```{r,echo=FALSE,error=FALSE,warning=FALSE}
work_last_data <- gss %>% select( worked_last_week, feelings_life) %>% filter(!is.na(feelings_life) & !is.na(worked_last_week)& worked_last_week != "Don't know")
work_last_data$feelings_life <-round(work_last_data$feelings_life)
group <- work_last_data %>%group_by(worked_last_week) %>% summarise(n=n())
 ggplot(work_last_data, aes(x= feelings_life, fill = worked_last_week)) +
  geom_bar() +
  labs( title="different feelings_life group")
```
Table 1
```{r}
table1 <- feel %>% summarise(min = min(feelings_life),
                             `1st Qu.` = quantile(feelings_life, 0.25),
                             median = median(feelings_life),
                             `3st Qu.` = quantile(feelings_life, 0.75),
                             max = max(feelings_life),
                             IQR = `3st Qu.` -`1st Qu.`,
                             sd = sd(feelings_life),
                             small_outliers = sum(feelings_life< `1st Qu.` -1.5*IQR),
                             large_outliers = sum(feelings_life> `3st Qu.` +1.5*IQR)) %>% 
  kable(caption = "general interested relevant groups") %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  kable_styling(latex_options = "HOLD_position")
```
Table 2
```{r,echo=FALSE, warning=FALSE}
table <-gss %>% filter(!is.na(income_respondent)&!is.na(education)&!is.na(self_rated_health)&!is.na(self_rated_mental_health)&!is.na(average_hours_worked)&!is.na(average_hours_worked)&!is.na(marital_status)& self_rated_health != "Don't know"& self_rated_mental_health != "Don't know"&average_hours_worked != "Don't know"& marital_status != "Don't know")
table1 <- table %>% count(income_respondent, sort = TRUE)
table2 <- table %>% count(education, sort = TRUE) 
table3 <- table %>% count(self_rated_health, sort = TRUE) 
table4 <-table %>% count(self_rated_mental_health, sort = TRUE) 
table5 <-table %>% count(average_hours_worked, sort = TRUE) 
table6 <-table %>% count(marital_status, sort = TRUE) 
list(table1,table6)%>% 
  kable(caption = "general interested relevant groups")%>%
  kable_classic(full_width = F) %>%
  kable_styling(latex_options = "HOLD_position") 
list(table3,table4,table5)%>% 
  kable() %>%
  kable_classic(full_width = F)%>%
  kable_styling(latex_options = "HOLD_position") 
table2%>% 
  kable() %>%
  kable_classic(full_width = F)%>%
  kable_styling(latex_options = "HOLD_position") 
```

