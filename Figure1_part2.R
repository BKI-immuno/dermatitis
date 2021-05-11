# antigen driven clustering 

library(tidyverse)
library(data.table)

all.tcr.n.cli.n<-read.csv("data.analysis.csv")

# down stream analysis for GLIPH output
gliph.derm<-read.csv("dermatitis.gliph.output.csv")
test<-gliph.derm %>% group_by(pattern,Sample,vb_score,length_score) %>% summarise(n=n(),freq=sum(Freq)) %>%
  filter(Sample=='10:LD'&vb_score<0.05&length_score<0.05&freq>10) %>% arrange(-freq)

gliph.derm$lable<-sub('.*:','',gliph.derm$Sample)

colnames(gliph.derm)[14]<-'CDR3.amino.acid.sequence'

gliph.derm.prop<-gliph.derm %>% left_join(all.tcr.n.cli.n[,c('lable','CDR3.amino.acid.sequence','productive')])
clustered<-gliph.derm.prop %>% filter(pattern=='SSQD') %>% arrange(-productive)


# cell proportion of clones with ssqd
clustered %>% group_by(lable)%>% summarise(sum=sum(productive)) %>%
  ggplot(aes(x=reorder(lable,-sum),y=100*sum)) +
  geom_bar(stat='identity')+
  theme_classic()+
  labs(x='',y='% among all T cells')
ggsave("ssqd.by.sample.tiff",width=6,height=4,dpi=300)



