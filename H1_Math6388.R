FlyAtlas <- read_tsv("http://nsmn1.uh.edu/rpmeisel/teaching/20090519all.txt")
colnames(FlyAtlas) <- make.names(colnames(FlyAtlas))
FlyAtlas

FlyAtlas %>%
  arrange(desc(BrainMean))

FlyAtlas %>%
  group_by(Oligo) %>%
  arrange(desc(BrainMean))
filter(FlyAtlas, BrainMean>3619)

FlyAtlas
FlyAtlas$Oligo

FlyAtlas %>% 
  summary()


##Number 7
FlyAtlas %>%
  arrange(desc(BrainMean)) %>%
  select(Oligo)

mean(c(FlyAtlas$BrainMean,FlyAtlas$HeadMean))

###Number 8 

FlyAtlas %>% mutate(HeadBrain = (BrainMean + HeadMean)/2) %>%
  select(Oligo, HeadBrain, HeadMean, BrainMean)

Atlas1 <- FlyAtlas %>% mutate(HeadBrain1 = BrainMean * n() + HeadMean * n()) %>% 
  mutate(HeadBrain1/n()) %>% select(Oligo, HeadBrain1, HeadMean, BrainMean)

Atlas1  

##Number 9
FlyAtlas %>%
  ggplot(aes(x = HeadMean, y = BrainMean)) + geom_point()

###Number 10
##Yes 
FlyAtlas %>%
  drop_na(HeadMean,BrainMean) %>%
  ggplot(aes(x = HeadMean, y = BrainMean)) + geom_point()


###Number 11
FlyAtlas %>%
  drop_na(HeadMean,BrainMean, head.vs.whole.fly....T.Test_Change.Direction) %>%
  ggplot(aes(x = HeadMean, y = BrainMean, color = head.vs.whole.fly....T.Test_Change.Direction)) + geom_point()

###Number 12

FlyAtlas %>%
  drop_na(HeadMean,BrainMean, head.vs.whole.fly....T.Test_Change.Direction) %>%
  ggplot(aes(x = HeadMean, y = BrainMean, color = head.vs.whole.fly....T.Test_Change.Direction)) + geom_point() + labs(color='Head v Fly')

### Number 13

FlyAtlas %>% 
  drop_na(HeadMean,BrainMean,head.vs.whole.fly....T.Test_Change.Direction,brain.vs.whole.fly...T.Test_Change.Direction) %>%
  ggplot(aes(x = HeadMean, y = BrainMean, color = head.vs.whole.fly....T.Test_Change.Direction)) + 
  geom_point() + facet_grid(~brain.vs.whole.fly...T.Test_Change.Direction) + labs(color='Head v Fly')

### Number 14
FlyAtlas %>%
  summarise(meanHead = mean(HeadMean, na.rm=TRUE))

## Number 15
FlyAtlas %>% 
  filter(head.vs.whole.fly....T.Test_Change.Direction == "Up") %>%
  summarise(meanHeadUP = mean(HeadMean, na.rm = TRUE))

### Number 16
FlyAtlas %>%
  group_by(head.vs.whole.fly....T.Test_Change.Direction) %>%
  summarise(meanBrain = mean(BrainMean))



