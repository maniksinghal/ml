install.packages("tidyverse")
library(tidyverse)
# Do smaller cars give better mileage?
?mpg
head(mpg)
class(mpg)
str(mpg)
dplyr::glimpse(mpg)
mpg
# Every plot has 3 components 
# data + aesthetic mappings + layer for rendering observations
ggplot(data=mpg, aes(x=displ,y=hwy))+
  geom_point()

ggplot(data=mpg, aes(x=displ,y=hwy, color = class))+
  geom_point()

# Blue ignored here
ggplot(data=mpg, aes(x=displ,y=hwy, color = "blue"))+
  geom_point()

# add color in geom_point
ggplot(data=mpg, aes(x=displ,y=hwy))+
  geom_point(color="blue")

# Color coding based on drv
ggplot(data=mpg, aes(x=displ,y=hwy, color = drv))+
  geom_point()

# Visualize more columns using shape attribute (can specify size too)
ggplot(data=mpg, aes(x=displ,y=hwy, color = class, shape = drv))+
  geom_point()

ggplot(mpg, aes(cty, hwy, size = displ)) + 
  geom_jitter()      # Give a shock to the points

table(mpg$class)

# Relationship between city and hwy
ggplot(mpg,aes(cty,hwy))+
  geom_point()


# model and manufacturer
unique(mpg$model)
unique(mpg$manufacturer)
# Not a good way to view models per manufacturer
ggplot(mpg,aes(model,manufacturer))+geom_point()

df <- mpg %>% 
  mutate("manuf_model"=paste(manufacturer,model,sep =" "))
str(df)

sort(table(df$manuf_model))
ggplot(df,aes(manuf_model))+
  geom_bar()+
  coord_flip()


# Take a guess on what you'd see
ggplot(mpg, aes(cty, hwy)) + geom_point()
ggplot(diamonds, aes(carat, price)) + geom_point()
ggplot(economics, aes(date, unemploy)) + geom_line()
ggplot(mpg, aes(cty)) + geom_histogram()

#Drive Vs City Economy? 


ggplot(mpg,aes(drv,cty))+
  geom_boxplot()+
  scale_x_discrete(limits = c("f","r","4"),
  labels=c("Front Wheel","Rear Wheel","Four Wheel"))


ggplot(mpg,aes(drv,cty))+
  geom_violin()+
  scale_x_discrete(limits = c("f","r","4"),
  labels=c("Front Wheel","Rear Wheel","Four Wheel"))

# Class Vs hwy?

ggplot(mpg, aes(reorder(class, hwy, FUN = median), hwy)) +
  geom_boxplot()

#Drive Vs engine size and class?

unique(mpg$drv)
unique(mpg$class)
unique(mpg$displ)

str(mpg)
sapply(c(mpg["drv"],mpg["class"],mpg["displ"]),unique)

sapply(mpg[c("drv","class","displ")],unique)

keyvars <- c("drv","class","displ")
sapply(mpg[keyvars],unique)

ggplot(mpg,aes(class,displ,color=drv))+
  geom_point()

ggplot(mpg,aes(class,displ,color=drv))+
  geom_jitter(width=0.05,height=0)

ggplot(mpg, aes(reorder(class, displ, FUN = median), 
                displ, colour = drv)) + 
  geom_jitter(width = 0.15, height = 0)

ggplot(mpg, aes(reorder(class, -displ, FUN = median), 
                displ, colour = drv)) + 
  geom_jitter(width = 0.15, height = 0)

# faceting 

ggplot(mpg, aes(displ, hwy)) +
  geom_point() +
  facet_wrap(~class)

# trends

ggplot(mpg, aes(displ, hwy)) +
  geom_point() +
  geom_smooth()

ggplot(mpg, aes(displ, hwy)) +
  geom_point() +
  geom_smooth(method="lm")
