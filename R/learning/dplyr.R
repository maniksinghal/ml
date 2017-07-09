library(nycflights13)
library(tidyverse)
?flights
flights
head(flights)
tail(flights)
View(flights)
# What we'll learn
# filter()
# arrange()
# select()
# mutate()
# summarize()
# group_by()

################## Filter ########################
# flights that started on jan 1
filter(flights,month==1,day==1)
fstjan <- filter(flights,month==1,day==1)
(fstjan <- filter(flights,month==1,day==1))

#flights that departed in Nov / Dec
filter(flights, month == 11 | month == 12)
nov_dec <- filter(flights, month %in% c(11, 12))

#flights thatweren't delayed (on arrival or departure)
#by more than two hours

filter(flights, !(arr_delay > 120 | dep_delay > 120))
filter(flights, arr_delay <= 120, dep_delay <= 120)

################## Arrange #######################
arrange(flights, year, month, day)
arrange(flights, desc(arr_delay))
#most delayed
arrange(flights, desc(dep_delay))
#least delay / departed early
arrange(flights, dep_delay)
# travelled longest and shortest
arrange(flights, desc(distance)) %>% 
  select(1:5, distance)
arrange(flights, distance) %>%
  select(1:5, distance)
################## Select ########################
select(flights, year, month, day)
select(flights, year:day)
select(flights, -(year:day))
select(flights, tail_num = tailnum)
rename(flights, tail_num = tailnum)
select(flights, time_hour, air_time, everything())
################## mutate ########################
flights_sml <- select(flights,  year:day, 
                      ends_with("delay"), 
                      distance,  air_time)

mutate(flights_sml,  
       gain = arr_delay - dep_delay, 
       speed = distance / air_time * 60)


mutate(flights_sml,  
       gain = arr_delay - dep_delay, 
       hours = air_time / 60,  
       gain_per_hour = gain / hours)

transmute(flights, 
          gain = arr_delay - dep_delay, 
          hours = air_time / 60,  
          gain_per_hour = gain / hours)

################## summarize ########################
summarize(flights, 
          delay = mean(dep_delay, na.rm = TRUE))

by_day <- group_by(flights, year, month, day)

by_day

summarize(by_day, 
          delay = mean(dep_delay, na.rm = TRUE))

flights %>% 
  group_by(year,month,day) %>% 
  summarise(delay = mean(dep_delay, na.rm = TRUE))


# relationship between distance and average delay
# for each location
(by_dest <- group_by(flights, dest))
delay <- summarize(by_dest,
  count = n(),
  dist = mean(distance, na.rm = TRUE),
  delay = mean(arr_delay, na.rm = TRUE)
)
delay <- filter(delay, count > 20, dest != "HNL")
delay


delays <- flights %>%
  group_by(dest) %>%
  summarize(
    count = n(),
    dist = mean(distance, na.rm = TRUE),
    delay = mean(arr_delay, na.rm = TRUE)
  ) %>%
  filter(count > 20, dest != "HNL")

not_cancelled <- flights %>%
  filter(!is.na(dep_delay), !is.na(arr_delay))

not_cancelled %>%
  group_by(year, month, day) %>%
  summarize(mean = mean(dep_delay))

delays <- not_cancelled %>%
  group_by(tailnum) %>%
  summarize(
    delay = mean(arr_delay)
  )

ggplot(data = delays, mapping = aes(x = delay)) +
  geom_freqpoly(binwidth = 10)


