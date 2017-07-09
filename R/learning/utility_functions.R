seq(start,stop,factor)  # gen sequence from start to stop adding factor 
rep(a, times=b)  # repeat a, b times
rep(a, each=b)  # repeat each element of a, b ti

round(a) # Round to nearest integer (3.7 => 4, 3.3 => 3) (Optinal ,n => round to that many decimals)
ceiling(a)  # If a has decimal, then smallest integer greater than a => a+1
floor(a)    # If a has decimal, largest integer less than a floor(-3.5) => -4
trunc(a)    # Strips of decimal
signif(a,n) #Rounds a to specified signif digits. signif(123.25,2) => 120

#> find.package("ggplot2")
#[1] "/usr/local/lib/R/3.4/site-library/ggplot2"
#> .libPaths()
#[1] "/usr/local/lib/R/3.4/site-library"                                     
#[2] "/usr/local/Cellar/r/3.4.0_1/R.framework/Versions/3.4/Resources/library"

find.package("ggplog2")
dir("path")


# Functions
#Annonymous function:
(function(x){3*x+2})(1)   # fn: 3x+2, input:1
(function(x){3*x+2})(c(-10:10))  # Now passing vector. Output is also vector

#Named function
fun1 <- function(a, b)
{
  print("This is a function")   # Only characters possible. By default adds new line.
  cat("a*b is",a*b,"\n")   # To print characters and numbers. No new-line by default
  print("hello")
}

#Variable arguments
fun3 <- function(a, b=2)     # fun3(1), fun3(1,) fun3(1,2), fun3(1,b=100) => all accepted
{                            # Takes upto 2 args. errors if passed more arguments.
  print(a+b)
}

#Ignore additional arguments
fun4 <- function(a=10, b=15,...)  #Ignore any additional arguments passed using ...
{
  #Global variable (using <<)
  my_glob <<- a*b
  print(a*b)
}

fun5 <- function(...)
{
  print(my_glob)   # Its a global declared in fun4 so can be refrenced.
  print(c(...))
}

#Passing function to function
fun6 <- function(df, fun=mean,...)  # If user passes na.rm=T (remove NAs, it gets passed to fun)
{
  for(i in 1:length(df)) {
    print(fun(df[[i]],...))
  }
}


lapply(df,mean,na.rm=T)
# You can pass any object to lapply (vector, list, matrix, DF)
# Apply function on each element
# Always returns list output
lapply(list(friends=10,degrees=c("BE","PGDM","FPM"),married=TRUE),nchar)
# Output: list => 2, (2 4 3), 4


#sapply (Similar to lapply, but returns a better formatted output)
sapply(df,mean,na.rm=T) #simplify apply
# Go through each element of df [[]], and gives it to mean
# Output format depnds on input
# If output for function on each element is not same length, then returns a list
# If returning a list, it returns a named list

# Similar to sapply, but takes additional argument, the output format
vapply(df,mean)

# Given a package name
# Find all packages in the directory where the package is located
ls_pkg_dir <- function(path)
{
  split <- strsplit(find.package(path), "/")[[1]]  # String-split using delimiter /
  #Remove the package from path
  split <- split[-length(split)]  # Remove last element from vector (the package name)
  dir(paste(split,collapse="/"))  # Collapse vector again to form directory path, and dir to it
}

# In IRIS object, find mean of Sepal.Length for each species
# separately.
iris_avg <- function()
{
  iris <- iris
  for(i in 1:length(unique(iris$Species)))   # Iterate over unique species
  {
    cat(as.character(unique(iris$Species)[i]))  # Print species name
    
    # Filter records having that specie
    # Find the Sepal.Length
    # Mean it
    print(mean(iris[iris$Species == unique(iris$Species)[i],]$Sepal.Length))
  }
  
  #Other way
  # tapply: Split - Apply - Combine
  # View Sepal.Length for every Species. Don't show everything, just the mean
  # tapply returns an array
  tapply(iris$Sepal.Length, iris$Species,mean)
  
  # Yet another way
  # aggregate: ~ represents a formula
  # Show Sepal.length as a function os Species, and aggregate as mean
  # Returns data frame
  aggregate(Sepal.Length~Species,data=iris,mean)
  
}

#Find mileage per gear and mean it
tapply(mtcars$mpg,mtcars$gear,mean)
aggregate(mpg~gear,data=mtcars,mean)

#Lets say we want to add the mean mileage now
sum(tapply(mtcars$mpg,mtcars$gear,mean))
sum(aggregate(mpg~gear,data=mtcars,mean))
# Note that sum, if passed a DF, unwinds it and then adds all elements
# So gears also got added when performing sum(aggregate..)
sum(aggregate(mpg~gear,data=mtcars,mean)$mpg)  # Correct way

library(nycflights13)
library(tidyverse)
flights
# dplyr package:
#==============
#First argument is always DF
# Output is also DF

#Filter is always row operation.
filter(flights,month==1,day==1)  # Filter records with month==1 and day==1 (, is AND, | is OR)
(fstjan <- filter(flights, month==1,day==1))    #( ) around to print the output as well

filter(flights,month >= 11,month <= 12)  # Flights departing in Nov or Dec
filter(flights, month == 11| month == 12)
filter(flights, month %in% c(11,12))

arrange(flights, year, month, day)   # Order rows by given columns 
arrange(flights, desc(arr_delay))
arrange(flights, arr_time)

arrange(flights, distance) %>%    
  select(arr_time, distance)       # Select only arr_time, distance columns
select(flights, year, month, day)
select(flights, year:day)  # All cols from year to day
select(flights, -(year:day))  # All cols except from year to day

select(flights, tail_num = tailnum)   # Select tailnum column and label as tail_num
rename(flights, YEAR = year)

select(flights, time_hour, air_time, everything())  # Show time_hour/air_time first, then reset
select(flights, year:day, ends_with("delay"), distance, air_time)
  # Select year to day, all cols that end with delay, then distance, air_time

flights_sml <- select(flights, year:day, ends_with("delay"), distance, air_time)

mutate(flights_sml, gain = arr_delay - dep_delay,
       speed = distance / air_time * 60)    # Create new variable and add more columns

transmute(flights, gain = arr_delay - dep_delay,  # Creates with only 3 columns, not rest of them.
            hours = air_time / 60,
            gain_per_hour = gain / hours)

summarize(flights, delay=mean(dep_delay, na.rm=TRUE))
View(flights)        #Excel form of data set



(by_day <- group_by(flights, year, month, day))
summarize(by_day, delay = mean(dep_delay, na.rm = TRUE))  # day is gone
# one more summarize: month would be gone, right to left

(by_dest <- group_by(flights, dest))
delay <- summarize(by_dest, count = n(), dist = mean(distance, na.rm=TRUE), delay = mean(arr_delay, na.rm=TRUE))

not_cancelled <- flights %>% filter(!is.na(dep_delay), !is.na(arr_delay))
not_cancelled
