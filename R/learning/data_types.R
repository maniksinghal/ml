#Data types
int1 <- 123L  # integer
class(int1)   # Tell class of the data type
val1 <- 25.12345  # numeric (default type)
c1 <- "learning R"  # character type  (Cannot be typecasted to integer)
c2 <- "l"
class(as.integer(c1))   # We also have mode() and typeof()

# Class precedence:
# => character => numeric => integer => logical

# Typecast one type to another
as.integer(val1)  
a <- as.integer("123")  #works fine, if character is numerical form
c <- as.character(200L)  # integer gets stringified
l <- as.character(TRUE)  # converts to "TRUE"
p <- as.logical("FALSE")  # possible (converts to FALSE)
is.integer(val1)  # validate class of the data type

#Indeterminate (NA) vs NULL
# => Cannot be detemined (typecasting character to integer)
c1 <- "some character"
i <- as.integer(c1)  # Cannot be converted (NA)
l <- as.logical(35)

#Factor variables
# Assign dummy values to character variables
gender <- c("male", "female", "female", "male")
genderF <- factor(gender)  # Assign values in alphabetical order
as.character(genderF)  # Gives back the character array

iq <- c("low", "high", "veryhigh", "medium")
iqF <- factor(iq)  # Values get assigned randomly (based on alphabetic order)
      # high = 1, low = 2, medium = 3, veryhigh =4
str(iqF)  #Factor w/ 4 levels "high","low","medium",..: 2 1 4 3

iqF <- factor(iq, ordered = T, levels=c("low", "medium", "high", "veryhigh"))
str(iqF)  #Ord.factor w/ 4 levels "low"<"medium"<..: 1 3 4 2

#Factoring integers
f = factor(c(1:10))
as.integer(f)  # Converts back
as.logical(f)  # Factor has unary values, cannot be converted to binary (NA NA NA ...)

#Vectors
gender <- c("male", "female", "male", "female")  #concatenate
class(gender)  #Tells data-type of element in vector
l = c(1:100)  # Vector of 100 elements
mixvec <-  c(1,2,3,T,'VISTA')  # Converted to strongest class -> character vector
#OR
mixvec <- vector(mode="logical",length=2)  # FALSE FALSE
mixvec <- vector(mode="character",length=2) # "" ""
mixvec <- vector(length=1)  # Default is logical
#Indexing:
numvec[2] # 2nd element of numvec
numvec[]  # all elements of numvec
intvec[1,2]   # Trying to reference 2D subset, but vector is 1D
intvec[c(1,2)]  #
intvec[1:2]   #Values from 1 to 2
intvec[c(T,F)]  # Omit 2nd/4th/6th.. element and return others. (T,F recycles T,F,T,F,T,F)
intvec[c(1,3)]   # Returns elements 1 and 3
intvec[-3]  # Returns vector with 3rd element removed
intvec * numvec  # Element wise multiply (dimensions should be same)
intvec * 2   # Elementwise multiply by 2
intvec == intvec2  # Comparison of vectors element wise (Returns vector of T,F)
rev(intvec)  # reverse 
append(intvec, numvec)  # numvec appended to intvec
length(intvec)  # Size of vector (num elements)
nchar(charvec)  # Vector of sizes of individual elements (number of digits for integer array)

intvec[1000]  # NA if actual size is small. Internally creates remaining elements, all NA
length(intvec) # Now becomes 1000, with all remaining elements set to NA

#named elements
names(intvec)  # by default returns NULL
names(intvec) <- c("guy1", "guy2", "guy3")

intvec = c(100, 101, 102, 103)
names(intvec) <- c("EL0", "EL1", "EL2", "EL3", "EL4")
intvec
#EL0 EL1 EL2 EL3 EL4 
#100 101 102 103 104
intvec["EL2"]   # same as intvec[2]

intvec <- c(5,3,6,8,10,2)
sort(intvec)  # Sorted vector    2 3 5 6 8 10
order(intvec)  # Position of sorted elements in original vector  6 2 1 3 4 5
rank(intvec) # Rank of each element in the list (OR position of element in sorted order)

table(intvec)
# Returns a named vector
# Names: Unique (sorted) elements of intvec
# Values: Frequency of the element in the intvec

unique(intvec)  # Unique values from the vector, in order of occurrence
match(intvec, u)  # Return position of each element in intvec in vector u



#List/Remove variables in memory


ls()   # Shows all variabels in a list
rm(var) # Remove/delete variable from memory
rm(list=ls())  #Remove all variables from memory


#20 %% 3 => remainder
# 20 %/% 3 => quotient

#Dates
#======
a <- Sys.Date()     # Stored as num days since Jan 1, 1970

datetext<-"2017-07-03"
d1 <- as.Date(datetext)

as.Date('12/26/2016', format="%m/%d/%Y")
#Day Options are %d only
# %d -	 Day of the month

#Month Options are %m,%b and %B
# %m -   Month (integer)
# %b -   Month (abbreviated)
# %B -   Month (full name)

#Year options are %y and %Y
# %y -	 Year (2 digits)
# %Y -	 Year (4 digits)

date2==date1  # Compare dates







