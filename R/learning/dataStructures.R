#Data-structures
#                     1-dimensional     2-dimensional         > 2D
# Homogeneous         Vector(1D array)   Matrix(2D array)     Array(3D/4D... array)

# Heterogeneous       List               Data frame (table of 
#                                         different columns)

# List is one dimensional in nature, but it can contain any objects
# a list, a matrix, ... etc.
l = list(c(1:3), "manik", c(20:40))
str(l)
#List of 3
#$ : int [1:3] 1 2 3
#$ : chr "manik"
#$ : int [1:21] 20 21 22 23 24 25 26 27 28 29 ...

#get
# Converts an object name to object itself
for (i in ls()) {
  print(get(i))  #i => obj-name, get(i) => obj
}



#Matrices
#========
matrix(r,c)
dim(matrix(1,3))  # Vector of num-rows, num-colums
t(x)   # Transpose
matrix(0,3,4)  # 3x4 matrix filled with 0s
mat1 <- matrix(vec, nrow=3)  # Re-shape vector into matrix. Values filled column-wise.
                             # If not exact multiple, filling starts over from begninning
matrix(vec, nrow=4,byrow=T)  #fill row-wise
matrix(vec, nrow=3,ncol=4)  # Use only 3x4 elements, discard rest while reshaping


intvec_10elem = c(1:10)
dim(intvec_10elem)  # is NULL, only 1D
dim(intvec_10elem) <- c(5,2)  # Reshape vec to 5 rows 2 columns
dim(intvec_10elem) <- c(3,4)  # Gives error, dimension mismatch. No recycling as with matrix()
                    <- c(3,3)  # Again gives error, dimension should exectly match.
class(mat)  # Returns 'matrix'
typeof(mat)  # integer

cbind(vec_a, vec_b) # attach vectors as columns into a matrix

rownames(mat) # Row names vector
colnames(mat) <- c("a", "b", ...)  # 
dimnames(mat) <- list(c("r1", "r2", "r3"), c("c1", "c2"))  #Assign both row and col names

mat[i]  # ith element of matrix, when parsed column-wise
mat[i,] # ith row
mat[,i] # ith column

mat_10elem[20]  # Vector view, so gives NA
mat_10elem[20,20]  # Does not exist, error.


# Lists
list1 <- vector(mode="list",length=3)
length(list1)  # 3
list1 <- list()  # Gets empty list

names(list1) # Names of elements
names(list1) <- c("one", "two", "three")

#Every element of a list is a list
list1 <- list(c(1,2), T, c("a", "b", "c"))  # 3 elements, 1st and 3rd being lists themselves.
list1[1]  # We get a list object of one element
list1[1][1][1][1][1]... # same as list1[1]
list[[1]] # This gives the actual element in the first position of the list
#$one
#[1] 1 2

list[[1]] # This gives the inner vector
# [1] 1 2

list1$one  # Called by name => equivalent to [[]]. Gives the inner element, not as a list
list[1:2]  # Elements 1 to 2 as a separate list
unlist(list1)  # Recursively unwind a list into a single vector of elements


# Data frames
# When printed, shows row-numbers and column-names




#==============
name1 = c("Antony", "Ben", "Chthy", "Diana")
age <- c(12, 18, 17, 19)
gender <- factor(c("male", "female", "male", "female"))  # Factor relevant vectors yourself
df <- data.frame(name1, age, gender, stringAsFactors=F)
df[]  # Everything
df[1]  # First column as a data frame (with row numbers and col name)
df["name1"]  # Same thing
df[1,]   # First row all columns as data frame (Mixed data type) - as data frame
df[,1]  #  Vector (First column, homogenous)
df[[1]]  # Vector
df$name1  #vector
df[1,2]  # Age(12), vector
head(df,n)  # First 6 (by default) or n rows of DF
tail(df)  # Last 6 values
summary(df)  # Prints summary (Min/max/median/mean/.. of each row)
