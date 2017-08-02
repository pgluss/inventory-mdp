# R for heatmaps -- Written by Paula Gluss, Dept. of Mathematics, Univ. of Wash.
# This code takes in a .csv file and outputs a heatmap with or without the optimal
# values at each state printed at their respective (x, y) coordinates.

# Supply the filename and full path
filename = "~/test1.csv"                        # Linux example
#filename = "C:\Users\pgluss\Desktop\ha3.csv" # Windows example

# Set desired bounds for the heatmap
# If you want the whole data set, set boundx to c(min(x),max(x)) and boundy to c(min(y),max(y))
boundx = c(0,10) 
boundy = c(0,10)

# Read in the data from the file and save x, y, f-vals, decisions into separate variables
data = read.csv(filename, sep = ",", header = T)
x = data[,1]                # Save x data (x data should be the first column in the CSV)
y = data[,2]                # Save y data (y data should be the second column in the CSV)
f = round(data[,3]*10) / 10 # Save and round f-values to one decimal point (should be third column of CSV)
dec = data[,4]              # Save decisions (decisions should be the fourth column of the CSV)

# Create a decision matrix where data will be stored
# A decision matrix is a two dimensional matrix that stores the decision for a given (x,y)
rows = max(x) - min(x) + 1  # Number of x coordinates
cols = max(y) - min(y) + 1  # Number of y coordinates
decM = matrix(nrow = rows, ncol = cols)  # Constructs the matrix

# Fill the matrix with decisions
# Iterates through each (x,y) combination and stores the associated decision in the matrix at
# the appropriate indices. Note that matrix indices start at 1 while (x,y) coordinates can be
# positive or negative. Thus, the following code takes that into account and adjusts as needed.
for (i in c(1:length(x))) {
  decM[x[i] + (1 - min(x)), y[i] + (1 - min(y))] = dec[i]
}

# Alter the size of the decision matrix to obtain the section you want to graph
# If you want the whole matrix, comment out the following line
decM = decM[(boundx[1] + 1 - min(x)) : (boundx[2] + 1 - min(x)), (boundy[1] + 1 - min(y)) : (boundy[2] + 1 - min(y))]

# Build the heatmap
xhm = c(1:nrow(decM))
yhm = c(1:ncol(decM))

image(xhm, yhm, decM, ylab = "Used Item Inventory", xlab = "New Item Inventory", xaxt = "n", yaxt = "n")
title(main = "Control Policy Heatmap")
legend('topright', c('Idle - Red', 'Produce from Raw - Orange', 'Produce from Used - Biege'), text.col = c('black', 'black', 'black'), bty= 'n')

# Add axis labels
axis(1, at = xhm, labels = c(boundx[1]:boundx[2]))
axis(2, at = yhm, labels = c(boundy[1]:boundy[2]))

# Add f-values to heat map. Comment out the whole for loop if the f-vals are not desired.
# To adjust the font size of the f-vals on the heatmap, change cex = XXX.
# cex = 0.1 is small, cex = 0.8 is large
for (i in c(1:length(f))) {
  text(x[i] + (1 - min(x)), y[i] + (1 - min(y)), labels = f[i], cex = 0.5)
}