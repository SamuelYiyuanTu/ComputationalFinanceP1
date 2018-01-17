# 405 HW

# Question 2 histogram
q2_data = read.table("/Users/samueltu/Desktop/Winter_2018/405_ComputFin/Week1/q2.csv")
hist(q2_data[,1])

# Question 3 histogram
q3_data = read.table("/Users/samueltu/Desktop/Winter_2018/405_ComputFin/Week1/q3.csv")
hist(q3_data[,1])
# Theoretical Value
sum(dbinom(x = 40:44, size = 44, prob = 0.64))

q4_data = read.table("/Users/samueltu/Desktop/Winter_2018/405_ComputFin/Week1/q4.csv")
hist(q4_data[,1])
