using Plots
using DelimitedFiles
temp = readdlm("./Output/MachDvals.txt")
dump(temp)


x = range(0, 1, 31)

plot(x[1:31], temp[1,:])

