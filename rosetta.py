from rosetta import rosetta, SoilData

data = [
        [30,30,40,1.5,0.3,0.1],
        [20,60,20],
        [55,25,20,1.1]
    ]

mean, stdev, codes = rosetta(3, SoilData.from_array(data))
print(mean)