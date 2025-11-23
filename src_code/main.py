def add2numbers(x, y):
    return x+y
print("KAshanHyderJaffri")
print("TAFENSW")
print( "KashanHyderJaffri")
print("TAFENSW")
print("Task1completed")
def myfunction (x,y,z):
    s=x+y+z
    return s
d = myfunction (1,2,3)
import math

# Coordinates
lat_syd = 33.8688
lon_syd = 151.2093
lat_melb = 37.8136
lon_melb = 144.9631

# Convert degrees to radians
lat_syd = math.radians(lat_syd)
lon_syd = math.radians(lon_syd)
lat_melb = math.radians(lat_melb)
lon_melb = math.radians(lon_melb)

# Haversine formula
dlat = lat_melb - lat_syd
dlon = lon_melb - lon_syd

a = math.sin(dlat / 2)**2 + math.cos(lat_syd) * math.cos(lat_melb) * math.sin(dlon / 2)**2
c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

# Radius of Earth in km (average)
R = 6371.0

# Distance
distance_km = R * c

print("Air distance between Sydney and Melbourne:", round(distance_km, 2), "km")