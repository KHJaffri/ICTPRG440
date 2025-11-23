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

def haversine_distance(lat1, lon1, lat2, lon2):
    # Radius of Earth in kilometers
    R = 6371.0

    # Convert degrees to radians
    lat1 = math.radians(lat1)
    lon1 = math.radians(lon1)
    lat2 = math.radians(lat2)
    lon2 = math.radians(lon2)

    # Differences
    dlat = lat2 - lat1
    dlon = lon2 - lon1

    # Haversine formula
    a = math.sin(dlat / 2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    # Distance in km
    distance = R * c
    return distance

# Sydney and Melbourne coordinates
lat_syd, lon_syd = 33.8688, 151.2093
lat_melb, lon_melb = 37.8136, 144.9631

# Calculate distance
distance_km = haversine_distance(lat_syd, lon_syd, lat_melb, lon_melb)

print("Haversine Distance between Sydney and Melbourne:", round(distance_km, 2), "km")