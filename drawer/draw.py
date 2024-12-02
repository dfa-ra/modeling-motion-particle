import matplotlib.pyplot as plt

print("start drawer")

with open("/media/ra/_work/ra/ITMO/PHISICS/projectC/drawer/data.txt", "r") as file:
    coordinates = [tuple(map(float, line.split())) for line in file]

print(coordinates)

x_values = [coordinate[0] for coordinate in coordinates]
y_values = [coordinate[1] for coordinate in coordinates]
z_values = [coordinate[2] for coordinate in coordinates]


plt.figure(figsize=(10, 10))
plt.plot(x_values, y_values,  linestyle='-')
# plt.xlim(-0.0000002, 0.0000002)
# plt.ylim(-0.0000001, 0.000001)

# plt.xlim(100000, 100000)
# plt.ylim(100000, 100000)

plt.title("График по точкам")
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.axhline(0, color='black', linewidth=0.5, ls="--")
plt.axvline(0, color='black', linewidth=0.5, ls="--")

plt.show()