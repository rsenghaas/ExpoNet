import numpy as np
import matplotlib.pyplot as plt


class pathmap():
    def __init__(self, x_limits, y_limits, resolution):   
        self.x_lim = x_limits
        self.y_lim = y_limits
        self.resolution = resolution
        self.x_range = self.x_lim[1] - self.x_lim[0]
        self.y_range = self.y_lim[1] - self.y_lim[0]
        self.data = self.__initialize_map()

    # The data array is an array of dynamic lists.
    # The data saved is a pair of numbers (step indices that hit a pixel)
    def __initialize_map(self):
        data = np.empty((self.resolution[0], self.resolution[1]), dtype=object)
        for i in range(self.resolution[0]):
            for j in range(self.resolution[1]):
                data[i,j] = []
        return data

    def number_to_coordinates(self, z):
        x = z.real
        y = z.imag
        m = int((x - self.x_lim[0])/(self.x_range) * self.resolution[0])
        n = int((y - self.y_lim[0])/(self.y_range) * self.resolution[1])
        return (m, n)

    def coordinates_to_number(self, c):
        m = c[0]
        n = c[1]
        x = self.x_lim[0] + m  * 1.0/ self.resolution[0]  * self.x_range
        y = self.y_lim[0] + n * 1.0 / self.resolution[1] * self.y_range
        return x + 1j * y


    def draw_line(self, z1, z2, t):
        c1 = self.number_to_coordinates(z1)
        c2 = self.number_to_coordinates(z2)
        s = abs(c2[0] - c1[0]) + abs(c2[1] - c1[1]) + 1
        for i in range(2 * s):
            m = int((c1[0] * i) / (2 * s) + c2[0] * (1 - i / (2 * s)))
            n = int((c1[1] * i) / (2 * s) + c2[1] * (1 - i / (2 * s)))
            if (m >= 0 and m < self.resolution[0]) and (n >= 0 and n < self.resolution[1]):
                if len(self.data[m,n]) == 0:
                    self.data[m,n].append([t,t])
                elif self.data[m,n][-1][-1] == t - 1:
                    self.data[m,n][-1][-1] = t
                elif self.data[m,n][-1][-1] != t:
                    self.data[m,n].append([t,t]) 

    def pixel_on_map(self, z):
        x = z.real
        y = z.imag
        if (x >= self.x_lim[0] and x < self.x_lim[1]) and (y >= self.y_lim[0] and y < self.y_lim[1]):
            return True
        else:
            return False

    def save_map(self, filename="map"):
        mp = np.zeros(self.resolution)
        for i in range(self.resolution[0]):
            for j in range(self.resolution[1]):
                if len(self.data[i,j]) == 0:
                    mp[i,j] = 1
        fig = plt.figure(figsize=(15, 10), dpi=200)
        plt.imshow(mp.T, cmap='Greys_r')
        fig.savefig("./graphics/" + filename + ".png", dpi=200)