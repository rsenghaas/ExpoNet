# Include modules
import src.pathmap as mp
import src.diffeq as dq
import src.sw_curve as sw

# Include default libraries
import numpy as np
# FAQ: Should we use a dynamic list here for x and y's?
# Pro: More flexible if we use a different cutoff condition than steps
# Cons: Performance could suffer.
class path():
    def __init__(self, length, label, gen, x_limits, y_limits, resolution, color='black'):
        self.x = np.zeros(length, dtype=complex)
        self.y_i = np.zeros(length, dtype=complex)
        self.y_j = np.zeros(length, dtype=complex)
        self.label = label
        self.gen = gen
        self.pathmap = mp.pathmap(x_limits, y_limits, resolution)
        self.color = color
        self.intersection_tracker = []
        self.active = True

    def evolve(self, step, dt, f, theta, expo):
        val_x = self.x[step - 1]
        val_y_i = self.y_i[step - 1]
        val_y_j = self.y_j[step - 1]
        dy = dq.rk4(f, [val_x, val_y_i, val_y_j], dt, theta, expo=expo)
        if dy[0] == 0:
             self.active = False
             return
        self.x[step] = val_x + dy[0]
        self.y_i[step] = val_y_i + dy[1]
        self.y_j[step] = val_y_j + dy[2]
        if abs(self.x[step] > 1000) or abs(self.x[step]) < 1/1000 or dy[1] > 10 or dy[2] > 10:
            self.active = False

    def update_map(self, step):
        self.pathmap.draw_line(self.x[step - 1], self.x[step], step)
    

    