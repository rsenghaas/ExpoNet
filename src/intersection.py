# Class to keep track of (possible intersections)

class intersection():
    def __init__(self, t, path_index, target_map, coordinates):
        self.steps = [t, t]
        self.path_index = path_index
        self.target_map = target_map
        self.coordinates = [coordinates]
        self.next = -1
        

