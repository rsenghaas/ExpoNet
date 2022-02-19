# Executable file
# This file should only contain function calls 
# imported from other files

# Include modules 
import src.sw_curve as sw
import src.path as pth
import src.network as nw
import src.pathmap as mp

# Include default libraries (ideally not much needed)
import numpy as np

# Main function
def main():
    dt = 1e-8
    steps = 50000
    theta = 0.1
    expo_net = nw.network(sw.H_c3, dt, steps, theta, expo=True)
    expo_net.start_paths()
    expo_net.evolve()
    print("Evolution_finished!")
    expo_net.plot_network(paths=[], steps=0, fix_axis=True, filename="expo_net_lq")
    expo_net.plot_all_networks()
    for i in range(len(expo_net.all_intersections)):
        print(expo_net.all_intersections[i].path_index)
        print(expo_net.all_intersections[i].target_map)
        print("Intersection {}".format(i))
        for j in range(len(expo_net.all_intersections[i].coordinates)):
            print(expo_net.all_intersections[i].steps, expo_net.all_intersections[i].coordinates[j])


if __name__ == '__main__':
    main()