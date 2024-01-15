import numpy as np
import matplotlib.pyplot as plt
import ast
import imageio
import os

def plotFromFile(filename):
    with open(filename) as f:
        lines = f.readlines()


    output_file = 'tl/plot_task2.gif'
    writer = imageio.get_writer(output_file, fps=5)
    for i, line in enumerate(lines): 
        arr = ast.literal_eval(line)
        arr_transpose = list(map(list, zip(*arr)))
        x, y, z = arr_transpose

        fig = plt.figure(figsize=(7,7))
        ax = plt.axes(projection='3d')

        ax.scatter3D(x,y,z, label='Particles')
        ax.set_xlim(-2,2)
        ax.set_ylim(-2,2)
        ax.set_zlim(-2,2)
        ax.legend()
        f_name = f"tl/plot_task2_{i}.png"
        plt.savefig(f_name)
        plt.close(fig)

        image = imageio.get_reader(f_name).get_data(0)
        writer.append_data(image)

        print ("Written file ", f_name)
        

    
    writer.close()
    print(f"Video created {output_file}")


plotFromFile("plot_task2.txt")
