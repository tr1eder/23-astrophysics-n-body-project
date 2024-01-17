import numpy as np
import matplotlib.pyplot as plt
import ast
import imageio
import os

def plotFromFile(filename):
    with open(filename) as f:
        lines = f.readlines()


    output_file = 'tl/plot_task2.gif'
    size = 2
    writer = imageio.get_writer(output_file, fps=20)
    for i, line in enumerate(lines): 
        arr = ast.literal_eval(line)
        arr_transpose = list(map(list, zip(*arr)))
        x, y, z = arr_transpose

        fig = plt.figure(figsize=(7,7))
        ax = plt.axes(projection='3d')

        ax.scatter3D(x,y,z,s=.1,c='k', label='Particles')
        ax.set_xlim(-size,size)
        ax.set_ylim(-size,size)
        ax.set_zlim(-size,size)
        ax.legend()
        f_name = f"tl/plot_task2_{i}.png"
        plt.savefig(f_name)
        plt.close(fig)

        image = imageio.get_reader(f_name).get_data(0)
        writer.append_data(image)

        print (f"Written file: {f_name:25}")
        

    
    writer.close()
    print(f"Video created {output_file}")


plotFromFile("plot_task2.txt")
