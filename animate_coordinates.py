from vpython import *
import time


#update these values based on the simulation you're animating
n_particles = 300
total_time = 1200#s
delta_t = 0.01#s
filename='coordinates300particles750.0V500Hz1200s.csv'


def create_coordinate_tensor(lines,n_steps,n_particles):
    #initialize an empty array to store the coordinates
    coords = []

    #fill the coordinate array as a tensor with a matrix of particle coordinates at each time step
    for t in range(n_steps):
        t_coords = []
        for n in range(n_particles):
            t_coords.append([float(v) for v in lines[t+n].split(',')])
        coords.append(t_coords)
    
    return coords

def spawn_background():
    #create background objects for the animation
    axis_length = 0.125
    x_axis = arrow(pos=vector(0,0,0),axis=vector(1,0,0),length=axis_length,color=color.green)
    y_axis = arrow(pos=vector(0,0,0),axis=vector(0,1,0),length=axis_length,color=color.yellow)
    z_axis = arrow(pos=vector(0,0,0),axis=vector(0,0,1),length=axis_length,color=color.red)
    floor = box(pos=vector(0,-0.01,0),length=0.5,height=0.01,width=0.5,color=color.blue)

def set_camera():
    #set the camera position and orientation
    scene.camera.pos = vector(0.4,0.4,0.4)
    scene.camera.axis = vector(-1,-1.3,-1)

def animation_loop(coords,ani_array):
    #iterating over each timestep
    for i,timestep in enumerate(coords):
        #iterating over each animated particle
        for j,ani in enumerate(ani_array):
            ani.pos=vector(
                timestep[j][0],
                timestep[j][1],
                timestep[j][2],
            )


def main(n_particles,total_time,delta_t,filename):
    #calculating values for the animation and tensor creation
    n_steps = int(total_time/delta_t)
    interval = (total_time/n_steps)

    #read in the file with the coordinates
    with open(filename) as f:
        lines = f.readlines()
        f.close()
    
    coords = create_coordinate_tensor(lines,n_steps,n_particles)
    
    spawn_background()

    set_camera()

    #initialize an array with all the particles to be animated based on t=0
    ani_array = [
        sphere(pos=vector(r[0],r[1],r[2]),radius=0.01,color=vector(0.5,0.5,0.5)) for r in coords[0]
    ]

    animation_loop(coords,ani_array)


main(n_particles,total_time,delta_t,filename)