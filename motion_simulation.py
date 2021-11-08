import numpy as np
import numpy.random as ran
import scipy.constants as con

###############################################################################

class particle():
    '''
    this class holds static data about each particle such as mass and charge,
    as well as position, velocity, acceleration, and net experienced electric field data.
    it also keeps an array with the history of each of the changing quantities
    '''

    def __init__(self,params):
        '''
        pass in array with all initial conditions
        [mass,charge,r0,v0,E,a0]
        where r0 = [x0,y0,z0],
        v0 = [vx0,vy0,vz0],
        E = [Ex,Ey,Ez], and
        a0 = [ax0,ay0,az0]
        use standard SI units for all values
        '''
        self.m = params[0]
        self.q = params[1]
        self.r = np.array(params[2])
        self.v = np.array(params[3])
        self.E_net = np.array(params[4])
        self.r_arr = [self.r]
        self.v_arr = [self.v]
        self.a_arr = [params[5]]
        self.E_arr = [self.E_net]
    
    def calc_a(self):
        '''
        this function calculates the instantaneous acceleration experienced by a particle
        based on its charge, mass, net experienced electric field, and gravity.
        note that y-axis is treated as verticle due to software that animates these motions.
        '''
        self.a = np.array([
            self.q*self.E_net[0]/self.m,
            (self.q*self.E_net[1]/self.m)-g,
            self.q*self.E_net[2]/self.m
        ])
    
    def calc_v(self,delta_t):
        '''
        this function calculates the velocity of a particle at the next time step
        based on its current velocity, acceleration, and time step size
        '''
        self.v = self.v + delta_t*self.a
    
    def calc_r(self,delta_t):
        '''
        this function calculates the position of a particle at the next time step
        based on its current position, velocity, acceleration, and time step size
        '''
        self.r = self.r + delta_t*self.v + 0.5*(delta_t**2)*self.a
        if self.r[1] < 0:#check for negative y value
            if abs(self.r[0])<0.25:#check if x value is on panel surface
                if abs(self.r[2])<0.25:#check if z value is on panel surface
                    self.floor_collision()
                else:
                    pass
            else:
                pass
        else:
            pass
    
    def floor_collision(self):
        '''
        this function is called in the update to the position vector.
        if the vertical position (y-direction) is found to be negative
        at the given time step, the particle's vertical
        position and velocity are set to zero
        '''
        self.r[1] = 0
        self.v[1] = 0

    def update(self,delta_t):
        '''
        this function runs all of the functions which update its motion parameters
        based on its current net experienced electric field, previous motion parameters,
        and time step size.
        pass in the desired time step size.
        at each time step, a loop in the motion simulation function will update each particle's
        net experienced electric field prior to running this function.
        finally, this function will append the updated motion parameters
        to arrays within the particle to reference later.
        '''
        self.calc_a()
        self.calc_r(delta_t)
        self.calc_v(delta_t)
        self.r_arr.append(self.r)
        self.v_arr.append(self.v)
        self.a_arr.append(self.a)



###############################################################
#the following functions are all used for physics calculations
def separation_vector(r_electrode,r):
    '''
    this finds the separation vector between a point in the yz
    plane and an electrode on the panel surface.
    because the electrodes are treated as infinite in the x direction,
    this is treated as a two dimensional vector problem
    and all x values are set to zero.
    please input numpy vectors for both.
    '''
    r_yz = np.array([0,r[1],r[2]])
    r_sep = r_yz - r_electrode
    r_sep_magnitude = np.linalg.norm(r_sep)
    r_sep_hat = r_sep/r_sep_magnitude
    return r_sep_magnitude,r_sep_hat

def field_from_electrode(V,r_electrode,r):
    '''
    takes in the electrode voltage,
    the position of the electrode (np array vector of form [x,0,0]),
    and position where field is being calculated.
    currently based on potential from an infinite line
    '''
    r_sep_magnitude,r_sep_hat = separation_vector(r_electrode,r)
    E_electrode = -(V/r_sep_magnitude)*r_sep_hat
    return E_electrode

def V(amplitude,freq,t,phi):
    '''
    function to calculate a sinusoidal voltage at time t
    takes in amplitude, frequency, time, and phase shift
    '''
    argument = freq*t-phi
    return amplitude*np.sin(argument)

def E_circuit_function(r,t,amplitude,frequency):
    '''
    function to calculate the electric field from the entire circuit
    at a location r and time t.
    takes in amplitude of circuit voltage and circuit frequency,
    as well as a np vector where the field is to be found, and the time
    '''
    d = 6*(10**(-3))#6 mm spacing between electrode (center to center)
    V1_locations = [(d/2)+(i-42)*d for i in range(82)]
    V1 = V(amplitude,frequency,t,0)

    V2_locations = [(d/2)+(i-41)*d for i in range(82)]
    V2 = V(amplitude,frequency,t,np.pi/3)#phase shift of 60 degrees

    V3_locations = [(d/2)+(i-40)*d for i in range(82)]
    V3 = V(amplitude,frequency,t,2*np.pi/3)#phase shift of 120 degrees
    
    E_circ = np.array([0.0,0.0,0.0])#initialize field of zero

    for z in V1_locations:
        E_circ += field_from_electrode(V1,np.array([0,-0.001,z]),r)
    for z in V2_locations:
        E_circ += field_from_electrode(V2,np.array([0,-0.001,z]),r)
    for z in V3_locations:
        E_circ += field_from_electrode(V3,np.array([0,-0.001,z]),r)
    return E_circ

def net_field(E_circuit,particle_array,excluded_index):
    '''
    calculates the net electric field from circuit and other particles
    at the position of a particluar particle.
    pass in E_circuit array with x,y,z components of field from circuit,
    array of all particles in simulation,
    and index of particle for which the net experienced field is being calculated
    '''
    net_field = E_circuit
    
    r_excluded = np.array([
        particle_array[excluded_index].r[0],
        particle_array[excluded_index].r[1],
        particle_array[excluded_index].r[2]
    ])
    
    for i,p in enumerate(particle_array):
        if i == excluded_index:
            pass
        else:
            r = r_excluded-p.r
            r_magnitude = np.linalg.norm(r)
            r_hat = r/r_magnitude
            net_field += k*p.q*r_hat/(r_magnitude**2)
            
    return net_field

def motion_simulation(particle_array,total_time,delta_t,amplitude,frequency,filename):
    '''
    this function runs the motion simulation for all particles in the particle array.
    pass in array of all particle objects,
    total time to simulate motion, time step size (delta_t),
    and the magnitude of the electric field from the circuit
    (the specifics of the circuit's electric field are subject to change)
    '''
    n_steps = int(total_time/delta_t)

    with open(filename,"a") as f:
        for i,p in enumerate(particle_array):
            f.write(str(p.r[0])+','+str(p.r[1])+','+str(p.r[2])+'\n')

    for i in range(n_steps):
        #update electric field at each particle location
        for j,p in enumerate(particle_array):
            p.E_net = net_field(E_circuit_function(p.r,delta_t*i,amplitude,frequency),particle_array,j)
        #update each particle's motion quantities
        for k,p in enumerate(particle_array):
            p.update(delta_t)
            with open(filename,"a") as f:
                f.write(str(p.r[0])+','+str(p.r[1])+','+str(p.r[2])+'\n')
        #save all particle's position data to file



#####################################################################################
#the following functions are used to help initialize particles with random parameters
def random_about_mean(mean,percent):
    '''
    takes in a mean value, and returns a float within
    +/- a percent of that mean.
    pass in percent as a decimal (ex: 5 percent as 0.05)
    assumes uniform distribution on the interval
    '''
    delta = percent*mean
    return ran.uniform(mean-delta,mean+delta)

def random_within_bounds(boundary):
    '''
    takes in a boundary and returns a float between +/- that boundary from 0
    assumes uniform distribution
    '''
    return ran.uniform(-boundary,boundary)

def random_r0(boundary):
    '''
    returns a 3D vector with 0 y component
    and randomized values for x and z (within +/- boundary)
    '''
    r = [
        random_within_bounds(boundary),
        0,
        random_within_bounds(boundary)
    ]
    return r


def random_particle(m_data,q_data,boundary,E_data):
    '''
    returns partice with randomized initial conditions
    pass in m_data=[m_avg,m_range]
    q_data=[q_avg,q_range]
    boundary for max distance panel extends from origin
    E_data=[amplitude,frequency]
    '''
    mass = random_about_mean(m_data[0],m_data[1])
    charge = random_about_mean(q_data[0],q_data[1])
    r0 = random_r0(boundary)
    v0 = [0,0,0]
    E = E_circuit_function(r0,0,E_data[0],E_data[1])
    a0 = [0,0,0]
    return particle([mass,charge,r0,v0,E,a0])



###############################################################################
###############################################################################

k = 1/(4*con.pi*con.epsilon_0)
g = 0.166*con.g#gravitational acceleration on lunar surface

mass_density = 2650#kg/m^3
d = 20#micro m
r = (d/2)*(10**(-6))#m
v = (4*con.pi/3)*((r)**3)#m^3
m_avg = mass_density*v#kg

charge_density = 20*(10**(-6))#C/kg
q_avg = charge_density*m_avg#C

def main():

    mass_data = [m_avg,0.1]
    charge_data = [q_avg,0.05]
    boundary = 0.25
    amplitude = 7.5*(10**(2))#V
    frequency = 500#Hz
    E_data = [amplitude,frequency]

    n_particles = 300

    particle_array = [
        random_particle(mass_data,charge_data,boundary,E_data) for i in range(n_particles)
    ]

    total_time = 1200#s
    delta_t = 0.01#s

    filename = 'coordinates'+str(n_particles)+'particles'+str(amplitude)+'V'+str(frequency)+'Hz'+str(total_time)+'s.csv'

    motion_simulation(particle_array,total_time,delta_t,amplitude,frequency,filename)

main()