import rebound
import numpy as np
sim = rebound.Simulation()
OMEGA = 0.00012843527     # [1/s]
sim.ri_sei.OMEGA = OMEGA
surface_density = 523.    # kg/m^2
particle_density = 400.   # kg/m^3
sim.G = 6.67428e-11       # N m^2 / kg^2
sim.dt = 1e-3*2.*np.pi/OMEGA
sim.softening = 0.1       # [m]
boxsize = 500.            # [m]
sim.configure_box(boxsize)
#sim.configure_ghostboxes(2,2,0)
sim.integrator = "sei"
sim.boundary   = "shear"
sim.gravity    = "tree"
sim.collision  = "tree"
sim.collision_resolve = "hardsphere"

def cor_bridges(r, v):
        eps = 0.32*pow(abs(v)*100.,-0.234)
        if eps>1.:
            eps=1.
        if eps<0.:
            eps=0.
        if abs(v)<0.01:
            eps=0.0001
                
        return eps
sim.coefficient_of_restitution = cor_bridges


def powerlaw(slope, min_v, max_v):
    y = np.random.uniform()
    pow_max = pow(max_v, slope+1.)
    pow_min = pow(min_v, slope+1.)
    return pow((pow_max-pow_min)*y + pow_min, 1./(slope+1.))

total_mass = 0.
while total_mass < surface_density*(boxsize**2):
    radius = powerlaw(slope=-2.85, min_v=1, max_v=5)  # [m]    
    mass = particle_density*4./3.*np.pi*(radius**3)
    x = np.random.uniform(low=-boxsize/2., high=boxsize/2.)
    sim.add(
        m=mass,
        r=radius,
        x=x,
        y=np.random.uniform(low=-boxsize/2., high=boxsize/2.),
        z=np.random.normal(),
        vx = 0.,
        vy = -3./2.*x*OMEGA, 
        vz = 0.)
    total_mass += mass

    import matplotlib.pyplot as plt
import matplotlib.patches as patches
def plotParticles(sim):
    font1 = {'size':20}
    fig = plt.figure(figsize=(8,8))
    ax = plt.subplot(111,aspect='equal')
    ax.set_xlabel("Radial coordinate [m]",fontdict=font1)
    ax.set_ylabel("Azimuthal coordinate [m]",fontdict=font1)
    ax.set_ylim(-boxsize/2.,boxsize/2.)
    ax.set_xlim(-boxsize/2.,boxsize/2.)

    for i, p in enumerate(sim.particles):
        circ = patches.Circle((p.x, p.y), p.r, facecolor='darkgray', edgecolor='black')
        ax.add_patch(circ)
#plt.ioff()
plotParticles(sim)
plt.close()
step=3
font2 = {'size':15}
font1 = {'size':20}


for i in range(int(6.*np.pi/(OMEGA*sim.dt*step))):
        sim.steps(step)   #2.*np.pi/OMEGA
        plotParticles(sim)
        t=i*step*sim.dt/(2*np.pi/(OMEGA))
        plt.title('t = %.3f [orb]' % t,fontdict=font2,loc='right')
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        if i == 152:
                plt.savefig('wakes{}'.format(i))
        

        plt.close()
print('done')
