import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.spatial.distance import cdist

# Parameters
particles = 40
positions = np.random.uniform(-50, 50, size=(particles, 2))
angles = np.random.uniform(0, 2 * np.pi, size=particles)
velocities = np.ones(particles) * 20
dt = 0.005
# Gradient function (attractive force towards the origin(food))
def gradient(x, y):
    r = np.clip(np.sqrt(x**2 + y**2), 0.1, None) #prevented division by zero
    return 10/r

fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(-50, 50)
ax.set_ylim(-50, 50)
ax.scatter(0, 0, color='red', s=100, label='food')
l = 1.5
w = 0.3
quiver =ax.quiver(positions[:, 0], positions[:, 1], l * np.cos(angles), l * np.sin(angles),
                    color='blue', scale=80,
                    headwidth = 3, headlength = 4,
                    headaxislength = 3, width = 0.0025)

# Update function for animation
def update(frame):
    global positions, angles, dt, velocities, particles, gradients
    #Brownian motion
    angles += np.random.normal(0, 0.2, particles)  
    
    gradients = gradient(positions[:, 0], positions[:, 1])
    #movement after saving current gradient values
    dx = np.cos(angles) * velocities * dt
    dy = np.sin(angles) * velocities * dt
    positions += np.column_stack((dx, dy))

    #collision detection and response
     #checking if the head of a particle is in the proximity of another particle
    head_points = np.column_stack((positions[:, 0] + (l/2 * np.cos(angles)), positions[:, 1] + (l/2 * np.sin(angles))))
    dist_matrix = cdist(positions, head_points)
    indices1 = np.argwhere((dist_matrix <= (l/2)) & (np.eye(particles) == 0))
    for centre, head in indices1:
        centres = positions[centre]
        heads = head_points[head]
        angles_centre = angles[centre]
        angles_head = angles[head]
        nx = heads[0] - centres[0]
        ny = heads[1] - centres[1]
        localx = nx * np.cos(-angles_centre) - ny * np.sin(-angles_centre)
        localy = nx * np.sin(-angles_centre) + ny * np.cos(-angles_centre)
        if -w/2 < localy < w/2 and -l/2 < localx < l/2:
            distance_collision = localx/(l/2)
        else:
            continue
        #if the head collides at centre of particle, it gets pushed translationally
        if distance_collision == 0:
            dp = 2 * dt * ((velocities[head]-10)/30) * np.cos(angles_head - angles_centre)
            positions[centre, 0] += dp * np.cos(angles_head)
            positions[centre, 1] += dp * np.sin(angles_head)
        #if the head collides other than centre of particle, it changes angle due to push
        else:
            angles[centre] += 2 * dt * (-1 * distance_collision) * ((velocities[head]-10)/30) * np.sin(angles_centre - angles_head)
       
       #the particle which is hitting gets slowed down
        velocities[head] -= velocities[head] * np.sin(angles_centre - angles_head) * 0.5
       #its angle changes slightly towards direction of particle that got hit 
        angles[head] += (((angles_centre - angles_head) + np.pi) % (2 * np.pi)) - np.pi
   
   #decision making based on gradient sensing
    changes = gradient(positions[:, 0], positions[:, 1]) - gradients
    inc_mask = changes > 0
    dec_mask = changes < 0
    angles[dec_mask] += np.random.uniform(-np.pi, np.pi, np.sum(dec_mask))
    velocities[inc_mask] += 0.7
    velocities[dec_mask] -= 0.7
    #clippping velocities according to the speed of a 1.5 micron bacteria
    velocities = np.clip(velocities, 10, 40)
    quiver.set_UVC(l * np.cos(angles), l * np.sin(angles))
    quiver.set_offsets(positions)
    return quiver,
ani = FuncAnimation(fig, update, frames=500, interval=60, blit=False)
plt.show()
