"""
Matplotlib Animation Example

author: Jake Vanderplas
email: vanderplas@astro.washington.edu
website: http://jakevdp.github.com
license: BSD
Please feel free to use and modify this, but keep the above information. Thanks!
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import seaborn
seaborn.set_style('white')
seaborn.set_style('ticks')
# First set up the figure, the axis, and the plot element we want to animate
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(13,7))
#ax1 = plt.axes(xlim=(-3+659, 3+659), ylim=(-.5, 10))
ax1.set_xlim(-3+659, 3+659)
ax1.set_ylim(-.5, 10)
ax2.set_xlim(-2, 2)
ax2.set_ylim(-1, 1)
ax2.set_axis_off()
ax1.set_xlabel('Frequency (Hz)')
ax1.set_ylabel('FFT Amplitude')
line1, = ax1.plot([], [], lw=2)
line2, = ax2.plot([], [], lw=2)
line = [line1,line2]

# initialization function: plot the background of each frame
def init():
    line[0].set_data([], [])
    line[1].set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
    x = np.linspace(-3+659, 3+659, 900)
    mu = .65*np.sin(2.*np.pi*(1-.03*i))+659
    y = 9.*np.exp(-np.power(x - mu, 2.) / (2 * np.power(.05, 2.)))
    #y = np.sin(2 * np.pi * (x - 0.01 * i))
    y += np.random.normal(scale=.1,size=len(x))

    xp = 1.5*np.sin(2.*np.pi*(1-.03*i)+np.pi)
    sc = .1-(.04*np.cos(2.*np.pi*(1-.03*i)+np.pi))
    x2 = [xp,xp+sc,xp,xp-sc,xp,xp,xp-sc,xp+sc]
    #y2 = 9.*np.exp(-np.power(x - mu, 2.) / (2 * np.power(.05, 2.)))
    #y2 += np.random.normal(scale=.1,size=len(x))
    yp = .07*np.cos(2.*np.pi*(1-.03*i)+np.pi)
    y2 = [yp+sc*.5,yp,yp-sc*.5,yp,yp+sc*.5,yp-sc*.5,yp,yp]
    line[0].set_data(x, y)
    line[1].set_data(x2, y2)
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=200, interval=10, blit=False)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
#anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()
