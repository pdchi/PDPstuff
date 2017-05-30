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
fig = plt.figure()
ax = plt.axes(xlim=(-3+659, 3+659), ylim=(-.5, 10))
plt.xlabel('Frequency (Hz)')
plt.ylabel('FFT Amplitude')
line, = ax.plot([], [], lw=2)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
    x = np.linspace(-3+659, 3+659, 900)
    #mu = .65*np.sin(2.*np.pi*(1-.03*i))+650
    mu = 659
    ets = .65*np.arcsin(1-.01*i)
    y = 3.*np.exp(-np.power(x - mu, 2.) / (2 * np.power(.05, 2.)))
    y += 3.*np.exp(-np.power(x - mu+ets, 2.) / (2 * np.power(.05, 2.)))
    y += 3.*np.exp(-np.power(x - mu-ets, 2.) / (2 * np.power(.05, 2.)))
    #y = np.sin(2 * np.pi * (x - 0.01 * i))
    y += np.random.normal(scale=.1,size=len(x))
    line.set_data(x, y)
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=200, interval=20, blit=False)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
#anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()
