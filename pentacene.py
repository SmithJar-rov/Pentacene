from vpython import *
import numpy as np
import random
import time
from math import *
import matplotlib.pyplot as plt

class space(object):
    def __init__(self,weight,penalty,aggDistance,aggCount,visual,gridSize):
        self.grid = []
        self.totalagg = []
        self.aggDistance = aggDistance
        self.aggCount = aggCount
        self.weight = weight
        self.penalty = penalty
        self.gridSize = gridSize
        self.pl = 0
        if visual:
            self.someObject = vMol
            self.scene=canvas(title = 'Pentacene',
                     height = 600,
                     width = 600,
                     range = gridSize * 2)
        else:
            self.someObject = mol

#    old pl function
#    def pl(self):
#        return self.weight*self.penalty*len(self.totalagg) + self.weight*(len(self.grid) - len(self.totalagg))

class mol(object):
    def __init__(self, gridSize):
        self.agg = []
        self.cord(vector(
            gridSize*random.random(),
            gridSize*random.random(),
            gridSize*random.random()))

    def cord(self, someVector):
        self.pos = someVector

    def color(self, someVector):
        self.color = someVector

class vMol(mol):
    def cord(self, someVector):
        self.pos = someVector
        self.bol = sphere(pos=self.pos,
                            radius=0.0025,
                            visible=True,
                            color=vec(1,1,1))
    def color(self, someVector):
        self.color = someVector
        self.bol.color = someVector

def funWrap(arg):
    #arg is equal to the maximum value of the function,
    #whether it be maximum pl or maximum decay chance
    #both decrease with respect to the amount of nearby molecules
    #so the function can be used for both
    def sigmoid(agg):
        return 1/2*arg + (arg)/(1.0+exp(len(agg)))
    return sigmoid

def populate(world, popSize, pl):
    #populate list containing molecules
    for i in range (0,popSize):
       world.grid = world.grid + [world.someObject(world.gridSize)]
       for j in world.grid:
           #check euclidian distance between molecule k and molecule j
           radius = mag(world.grid[i].pos - j.pos)
           #check proximity via border wrapping
           if world.grid[i].pos.x > j.pos.x:
               rx = mag(world.grid[i].pos - (j.pos+vector(world.gridSize,0,0)))
           else:
               rx = mag(world.grid[i].pos - (j.pos-vector(world.gridSize,0,0)))
           if world.grid[i].pos.y > j.pos.y:
               ry = mag(world.grid[i].pos-(j.pos+vector(0,world.gridSize,0)))
           else:
               ry = mag(world.grid[i].pos-(j.pos-vector(0,world.gridSize,0)))
           if world.grid[i].pos.z > j.pos.z:
               rz = mag(world.grid[i].pos - (j.pos+vector(0,0,world.gridSize)))
           else:
               rz = mag(world.grid[i].pos - (j.pos-vector(0,0,world.gridSize)))
           if (radius<=world.aggDistance or rx<=world.aggDistance or ry<=world.aggDistance or rz<=world.aggDistance) and radius!=0:
               world.grid[i].agg=world.grid[i].agg + [j]
               j.agg = j.agg + [world.grid[i]]

    #populate list containing aggregates and calculate initial pl value
    for i in world.grid:
        if len(i.agg)>=world.aggCount:
            world.totalagg = world.totalagg + [i]
        world.pl = world.pl + pl(i.agg)
    return world

def decay(world, chance, pl):
    world.pl = 0
    for i in world.grid:
        roll=100*random.random()
        if roll <= chance(i.agg):
            if len(i.agg)>=world.aggCount:
                world.totalagg.remove(i)
            for k in i.agg:
                before = len(k.agg)
                k.agg.remove(i)
                after = len(k.agg)
                if before >= world.aggCount and after < world.aggCount:
                    world.totalagg.remove(k)
            i.color(color.red)
            world.grid.remove(i)
        else:
            world.pl = world.pl + pl(i.agg)
    return world

def fPlot(tpoints, plpoints):
   plt.plot(tpoints,plpoints)
   plt.yscale('log')
   plt.xlabel('time (ms)')
   plt.ylabel('Pl')
   return plt

def main(weight, penalty, aggDistance, aggCount,
         visual = False, plot = False, save = False, popSize = 200, gridSize = 0.1):
    #time in milliseconds
    t=0
    #timestep in milliseconds
    dt = 100
    #chance of decay per millisecond per molecule
    k = 1e-2

    chance = funWrap(k*dt)
    pl = funWrap(weight)

    #create an populate the world
    world = space(weight, penalty, aggDistance, aggCount, visual, gridSize)
    world = populate(world, popSize, pl)

    #array containing time
    tpoints=np.array(t)
    #containing pl
    plpoints=np.array(world.pl)

    while (len(world.grid)!=0):
        world = decay(world, chance, pl)
        #print("Grid Total = ",len(world.grid))
        t = t + dt
        tpoints = np.append(tpoints,t)
        plpoints = np.append(plpoints, world.pl)

    if plot:
        plt=fPlot(tpoints, plpoints)
        plt.show()

    if save:
        plt=fPlot(tpoints,plpoints)
        plt.savefig(path)

if __name__ == '__main__':
    path = "/home/jared/projets/pentacene/plot.png"
    main(10,1/7,0.01,2,True,True,False)
