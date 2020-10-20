from vpython import *
import numpy as np
import random
import time
from math import *
import matplotlib.pyplot as plt

class space(object):
    def __init__(self,weight,aggWeight,aggDistance,aggCount,visual,gridSize):
        self.grid = []
        self.totalagg = []
        self.aggDistance = aggDistance
        self.aggCount = aggCount
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
#        return self.weight*self.aggWeight*len(self.totalagg) + self.weight*(len(self.grid) - len(self.totalagg))

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

def funWrap(argMax, argMin, aggCount):
    def linear(agg):
        if len(agg) < aggCount:
            return (argMin-argMax)/aggCount * len(agg) + argMax
        else:
            return argMin
    return linear

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

def main(aggWeight, aggDistance, aggCount, aggChance,
         visual = False, returnPlot = True):
    #population size in grid is 200 molecules
    popSize = 200
    #the grid size is 0.1
    gridSize = 0.1
    #light emission of a non-aggregate molecule
    weight = 10
    #time in milliseconds
    t=0
    #timestep in milliseconds
    dt = 100
    #chance of decay per millisecond per molecule
    k = 1e-2

    chance = funWrap(k*dt, k*dt*aggChance, aggCount)
    pl = funWrap(weight, weight*aggWeight, aggCount)

    #create an populate the world
    world = space(weight, aggWeight, aggDistance*gridSize, aggCount, visual, gridSize)
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

    if returnPlot:
        return(fPlot(tpoints, plpoints))

if __name__ == '__main__':
    plt=main(0.25,0.01,1,0.5,True)
    plt.show()
