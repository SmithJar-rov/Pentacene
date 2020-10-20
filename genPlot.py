import pentacene as pen
import numpy as np
import math
import os

def main(aggWeight, aggDistance, aggCount, aggChance):
    #after some testing, it seems like the desired effect does
    #not happen until aggDistance > 0.20
    #so that is why this is our lower bound
    for b in np.linspace(0.20, aggDistance, 6):
        folder = "distance=%s/"%(b)
        path = os.path.join('./plots/',folder)
        os.mkdir(path)
        for a in np.linspace(1/10*aggWeight, aggWeight, 5):
            for c in np.linspace(1,aggCount,5):
                for d in np.linspace(1/10*aggChance,aggChance,5):
                    print(a,b,c,d)
                    plt=pen.main(a,b,math.floor(c),d)
                    plt.savefig(os.path.join(path,"aggWeight=%s,count=%s,chance=%s.png"%(a,c,d)))
                    plt.close()

if __name__ == '__main__':
    main(1,0.70,25,1)
