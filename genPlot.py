import pentacene as pen
import numpy as np
import math
import os

def main(penalty, aggDistance, aggCount, minChance):
    for b in np.linspace(0.01, aggDistance, 10):
        folder = "distance=%s/"%(b)
        path = os.path.join('./plots/',folder)
        os.mkdir(path)
        for a in np.linspace(1/2*penalty, penalty, 3):
            for c in np.linspace(1,25,5):
                for d in np.linspace(1/2*minChance,minChance,3):
                    print(a,b,c,d)
                    plt=pen.main(a,b,math.floor(c),d)
                    plt.savefig(os.path.join(path,"penalty=%s,count=%s,chance=%s.png"%(a,c,d)))
                    plt.close()

if __name__ == '__main__':
    main(0.5,0.25,25,0.5)
