# An unofficial "Detecting Bias in Monte Carlo Renderers using Welch’s t-test"[1] implementation.
# This code fairly mimics the official example code.
# [1] A. Jung, J. Hanika, and C. Dachsbacher, “Detecting Bias in Monte Carlo Renderers using Welch’s t-test,” Journal of Computer Graphics Techniques (JCGT), vol. 9, no. 2, pp. 1–25, Jun. 2020, [Online]. Available: http://jcgt.org/published/0009/02/01/

import math
import argparse
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import imageio

#
def tcdf(t, v):
    b = v / (v + t * t)
    c = 1.0
    s = 1.0
    ioe = v % 2
    k = 2 + ioe
    if v < 1 : 
        return 0.0
    if v >= 4:
        while k <= v - 2:
            c *= b - b / k
            s += c
            k += 2
    c = t / math.sqrt(v)
    if ioe != 1 : 
        return 0.5 + 0.5 * math.sqrt(b) * c * s
    tmp = 0 if v == 1 else b * c * s
    return 0.5 + (tmp + math.atan(c)) / math.pi

#
def to255(x):
    x = max(min(x,1.0),0.0)
    return int(x * 255)

#
def viridis_quintic(x):
    x = min(1.0, x)
    x2 = x * x
    x3 = x2 * x
    x4 = x2 * x2
    x5 = x3 * x2
    return (
        to255(+0.280268003 - 0.143510503 * x + 2.225793877  * x2 - 14.815088879 * x3 + 25.212752309 * x4 - 11.772589584 * x5),
        to255(-0.002117546 + 1.617109353 * x - 1.909305070  * x2 +  2.701152864 * x3 -  1.685288385 * x4 +  0.178738871 * x5),
        to255(+0.300805501 + 2.614650302 * x - 12.019139090 * x2 + 28.933559110 * x3 - 33.491294770 * x4 + 13.762053843 * x5)
    )

#
def diffMain(img0, img1, pvimg, histimg, tileSize=32):
    #
    img0 = imageio.imread(img0)
    img1 = imageio.imread(img1)
    (w,h,_) = img0.shape
    numTilePixels = tileSize * tileSize
    wt = w//tileSize
    ht = h//tileSize
    img = Image.new("RGB", (wt, ht)) 
    pix = img.load()
    pvs = []
    # p-value image
    for ty in range(wt):
        for tx in range(wt):
            xb = (tx+0)*tileSize
            xe = (tx+1)*tileSize
            yb = (ty+0)*tileSize
            ye = (ty+1)*tileSize
            #
            pv3ch = []
            for ch in range(3):
                t0 = img0[yb:ye,xb:xe,ch].flatten()
                t1 = img1[yb:ye,xb:xe,ch].flatten()
                v0 = np.var(t0,ddof=1)
                v1 = np.var(t1,ddof=1)
                m0 = np.average(t0)
                m1 = np.average(t1)
                # t-value
                tmp = v0/numTilePixels + v1/numTilePixels
                if tmp == 0.0 :
                    pv3ch.append(1.0)
                    continue
                tv = -abs((m0-m1)/math.sqrt(tmp))
                # degrees of freedom
                denm0 = ((numTilePixels - 1.0) * numTilePixels * numTilePixels )
                denm1 = ( (v0 * v0) / denm0 + (v1 * v1) / denm0)
                nu = round( (tmp * tmp) / denm1 )
                # p-value
                pv = 2.0 * tcdf(tv, nu)
                pv3ch.append(pv)
            # 
            pv = min(pv3ch)
            pvs.append(pv)
            confidence = min(max(1.0 - pv,0.0),1.0)
            pix[tx,ty] = viridis_quintic(confidence)

    img = img.resize((w, h), Image.NEAREST)
    img.save(pvimg)

    # p-value histgram
    mu = sum(pvs)/float(len(pvs))
    fig  = plt.figure(figsize=(16,16))
    ax = plt.axes()
    ax.set_xlim([0.0,1.0])
    ax.plot(mu, 0.0, "ro")
    ax.hist(pvs,bins=32, range=(0.0, 1.0))
    plt.grid(which='major', ls='solid', color='lightgray')
    plt.savefig(histimg)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('img0', type=str, help='input image 0(*.hdr)')
    parser.add_argument('img1', type=str, help='input image 1(*.hdr)')
    parser.add_argument('--pvimg', type=str, default='pvalue.png',  help='output pvalue image(*.png)')
    parser.add_argument('--histimg', type=str, default='hist.png',  help='output histgram image(*.png)')
    parser.add_argument('--tilesize', type=int, default=32,  help='tile size')
    args = parser.parse_args()
    diffMain(args.img0, args.img1, pvimg=args.pvimg, histimg=args.histimg, tileSize=args.tilesize)

if __name__ == "__main__":
    main()
