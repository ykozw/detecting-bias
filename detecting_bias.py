# %%
import imageio
import numpy as np
import math
from PIL import Image, ImageFilter

# %%
def gauss_cdf(t):
    tt = abs(t)
    tt /= math.sqrt(2.0)
    x = 1.0 / (1.0 + 0.47047 * tt)
    erf = 1.0 - x * (0.3480242 + x * (-0.0958798 + 0.7478556 * x)) * math.exp(-tt * tt)
    return 1.0 - erf

# %%
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

# %%
def to255(x):
    x = max(min(x,1.0),0.0)
    return int(x * 255)

# %%
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

def main():
    # 1ピクセルをウェルチ用の1サンプルとして、
    # 32x32(=1024ピクセル)のタイルを1024サンプルとして使う
    img0 = imageio.imread('./bin/bpt.hdr')
    img1 = imageio.imread('./bin/pt.hdr')
    (w,h,_) = img0.shape
    tileSize = 32
    numTilePixels = tileSize * tileSize
    wt = w//tileSize
    ht = h//tileSize
    img = Image.new("RGB", (wt, ht)) 
    pix = img.load() #
    #print(img.shape)
    for ty in range(wt):
        for tx in range(wt):
            xb = (tx+0)*tileSize
            xe = (tx+1)*tileSize
            yb = (ty+0)*tileSize
            ye = (ty+1)*tileSize
            # 3チャンネル分
            pvs = []
            for ch in range(3):
                #ch = 2
                t0 = img0[yb:ye,xb:xe,ch].flatten()
                t1 = img1[yb:ye,xb:xe,ch].flatten()
                v0 = np.var(t0,ddof=1)
                v1 = np.var(t1,ddof=1)
                m0 = np.average(t0)
                m1 = np.average(t1)
                # T値を出す
                tmp = v0/numTilePixels + v1/numTilePixels
                if tmp == 0.0 :
                    pvs.append(1.0)
                    continue
                tv = -abs((m0-m1)/math.sqrt(tmp))
                # 自由度を出す
                denm0 = ((numTilePixels - 1.0) * numTilePixels * numTilePixels )
                denm1 = ( (v0 * v0) / denm0 + (v1 * v1) / denm0)
                if denm1 != 0.0:
                    nu = round( (tmp * tmp) / denm1 )
                else:
                    #%tb
                    #sys.exit(-1)
                    nu = 1.0
                # TODO: P値を出す
                replace_with_gauss = False
                gauss_nu_threshold = 0.0
                pv = gauss_cdf(tv) if (replace_with_gauss and nu >= gauss_nu_threshold) else 2.0 * tcdf(tv, nu)
                pvs.append(pv)

            # 3チャンネルで一番低い値
            pv = min(pvs)
            confidence = min(max(1.0 - pv,0.0),1.0)
            # 
            pix[tx,ty] = viridis_quintic(confidence)

    img = img.resize((w, h), Image.NEAREST)
    img.save('pvalue.png')

if __name__ == "__main__":
    main()


# %%
a = np.array([[1,2,3],[4,5,6],[7,8,9]])
print(a)
print(a[:3,:2])


