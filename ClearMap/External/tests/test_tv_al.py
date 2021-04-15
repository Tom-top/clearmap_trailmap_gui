import numpy as np


from ClearMap.External.pydeconv import deconv_tv_al, deconv_wiener, deconv_rl

from ClearMap.External.utils import myconvolve, psf


if __name__ == '__main__':
    #2d
    from matplotlib.pyplot import imread

    s = .1
    x = imread("./ClearMap/External/tests/data/usaf.png")

    x*= 255
    
    h = psf(x.shape,(3.,11))


    h = psf(x.shape,(5.,5.))

    h2 = psf(x.shape,(11,3))
    
    y = myconvolve(x,h)+s*np.amax(x)*np.random.uniform(0,1,x.shape)

    y2 = myconvolve(x,h2)+s*np.amax(x)*np.random.uniform(0,1,x.shape)

    # u = deconv_tv_al([y,y2],[h,h2])

    
    #u = deconv_tv_al(y,h, mu = 10., rho = 1., n_threads = 24)

    #u2 = deconv_wiener(y,h,0.1);
    u = deconv_rl(y,h, gamma = 10e-7, mult_mode='root')
    
    import matplotlib.pyplot as plt
    plt.figure(1); plt.clf();
    plt.subplot(1,3,1);
    plt.imshow(x)
    plt.subplot(1,3,2);
    plt.imshow(y);
    plt.subplot(1,3,3);
    plt.imshow(u)
    
    
    
    # test on 3d data
    
    img = data_raw[:96,:96,:96];
        
    imgd = deconv_tv_al(img,h,mu = 20., rho = 1., n_threads = 24)

    dv.DataViewer([img])
    dv.DataViewer([imgd])
    
    #3d
    from spimagine import read3dTiff

    x = read3dTiff("data/usaf3d.tif")[100:228,100:228,100:228]
    h = psf(x.shape,(3.,3.,11))

    h2 = psf(x.shape,(11,3,3))
    h3 = psf(x.shape,(3,11,3))
    
    y = myconvolve(x,h)+.1*np.amax(x)*np.random.uniform(0,1,x.shape)
    y2 = myconvolve(x,h2)+.1*np.amax(x)*np.random.uniform(0,1,x.shape)
    y3 = myconvolve(x,h3)+.1*np.amax(x)*np.random.uniform(0,1,x.shape)

    u = deconv_tv_al([y,y2,y3],[h,h2,h3],100.,1.)

    u2 = deconv_wiener([y,y2,y3],[h,h2,h3],0.01)




    
