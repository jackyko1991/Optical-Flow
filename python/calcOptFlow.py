# !/usr/bin/env python
"""
Michael Hirsch
Example calculations of optical flow, starting with Horn Schunk Optical Flow using OpenCV
"""
try:
    import cv2
    from cv2 import cv
except Exception:
    pass
from numpy import asarray,dstack
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import os

def optflowHornSchunk(new,ref,uv,smoothing=0.01):
    """
    http://docs.opencv.org/modules/legacy/doc/motion_analysis.html
    ***************************
    Note that smoothness parameter for cv.CalcOpticalFlowHS needs to be SMALLER than matlab
    to get similar result. Useless when smoothness was 1 in python, but it's 1 in Matlab!
    *****************************
    """
    cvref = cv.fromarray(ref)
    cvnew = cv.fromarray(new)
    #result is placed in u,v
    # matlab vision.OpticalFlow Horn-Shunck has default maxiter=10, terminate=eps, smoothness=1
    cv.CalcOpticalFlowHS(cvref, cvnew, False, uv[0], uv[1],
                         smoothing,
                         (cv.CV_TERMCRIT_ITER | cv.CV_TERMCRIT_EPS, 8, 0.1))

    # reshape to numpy float32, xpix x ypix x 2
    return dstack((asarray(uv[0]), asarray(uv[1])))

def setupuv(rc):
    """
    Horn Schunck legacy OpenCV function requires we use these old-fashioned cv matrices, not numpy array
    """
    (r,c) = rc
    u = cv.CreateMat(r, c, cv.CV_32FC1)
    v = cv.CreateMat(r, c, cv.CV_32FC1)
    return (u, v)

def calcofhs(new,ref,smoothing):
    uv = setupuv(new.shape)
    return optflowHornSchunk(new,ref,uv,smoothing)

if __name__ == '__main__':
    data_folder = 'E:\\Cloud\Google Drive\\optical flow\\code\\matlab\\Pre\\'
    save_folder = 'E:\\Cloud\Google Drive\\optical flow\\code\\matlab\\pre_cvOFHS\\'
    for i in xrange(1):
        filename0 = str(i+1) + '.jpg'
        filename1 = str(i+2) + '.jpg'
        img0 = cv2.imread(data_folder + filename0)
        img1 = cv2.imread(data_folder + filename1)

        # convert to gray color image
        img0 = cv2.cvtColor(img0, cv2.COLOR_BGR2GRAY)
        img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2GRAY)
        
        imgofhs = calcofhs(img0,img1,0.05)
        horz = imgofhs[...,0]
        vert = imgofhs[...,1]
        # horz = cv2.normalize(imgofhs[...,0], None, -255, 255, cv2.NORM_MINMAX)     
        # vert = cv2.normalize(imgofhs[...,1], None, -255, 255, cv2.NORM_MINMAX)
        # horz = horz.astype('uint8')
        # vert = vert.astype('uint8')

        horzNP = np.asarray(horz)
        vertNP = np.asarray(vert)

        w,h = horzNP.shape
        x = np.arange(0, w, 1)
        y = np.arange(-h, 0, 1)

        # X, Y = np.meshgrid(np.linspace(0, w, 10), np.linspace(0, h, 10))
        fu = interpolate.interp2d(np.arange(0, w, 1), np.arange(0, h, 1), horzNP,kind='cubic');
        fv = interpolate.interp2d(np.arange(0, w, 1), np.arange(0, h, 1), vertNP,kind='cubic');
        # # v = interpolate.interp2d(vert, X, Y)
        # # f = interpolate.interp2d(x, y, z, kind='cubic')
        u = fu(np.arange(0, w, 30), np.arange(0, h, 30))
        v = fv(np.arange(0, w, 30), np.arange(0, h, 30))

        # x, y = np.mgrid[0:h:500j, 0:w:500j]

        mag = np.sqrt(horz**2 + vert**2)

        cv2.imshow('HS method',mag)
        # cv2.imshow('velocity magnitude',vert)
        plt.figure()
        img = plt.imshow(mag)
        # plt.quiver(np.arange(0, w, 30), np.arange(0, h,30), u, v, scale=10,pivot='mid', color='r')
        plt.axis('off')

        plt.savefig(save_folder+ str(i) + '.jpg',bbox_inches='tight',pad_inches = 0)
        plt.show()

        # cv2.imshow('HS method',horz)

        # k = cv2.waitKey(0)
        # cv2.imwrite(save_folder+ str(i) + '_u.jpg',horz)
        # cv2.imwrite(save_folder+ str(i) + '_v.jpg',vert)
