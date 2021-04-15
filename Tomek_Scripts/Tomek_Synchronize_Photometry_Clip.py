#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 16:21:34 2019

@author: thomas.topilko
"""

import os;
import pandas as pd;
import matplotlib.pyplot as plt;
import matplotlib.animation as animation;
import numpy as np;

import moviepy.editor as mpy;
from moviepy.video.io.bindings import mplfig_to_npimage;
from moviepy.editor import VideoFileClip, VideoClip, clips_array,ImageSequenceClip;

data = "/home/thomas.topilko/Desktop/D1default_2.csv";
videoFile = "/home/thomas.topilko/Desktop/Raw_Video_Mice_808_12-11-2019_15-15-10.avi";

sheet = pd.read_csv(data, header=None, usecols=np.arange(0,9));

samplingRate = 12000;
resolution = 1;

X = np.array([ np.mean( [  float(i) for i in sheet[0][1:][ int(n) : int(n+(samplingRate*resolution)) ] ] ) for n in np.arange(0, len(sheet[0][1:]), (samplingRate*resolution)) ]);
Y = np.array([ np.mean( [  float(i) for i in sheet[1][1:][ int(n) : int(n+(samplingRate*resolution)) ] ] ) for n in np.arange(0, len(sheet[1][1:]), (samplingRate*resolution)) ]);

dataDuration = 200;
Xnew = X[:dataDuration];
Ynew = Y[:dataDuration];

vidClip = mpy.VideoFileClip(videoFile).subclip(t_start=0, t_end=dataDuration).crop(x1=931,y1=153,x2=1475,y2=410);
threshDisplay = int(2*vidClip.fps);

def LivePhotometryTrack(vidClip, X, Y, acceleration=1, showHeight=True, showTrace=True) :
    
    global threshDisplay
    
    def make_frame_mpl(t):
                
        i = int(t);
        
        if i < threshDisplay :
            
            try :
            
                trajectoryGraph.set_data(X[0:i],Y[0:i]);
#                trajectoryGraph.set_data(list(zip(*Y[0:i]))[0],list(zip(*Y[0:i]))[1]);
                
            except :
                
                pass;
    
            last_frame = mplfig_to_npimage(liveFig);
            return last_frame;
        
        else :
            
            delta = i - threshDisplay;
        
            liveAx2.set_xlim(delta,i);
            
            try :
            
                trajectoryGraph.set_data(X[0:i],Y[0:i]);
#                trajectoryGraph.set_data(list(zip(*Y[0:i]))[0],list(zip(*Y[0:i]))[1]);
                
            except :
                
                pass;
    
            last_frame = mplfig_to_npimage(liveFig);
            return last_frame;
    
    _FrameRate = vidClip.fps;
    _Duration = vidClip.duration;
    
    liveFig = plt.figure(figsize=(5,3), facecolor='white');
    
    gs = liveFig.add_gridspec(2, 2);
    
    liveAx2 = liveFig.add_subplot(gs[0:, 0:]);
    
    liveAx2.set_title("mV (Convoluted signal 465nm+450nm) over time");
    liveAx2.set_xlim([0, threshDisplay]);
    liveAx2.set_ylim([min(Y), max(Y)]);
#    liveAx2.set_aspect("equal");

    trajectoryGraph, = liveAx2.plot(X[0],Y[0],'-o',color="blue",alpha=0.8,ms=1.);
    
#    plt.gca().invert_yaxis();
#        plt.gca().invert_xaxis();
    plt.tight_layout();
    
    _acceleration = acceleration;
    
    anim = mpy.VideoClip(make_frame_mpl, duration=(_Duration*1));
    
    finalClip = clips_array([[clip.margin(2, color=[255,255,255]) for clip in
                    [(vidClip.resize(1.).speedx(_acceleration)), anim.speedx(_acceleration)]]],
                    bg_color=[255,255,255]); #.speedx(self._FrameRate)
    
    finalClip.write_videofile(os.path.join("/home/thomas.topilko/Desktop",'PhotoMetry_Tracking.mp4'), fps=10);
    
LivePhotometryTrack(vidClip, Xnew, Ynew, acceleration=5)


        
