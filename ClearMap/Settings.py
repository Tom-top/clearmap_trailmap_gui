# -*- coding: utf-8 -*-
"""
Module to set *ClearMap's* internal parameter and paths to external programs.

Notes:
    Edit the :func:`setup` routine to point to the ilastik and elastix paths 
    for specific hosts

See Also:
    * :const:`IlastikPath`
    * :const:`ElastixPath`
    * :mod:`~ClearMap.Parameter`
"""

import os
import socket

###############################################################################
#%% ClearMap paths
###############################################################################

def clearMapPath():
    """Returns root path to the ClearMap software
    
    Returns:
        str: root path to ClearMap
    """
    fn = os.path.split(__file__)
    fn = os.path.abspath(fn[0]);
    return fn;

ClearMapPath = clearMapPath();
"""str: Absolute path to the ClearMap root folder"""

ResourcesPath = os.path.join(clearMapPath(), 'Resources');
"""str: Absolute path to the ClearMap resources folder"""

ExternalPath = os.path.join(clearMapPath(), 'External');
"""str: Absolute path to the ClearMap external program folder"""


###############################################################################
#%% Paths to external programs and resources
###############################################################################

IlastikPath = None;
"""str: Absolute path to the Ilastik installation

Notes:
   `Ilastik Webpage <http://ilastik.org/>`_
   
   `Ilastik Download <http://old.ilastik.org/>`_
"""

#path to eastix installation
ElastixPath = os.path.join(ClearMapPath,"External/Elastix");
"""str: Absolue path to the elastix installation

Notes:
    `Elastix Webpage <http://elastix.isi.uu.nl/>`_
"""


#path to ImageJ/Fiji installation
ImageJPath = None;
"""str: Absolue path to the ImageJ/Fiji installation

Notes:
    `ImageJ/Fiji Webpage <https://fiji.sc/>`_
"""


#path to TeraSticher installation
#TeraStitcherPath = '/home/thomas.topilko/Programs/TeraStitcher-standalone-1.9.63-Linux';
#TeraStitcherPath = '/home/nicolas.renier/Programs/TeraStitcher-custom';
TeraStitcherPath = os.path.join(ClearMapPath, 'External/TeraStitcher-custom');
"""str: Absolue path to the TeraStitcher installation

Notes:
    `TeraSticher Webpage <http://abria.github.io/TeraStitcher/>`_
"""


def setup():
    """Setup ClearMap for specific hosts
    
    Notes:
        Edit this routine to include special setttings for specific hosts
        
    See Also:
        :const:`IlastikPath`, :const:`ElastixPath`, :const:`TeraStitcherPath`,
    """
    global IlastikPath, ElastixPath, TeraStitcherPath, ImageJPath
    
    hostname = socket.gethostname();
    
    if hostname == 'ChristophsComputer.rockefeller.edu':  #Christoph's Laptop 
        IlastikPath      = '/home/ckirst/Programs/ilastik';
        ElastixPath      = os.path.join(ExternalPath, 'elastix/build');
        TeraStitcherPath = os.path.join(ExternalPath, 'TeraStitcher'); 
        ImageJPath       = '/home/ckirst/Programs/Fiji.app';
    
    elif hostname == 'mtllab-Ubuntu': #MTL workstation
        IlastikPath      = '/usr/local/ilastik-1.1.9-Linux';
        ElastixPath      = '/usr/local/elastix'; 
        TeraStitcherPath = '/home/mtllab/Programs/TeraStitcher';
        ImageJPath       = None;
    
    elif hostname == 'ChristophsLaptop.rockefeller.edu':  #Christoph's Desktop 
        IlastikPath      = '/home/ckirst/Programs/ilastik-1.1.9-Linux';
        ElastixPath      = os.path.join(ExternalPath, 'elastix/build');
        TeraStitcherPath = os.path.join(ExternalPath, 'TeraStitcher'); 
        ImageJPath       = '/home/ckirst/Programs/Fiji.app';
      
    
    
    ## insert your hostname specific settings here ##
    #elif hostname == 'your-host-name':
    #    IlastikPath = 'path-to-ilastik';
    #    ElastixPath = 'path-to-elastix';   
    ##

    # check existence:
    if not ElastixPath is None:
        if not os.path.exists(ElastixPath):
            #raise RuntimeWarning('Settings: elastix path %s does not exists, cf. Settings.py or type help(Settings) for details.' % ElastixPath);
            print('Settings: elastix path %s does not exists, cf. Settings.py or type help(Settings) for details.' % ElastixPath);
            ElastixPath = None;
    
    if not IlastikPath is None:
        if not os.path.exists(IlastikPath):
            #raise RuntimeWarning('Settings: ilastik path %s does not exists, cf. Settings.py or type help(Settings) for details.' % IlastikPath);
            print('Settings: ilastik path %s does not exists, cf. Settings.py or type help(Settings) for details.' % IlastikPath);
            IlastikPath = None;
            
    if not TeraStitcherPath is None:
        if not os.path.exists(TeraStitcherPath):
            #raise RuntimeWarning('Settings: TeraStitcher path %s does not exists, cf. Settings.py or type help(Settings) for details.' % TeraStitcherPath);
            print('Settings: TeraStitcher path %s does not exists, cf. Settings.py or type help(Settings) for details.' % TeraStitcherPath);
            TeraStitcherPath = None;
            
    if not ImageJPath is None:
        if not os.path.exists(ImageJPath):
            #raise RuntimeWarning('Settings: TeraStitcher path %s does not exists, cf. Settings.py or type help(Settings) for details.' % TeraStitcherPath);
            print('Settings: ImageJ path %s does not exists, cf. Settings.py or type help(Settings) for details.' % ImageJPath);
            ImageJPath = None;

setup();



