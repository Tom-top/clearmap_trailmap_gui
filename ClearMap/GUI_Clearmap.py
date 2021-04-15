#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 18:34:51 2018

@author: tomtop
"""
		
import sys, os;
from natsort import natsorted;
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton,\
                            QGroupBox, QDialog,\
                            QVBoxLayout, QGridLayout, QLabel,\
                            QLineEdit, QTabWidget;
from PyQt5.QtCore import QCoreApplication;
import ClearMap.Tomek_Utilities as ut;
import Parameters_Raw as p;

class createObject() :
    
    def __init__(self) :
        pass;
        
class selectionMenu(QDialog):
    
   def __init__(self, folders, parent=None):
      super(selectionMenu, self).__init__(parent)
      
      self.title = 'Clearmap GUI v1.0'
      self.right = 10;
      self.left = 10;
      self.width = 320;
      self.height = 100;
      self.toggle = None;
      
      self.folders = folders;
      self.number = len(self.folders);
      
      self.buttons2Create = ["Stiching","Resamp_Auto","Resamp_cFOS",\
                            "Align_Auto","Align_Template","Cell_Detect"\
                            ,"Heatmap","cFOS_for_Annotation","all","none"];
                             
      self.textboxes2Create = ["orient","bgthresh","csthresh","cfosX","cfosY","cfosZ"];
      
      self.initUI();
      
   def initUI(self) :
       
      self.setWindowTitle(self.title);
      self.setGeometry(self.left,self.right,self.width,self.height);
      
      self.mainLayout = QVBoxLayout();
      
      self.tabs = QTabWidget();
      
      self.tab1 = QWidget();
      self.tabs.addTab(self.tab1,"Operations");
      
      self.createTab1Layout();
      self.tab1.layout = QVBoxLayout();
      self.tab1.layout.addWidget(self.GroupBoxTab1);
      self.tab1.setLayout(self.tab1.layout)
      
      self.tab2 = QWidget();
      self.tabs.addTab(self.tab2,"Parameters");
      
      self.createTab2Layout();
      self.tab2.layout = QVBoxLayout();
      self.tab2.layout.addWidget(self.GroupBoxTab2);
      self.tab2.setLayout(self.tab2.layout)
      
      self.mainLayout.addWidget(self.tabs);
      self.setLayout(self.mainLayout);
      
      self.show();
      
   def createTab1Layout(self) :
       
      self.GroupBoxTab1 = QGroupBox("Folders in the Working Directory");
      layout = QGridLayout();
      
      self.textTab1 = {};
      self.buttons = {};
      
      self.checked = {};
      
      y = 0;
      
      for n,folder in natsorted(enumerate(self.folders)) :

          self.textTab1[folder] = createObject();
          self.textTab1[folder].text = QLabel(folder);
          layout.addWidget(self.textTab1[folder].text, n, y);
          self.buttons[folder] = createObject();
          self.buttons[folder].row = n;
          
          for button in self.buttons2Create :
              
              setattr(self.buttons[folder], button, QPushButton(button));
              attr = getattr(self.buttons[folder], button);
              attr.setAutoDefault(False);
              attr.setCheckable(True);
              
              if button == 'all' :
                  attr.clicked.connect(lambda state, x=folder, y=button\
                                       : self.toggleAllHandler(x,y));
                  
              elif button == 'none' :
                  attr.clicked.connect(lambda state, x=folder, y=button\
                                       : self.toggleAllHandler(x,y));
              
              layout.addWidget(attr, n, y+1);
              
              y+=1;
              
          y = 0;
          
      self.ApplyButton = QPushButton("Apply");
      self.ApplyButton.setAutoDefault(False);
      self.ApplyButton.setCheckable(True);
      layout.addWidget(self.ApplyButton, self.number+1, len(self.buttons2Create));
      self.ApplyButton.clicked.connect(self.applyChanges);
      
      self.QuitButton = QPushButton("Quit");
      self.QuitButton.setAutoDefault(False);
      self.QuitButton.setCheckable(True);
      layout.addWidget(self.QuitButton, self.number+1, len(self.buttons2Create)-1);
      self.QuitButton.clicked.connect(self.quitApplication);
          
      self.GroupBoxTab1.setLayout(layout);
      
   def createTab2Layout(self) :
       
      self.GroupBoxTab2 = QGroupBox("Folders in the Working Directory");
      layout = QGridLayout();
      
      self.textTab2 = {};
      self.textbox = {};
      
      self.bgthresholds = {};
      self.csthresholds = {};
      self.orients = {};
      self.cfosX = {};
      self.cfosY = {};
      self.cfosZ = {};
      
      for n,t in enumerate(self.textboxes2Create) :
          
          label = QLabel(t);
          
          layout.addWidget(label, 0, n+1);
          layout.setColumnStretch(0, 1);
      
      y = 0;
      
      for n,folder in natsorted(enumerate(self.folders)) :

          self.textTab2[folder] = createObject();
          self.textTab2[folder].text = QLabel(folder);
          layout.addWidget(self.textTab2[folder].text, n+1, y);
          self.textbox[folder] = createObject();
          self.textbox[folder].row = n+1;
          
          for textbox in self.textboxes2Create :
              
              setattr(self.textbox[folder], textbox, QLineEdit());
              attr = getattr(self.textbox[folder], textbox);
              
              if textbox == "bgthresh" :
                  
                  attr.setFixedWidth(40);
                  attr.setText("{0}".format(p.detectCellShapeParameter["threshold"]));
                  attr.setToolTip("Background Threshold");
                  
              elif textbox == "csthresh" :
                  
                  attr.setFixedWidth(60);
                  attr.setText("{0}".format(p.thresholdPointParameter["threshold"]));
                  attr.setToolTip("Cell Size Threshold");
                  
              elif textbox == "orient" :
                  
                  attr.setFixedWidth(60);
                  attr.setText("{0}".format(p.FinalOrientation));
                  attr.setToolTip("Brain Orientation");
                  
              elif textbox == "cfosX" :
                  
                  attr.setFixedWidth(100);
                  attr.setText("{0}".format(p.cFosDetectionRange["x"]));
                  attr.setToolTip("cFOS X");
                  
              elif textbox == "cfosY" :
                  
                  attr.setFixedWidth(100);
                  attr.setText("{0}".format(p.cFosDetectionRange["y"]));
                  attr.setToolTip("cFOS Y");
                  
              elif textbox == "cfosZ" :
                  
                  attr.setFixedWidth(100);
                  attr.setText("{0}".format(p.cFosDetectionRange["z"]));
                  attr.setToolTip("cFOS Z");

              layout.addWidget(attr, n+1, y+1);
              
              y+=1;
              
          y = 0;
          
          layout.setColumnStretch(n+1, 1);
     
      self.GroupBoxTab2.setAlignment(0);
      self.GroupBoxTab2.setLayout(layout);
      
   def toggleAllHandler(self,folder,button) :
       
      attr = getattr(self.buttons[folder],button);
      row = self.buttons[folder].row;
       
      if attr.isChecked():
          
          if button == "all" :
              self.toggle = True;
          elif button == "none" :
              self.toggle = False;
          self.toggleAll(row);
          
      else:
          pass;
          
   def toggleAll(self,row) :
       
      for key in self.buttons.keys() :
          if self.buttons[key].row == row :
              for button in self.buttons2Create :
                  attr = getattr(self.buttons[key],button);
                  if button != "all" and button != "none" :
                       attr.setChecked(self.toggle);
                  if button == "all" or button == "none" :
                       attr.setChecked(False);
                       
   def applyChanges(self) :
      
      if self.ApplyButton.isChecked() :
          
          for key,value in zip(self.textbox.keys(),self.textbox.values()) :
              
              attr = getattr(value,"bgthresh");
              bgthresh = attr.text();
              self.bgthresholds[key] = eval(bgthresh);
              
              attr = getattr(value,"csthresh");
              csthresh = attr.text();
              self.csthresholds[key] = eval(csthresh);
              
              attr = getattr(value,"orient");
              orient = attr.text();
              self.orients[key] = eval(orient);
              
              attr = getattr(value,"cfosX");
              cfosX = attr.text();
              if cfosX == "<built-in function all>" :
                  self.cfosX[key] = "all";
              else :
                  self.cfosX[key] = eval(cfosX);
              
              attr = getattr(value,"cfosY");
              cfosY = attr.text();
              if cfosY == "<built-in function all>" :
                  self.cfosY[key] = "all";
              else :
                  self.cfosY[key] = eval(cfosY);
              
              attr = getattr(value,"cfosZ");
              cfosZ = attr.text();
              if cfosZ == "<built-in function all>" :
                  self.cfosZ[key] = "all";
              else :
                  self.cfosZ[key] = eval(cfosZ);
          
          for key,value in zip(self.buttons.keys(),self.buttons.values()) :
              temp = "";
              for buttons in self.buttons2Create :
                  if buttons != "all" and buttons != "none" and buttons != "thresh" :
                      attr = getattr(value,buttons);
                      if attr.isChecked() :
                          temp+="T";
                      else :
                          temp+="F";       
              self.checked[key] = temp;
              temp = "";
              
          QCoreApplication.instance().quit();
          
   def quitApplication(self) : 
        
      if self.QuitButton.isChecked() :
          QCoreApplication.instance().quit();
      
def main(folders):
    
    if not QApplication.instance():
        app = QApplication(sys.argv)
    else:
        app = QApplication.instance() 
    
    GUI = selectionMenu(folders);
    
    clickedFolders = GUI.checked;
    bgthresholds = GUI.bgthresholds;
    csthresholds = GUI.csthresholds;
    orients = GUI.orients;
    cfosX = GUI.cfosX;
    cfosY = GUI.cfosY;
    cfosZ = GUI.cfosZ;
    
    parameters = {
            "operations" : clickedFolders,
            "bgthresholds" : bgthresholds,
            "csthresholds" : csthresholds,
            "orients" : orients,
            "cfosX" : cfosX,
            "cfosY" : cfosY,
            "cfosZ" : cfosZ,
            };
    
    return app.exec_(), parameters;
#    sys.exit(app.exec_());
    
def launchGUI(fd,experiment) :
  
    print(ut.coloredMessage('[INFO] Automatic mode selected, GUI_Clearmap used','darkgreen'));

    folders = [];
  
    for folder in natsorted(os.listdir(fd)) :
        if folder.split('-')[0] == experiment or folder.split('_')[0] == experiment :
            if os.path.isdir(os.path.join(fd,folder)) :
                folders.append(folder);
        
  

    _,parameters = main(folders);
  
    return parameters;

#if __name__ == '__main__':
#    main();