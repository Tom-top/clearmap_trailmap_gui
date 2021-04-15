#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 18:34:51 2018

@author: tomtop
"""
		
import sys, os;
import numpy as np;
from natsort import natsorted;
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton,\
                            QHBoxLayout, QGroupBox, QDialog,\
                            QVBoxLayout, QGridLayout, QLabel;
from PyQt5.QtGui import QIcon;
from PyQt5.QtCore import pyqtSlot, QCoreApplication;

class createObject() :
    
    def __init__(self) :
        pass;
        
class selectionMenu(QDialog):
    
   def __init__(self, folders, experiment, allGroups, parent=None):
      super(selectionMenu, self).__init__(parent)
      
      self.title = 'Analysis GUI v1.0'
      self.right = 10;
      self.left = 10;
      self.width = 320;
      self.height = 100;
      self.toggle = None;
      self.experiment = experiment;
      self.groups = allGroups;
      
      flds = [];
      
      for fld in folders :
        
        if fld.split('_')[0] != fld :
          split = fld.split('_');
          exp = fld.split('_')[0];
        elif fld.split('-')[0] != fld :
          split = fld.split('-');
          exp = fld.split('-')[0];
        else :
          raise RuntimeError('No folder(s) with the name {} were found in {} \n'.format(experiment,homedir));
          
        if exp == self.experiment :
          flds.append(split[-1]);
          
      flds.append('all');
      flds.append('none');
      
      self.buttons2Create = flds;
      self.number = len(self.buttons2Create);
      
      self.initUI();
      
   def initUI(self) :
       
      self.setWindowTitle(self.title);
      self.setGeometry(self.left,self.right,self.width,self.height);
      
      self.createGridLayout();
      
      windowLayout = QVBoxLayout();
      windowLayout.addWidget(self.horizontalGroupBox);
      self.setLayout(windowLayout);
      
      self.show();
      
   def createGridLayout(self) :
       
      self.horizontalGroupBox = QGroupBox("Animals to be analyzed");
      layout = QGridLayout();
      
      self.text = {};
      self.buttons = {};
      self.checked = {};
      y = 0;
      
      for n,group in natsorted(enumerate(self.groups)) :

          self.text[group] = createObject();
          self.text[group].text = QLabel(group);
          layout.addWidget(self.text[group].text, n, y);
          self.buttons[group] = createObject();
          self.buttons[group].row = n;
          
          for button in self.buttons2Create :
              
              setattr(self.buttons[group], button, QPushButton(button));
              attr = getattr(self.buttons[group], button);
              attr.setAutoDefault(False);
              attr.setCheckable(True);
              
              if button == 'all' :
                  attr.clicked.connect(lambda state, x=group, y=button\
                                       : self.toggleAllHandler(x,y));
                  
              elif button == 'none' :
                  attr.clicked.connect(lambda state, x=group, y=button\
                                       : self.toggleAllHandler(x,y));
              
              layout.addWidget(attr, n, y+1);
              
              y+=1;
              
          y = 0;
          
      self.ApplyButton = QPushButton("Apply");
      self.ApplyButton.setAutoDefault(False);
      self.ApplyButton.setCheckable(True);
      layout.addWidget(self.ApplyButton, self.number+1, self.number);
      self.ApplyButton.clicked.connect(self.applyChanges);
      
      self.QuitButton = QPushButton("Quit");
      self.QuitButton.setAutoDefault(False);
      self.QuitButton.setCheckable(True);
      layout.addWidget(self.QuitButton, self.number+1, self.number-1);
      self.QuitButton.clicked.connect(self.quitApplication);
          
      self.horizontalGroupBox.setLayout(layout);
      
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
          for key,value in zip(self.buttons.keys(),self.buttons.values()) :
              temp = "";
              for buttons in self.buttons2Create :
                  if buttons != "all" and buttons != "none" :
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
      
def main(folders,experiment,allGroups):
  
   app = QApplication(sys.argv);
   GUI = selectionMenu(folders,experiment,allGroups);
   clickedFolders = GUI.checked;
   return app.exec_(),clickedFolders;
   #sys.exit(app.exec_());

if __name__ == '__main__':
   main();
