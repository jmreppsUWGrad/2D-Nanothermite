
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 02 13:32:07 2018

@author: Joseph

This contains classes for reading and writing files in good format

"""

class FileOut():
    def __init__(self, filename, isBin):
        self.name=filename
        if isBin:
            write_type='wb'
        else:
            write_type='w'
        self.fout=open(filename+'.txt', write_type)
    
    # Write a single string with \n at end
    def Write_single_line(self, string):
        self.fout.write(string)
        self.fout.write('\n')
    
    # Header with information about file
    def header_cond(self, title='Run'):
        self.Write_single_line('######################################################')
        self.Write_single_line('#             2D Heat Conduction Solver              #')
        self.Write_single_line('#              Created by J. Mark Epps               #')
        self.Write_single_line('#          Part of Masters Thesis at UW 2018-2020    #')
        self.Write_single_line('######################################################\n')
        self.Write_single_line('############### '+title+' FILE #########################')
        self.Write_single_line('##########'+self.name+'##################\n')
    
    def input_writer_cond(self, settings, BCs, T):
        self.Write_single_line('Settings:')
        keys=['Length','Width','Nodes_x','Nodes_y','k','Cp','rho']
        for i in keys:
            self.fout.write(i)
            self.fout.write(':')
            self.Write_single_line(str(settings[i]))
#            self.fout.write('\n')
        
        self.Write_single_line('\nSource Terms:\n')
        keys=['Source_Uniform']
        for i in keys:
            self.fout.write(i)
            self.fout.write(':')
            self.Write_single_line(str(settings[i]))

        self.Write_single_line('\nMeshing details:')
        keys=['bias_type_x','bias_size_x','bias_type_y','bias_size_y']
        for i in keys:
            self.fout.write(i)
            self.fout.write(':')
            self.Write_single_line(str(settings[i]))
#            self.fout.write('\n')
        
        self.Write_single_line('\nTime advancement:')
        keys=['Fo','total_time_steps', 'Time_Scheme','Convergence','Max_iterations']
        for i in keys:
            self.fout.write(i)
            self.fout.write(':')
            self.Write_single_line(str(settings[i]))
#            self.fout.write('\n')
        
        self.Write_single_line('\nBoundary conditions:')
        keys=['bc_left','bc_left_rad',\
              'bc_right','bc_right_rad',\
              'bc_south','bc_south_rad',\
              'bc_north','bc_north_rad']
        for i in keys:
            self.fout.write(i)
            self.fout.write(':')
            self.Write_single_line(str(BCs[i]))
#            self.fout.write('\n')
        
        self.fout.write('\nInitial conditions:\n')
        self.Write_single_line('T')
        for i in range(len(T[:,0])):
            self.Write_single_line(str(T[i,:]))
        
        self.fout.write('\n')
        
    def Write_timestep_data(self, timeStep, dt):
        self.fout.write('Time step: '+timeStep+'\n')
        self.fout.write('Time step size: '+dt+'\n\n')
        
    
    def close(self):
        self.fout.close()
        
class FileIn():
    def __init__(self, filename, isBin):
        self.name=filename
        if isBin:
            read_type='rb'
        else:
            read_type='r'
        self.fout=open(filename+'.txt', read_type)