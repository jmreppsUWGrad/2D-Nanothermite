
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 02 13:32:07 2018

@author: Joseph

This contains classes for reading and writing files in good format

"""

keys_Settings=['Length','Width','Nodes_x','Nodes_y','k','Cp','rho',\
               'bias_type_x','bias_size_x','bias_type_y','bias_size_y']
               
keys_Sources=['Source_Uniform','Source_Kim','Ea','A0','dH']
keys_Time_adv=['Fo','dt','total_time_steps', 'total_time','Time_Scheme',\
               'Convergence','Max_iterations','Output_directory','Number_Data_Output']
keys_BCs=     ['bc_left','bc_left_rad',\
              'bc_right','bc_right_rad',\
              'bc_south','bc_south_rad',\
              'bc_north','bc_north_rad']

import string as st

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
        
        self.Write_single_line('\nSource Terms:')
        keys=['Source_Uniform','Source_Kim','Ea','A0','dH']
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
        keys=['Fo','dt','total_time_steps', 'Time_Scheme','Convergence',\
              'Max_iterations','Output_directory']
        for i in keys:
            self.fout.write(i)
            self.fout.write(':')
            self.Write_single_line(str(settings[i]))
#            self.fout.write('\n')
        
        self.Write_single_line('\nBoundary conditions:')
        for i in keys_BCs:
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
        self.fin=open(filename+'.txt', read_type)
        
    def Read_Input(self, settings, Sources, BCs):
        for line in self.fin:
            if st.find(line, ':')>0 and st.find(line, '#')!=0:
                line=st.split(line, ':')
                # Domain settings
                if line[0] in keys_Settings:
                    if line[0]=='Nodes_x' or line[0]=='Nodes_y':
                        settings[line[0]]=int(line[1])
                    elif line[1]=='None\n':
                        settings[line[0]]=st.split(line[1], '\n')[0]
                    else:
                        settings[line[0]]=float(line[1])
                # Source term info
                elif line[0] in keys_Sources:
                    if line[1]=='None\n' or line[1]=='True\n':
                        Sources[line[0]]=st.split(line[1], '\n')[0]
                    else:
                        Sources[line[0]]=float(line[1])
                # Time advancement details
                elif line[0] in keys_Time_adv:
                    if line[0]=='Time_Scheme' or line[1]=='None\n':
                        settings[line[0]]=st.split(line[1], '\n')[0]
                    elif line[0]=='total_time_steps' or line[0]=='Max_iterations'\
                        or line[0]=='Number_Data_Output':
                        settings[line[0]]=int(line[1])
                    elif line[0]=='Output_directory':
                        settings[line[0]]=line[1]+':'+st.split(line[2], '\n')[0]
                    else:
                        settings[line[0]]=float(line[1])
                # Boundary conditions
                elif line[0] in keys_BCs:
                    BC_info=st.split(line[1], ',')
                    BCs[line[0]]=[]
                    # Radiation BCs
                    if line[0]=='bc_left_rad' or line[0]=='bc_right_rad'\
                        or line[0]=='bc_north_rad' or line[0]=='bc_south_rad':
                        try:
                            BCs[line[0]]=[float(BC_info[0])]
                            BCs[line[0]]+=[float(BC_info[1])]
                        except:
                            BCs[line[0]]=st.split(BC_info[0], '\n')[0]
                        del BC_info[0]
                    # All other BCs
                    i=0
                    while len(BC_info)>1:
                        BCs[line[0]]+=[BC_info[0]]
                        del BC_info[0]
                        # Temp or flux BC
                        if BCs[line[0]][3*i]=='T' or BCs[line[0]][3*i]=='F':
                            # Value into a float
                            BCs[line[0]]+=[float(BC_info[0])]
                            del BC_info[0]
                            BCs[line[0]]+=[(int(BC_info[0]),int(BC_info[1]))]
                            del BC_info[1], BC_info[0]
                        # Convective BC
                        elif BCs[line[0]][3*i]=='C':
                            BCs[line[0]]+=[(float(BC_info[0]),float(BC_info[1]))]
                            del BC_info[1], BC_info[0]
                            BCs[line[0]]+=[(int(BC_info[0]),int(BC_info[1]))]
                            del BC_info[1], BC_info[0]
                        # Radiation BC
                        else:
                            try:
                                BCs[line[0]]=[float(BCs[line[0]][0])]
                                BCs[line[0]]+=[float(BC_info[0])]
                            except:
                                continue
                            del BC_info[0]
                        i+=1
                    
        self.fin.close()