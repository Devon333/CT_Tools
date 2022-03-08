#**************************************************************************
# Molecule class Reads ADF output, organizes MOs, organizes excited state 
# info by symmetry 
# Alva Dillon, Hanwen Liu Spring 2022
# 
#**************************************************************************
import sys

class Molecule:
    filename = None
    Orbital_Occ = {} # orbital label: occupation
    Orbital_Energy = {} # orbital label: MO Energy
    Orbital_Localized_Character = {} # orbital label: [percent metal , percent molecule]
    Orbital_Character = {} # orbital label:[s character, p character, d character]
    metal_type= None
    Excited_States= {} # {Symm:{exc_num: [energy, osc_str, tdm_x, tdm_y, tdm_z]} } 
    Excited_State_Decomp = {} # {Symm: {Exc_num:{transition label: [weight, tdm_x, tdm_y, tdm_z]} } }

    def __init__(self, ADF_Outname, metal="Ag"):
        self.filename = str(ADF_Outname)
        self.metal_type = metal 
    

    def make_latex_table_header(self, outfile):
        header='''\\documentclass[10pt,a4paper]{article}
\\usepackage[utf8]{inputenc}
\\usepackage[T1]{fontenc}
\\usepackage{subfig}
\\usepackage{amsmath}
\\usepackage{float}
\\usepackage{amsfonts}
\\usepackage{amssymb}
\\usepackage{graphicx}
\\usepackage{ragged2e}
\\usepackage[margin=.75in]{geometry}
\\usepackage[english]{babel}
\\usepackage{comment}
\\usepackage{longtable}
\\begin{document}
%\\begin{table}[H]
 \\centering
 \\begin{longtable}{c|c|c|c|c|c}
   \\caption{Excited States and Electron Configurations Table}
   \\hline
     Energy & Osc. Str. & Transition & Weight & Transition Dipole Moment(x,y,z) & Orb Char (s,p,d)\\\\
    \\hline \n 
   \\endfirsthead 
   \\hline
   \\endhead
   \\hline
   \\endfoot
   \\endlastfoot 
'''
        out = open(outfile,'w')
        out.write(header)
        out.close()

    def make_latex_table_end(self,outfile):
        end=""" %\hline \end{tabular} \n \end{longtable} \n \end{document} """
        out = open(outfile,'a')
        out.write(end)
        out.close()

    def fill_orbital_occ(self,split_line):
        orbital = str(split_line[2])+str(split_line[3])
        orbital = orbital.split(":")
        orbital = orbital[0].lower()
        self.Orbital_Occ[orbital] = str(split_line[1])

    def fill_orbital_energy(self,split_line):
        orbital = str(split_line[2])+str(split_line[3])
        orbital = orbital.split(":")
        orbital = orbital[0].lower()
        self.Orbital_Energy[orbital] = str(split_line[0])

    def fill_exc_states(self, symmetry, split_line):
        self.Excited_States[symmetry][split_line[0]] = [ float(split_line[1]), float(split_line[2]), float(split_line[3]), float(split_line[4]), float(split_line[5]) ]
     
    def fill_exc_state_decomp(self, symmetry, split_line):
        state_num = split_line[0].split(":")
        state_num = state_num[0]
        transition =split_line[1:4]
        transition = f"{transition[0]} $\\rightarrow$ {transition[2]}"
        #self.Excited_State_Decomp[symmetry][state_num] = {}
        self.Excited_State_Decomp[symmetry][state_num][transition] = [ float(split_line[4]), float(split_line[5]), float(split_line[6]), float(split_line[7]) ]
        #print(self.Excited_State_Decomp)
 
    def character_type_sum(self,split_line):
        char=[0,0,0]
        if len(split_line) == 11:
            percent = split_line[4]
            if "D" in split_line[6]:
                char[2]+= float(percent[:-1])
            if "P" in split_line[6]:
                char[1] += float(percent[:-1])
            if "S" in split_line[6]:
                char[0] += float(percent[:-1])
        if len(split_line) == 7:
            percent = split_line[0]
            if "D" in split_line[2]:
                char[2] += float(percent[:-1])
            if "P" in split_line[2]:
                char[1] += float(percent[:-1])
            if "S" in split_line[2]:
                char[0] += float(percent[:-1]) 
        return char


    def print_by_transition_dipole_moment(self, symmetry, component_indx, component_cutoff, outputfilename):
        states = self.Excited_States
        state_decomp = self.Excited_State_Decomp
        orbs = self.Orbital_Character
        self.make_latex_table_header(outputfilename)
        output = open(outputfilename,'a')
        for exc in states[symmetry]:
            if abs(float(states[symmetry][exc][2+component_indx])) > component_cutoff:
                output.write(f"{states[symmetry][exc][0]:.2f} &\t {states[symmetry][exc][1]:.2f} &\t &\t &\t {states[symmetry][exc][2]:.2f}, {states[symmetry][exc][3]:.2f}, {states[symmetry][exc][4]:.2f} &\t   \\\\ \n  ")
                for tran in state_decomp[symmetry][exc]:
                    t_split = tran.split("$\\rightarrow$")
                    occ = t_split[0].replace(" ","")
                    unocc = t_split[1].replace(" ","")
                    #print(orbs) 
                    #print(f"occupied orb {occ} unoccupied orb {unocc}")
                    try:
                        output.write(f"\t& \t& {tran} &\t {state_decomp[symmetry][exc][tran][0]:.4f} &\t {state_decomp[symmetry][exc][tran][1]:.3f}, {state_decomp[symmetry][exc][tran][2]:.3f}, {state_decomp[symmetry][exc][tran][3]:.3f}  &\t {float(orbs[occ][0]):.2f}, {float(orbs[occ][1]):.2f}, {float(orbs[occ][2]):.2f} $\\rightarrow$ {float(orbs[unocc][0]):.2f}, {float(orbs[unocc][1]):.2f}, {float(orbs[unocc][2]):.2f} \\\\ \n ") 
                    
                    #    output.write(f"\t& \t& {tran} &\t {state_decomp[symmetry][exc][tran][0]} &\t {state_decomp[symmetry][exc][tran][1]}, {state_decomp[symmetry][exc][tran][2]}, {state_decomp[symmetry][exc][tran][3]}  &\t {orbs[occ]} $\\rightarrow$ {orbs[unocc]}  \\\\ \n ") 
                    #try:
                    #    output.write(f"\t& \t& {tran} &\t {state_decomp[symmetry][exc][tran][0]} &\t {state_decomp[symmetry][exc][tran][1]}, {state_decomp[symmetry][exc][tran][2]}, {state_decomp[symmetry][exc][tran][3]}  &\t {orbs[occ]} $\\rightarrow$ ?  \\\\ \n ") 
                    except:
                        output.write(f"\t& \t& {tran} &\t {state_decomp[symmetry][exc][tran][0]:.4f} &\t {state_decomp[symmetry][exc][tran][1]:.3f}, {state_decomp[symmetry][exc][tran][2]:.3f}, {state_decomp[symmetry][exc][tran][3]:.3f}  &\t ? $\\rightarrow$ ?  \\\\ \n ") 
        output.close()
        self.make_latex_table_end(outputfilename)
    def localized_orbital_sum(self, split_line):
        char=[0,0]
        if len(split_line) == 11:
            percent = split_line[4]
            if self.metal_type in split_line[10]:
                char[0] += float(percent[:-1])
            else:
                char[1] += float(percent[:-1])
        if len(split_line) == 7:
            percent = split_line[0]
            if self.metal_type in split_line[6]:
                char[0] += float(percent[:-1])
            else:
                char[1] += float(percent[:-1])
        return char 


    def fill_orbital_characters(self,label, s_percent, p_percent, d_percent):
        self.Orbital_Character[label.lower()] = [s_percent, p_percent, d_percent]

    def fill_localized_characters(self,label,sys1,sys2):
        self.Orbital_Localized_Character[label.lower()]=[sys1,sys2]


    def get_exc_states(self,symmetry):
        ADF_Out = open(self.filename,'r')
        line = ADF_Out.readline()
        self.Excited_States[symmetry]={}
        findline= f"Symmetry {symmetry}"
        findline2="Transition dipole moments mu"
        while line:
            if findline in line:
                while findline2 not in line:
                    line = ADF_Out.readline()    
                for _ in range(5):
                    line = ADF_Out.readline()
                sp_line = line.split()
                #print(line)
                while len(sp_line) > 2:
                    #print(f"{sp_line}")
                    self.fill_exc_states(symmetry,sp_line)
#                    self.Excited_States[symmetry][sp_line[0]] = [ sp_line[1], sp_line[2], sp_line[3], sp_line[4], sp_line[5] ]
                    line = ADF_Out.readline()
                    sp_line = line.split()
#                while findline2 not in line:
#                    line = ADF_Out.readline()    
            line = ADF_Out.readline()
        ADF_Out.close()
        print(f"Excited States \n {self.Excited_States}")


    def get_exc_decomp(self, symmetry):
        ADF_Out = open(self.filename,'r')
        line = ADF_Out.readline()
        self.Excited_State_Decomp[symmetry]={}
        findline= f"Symmetry {symmetry}"
        findline2="Major MO -> MO transitions for the above excitations"
        stopline="Eigenvalues of small (approximate) problem"
        while line:
            if findline in line:
                while findline2 not in line:
                    line = ADF_Out.readline()
                for _ in range(9):
                    line = ADF_Out.readline()
                    sp_line = line.split()
                self.Excited_State_Decomp[symmetry][sp_line[0].split(":")[0]]={}                    
                while stopline not in line:
                    #print(line)
                    sp_line = line.split()
                    if line == "\n":
                        line = ADF_Out.readline()
                        line = ADF_Out.readline()
                        #print(f"printed line{line}")
                        sp_line = line.split()
                        self.Excited_State_Decomp[symmetry][sp_line[0].split(":")[0]]={}                    
                    #print(f"printed line{line}")
                    #print(f"printed split{sp_line}")
                    if sp_line == []:
                        break
                    self.fill_exc_state_decomp(symmetry,sp_line)
                    line = ADF_Out.readline()
            line = ADF_Out.readline()
        ADF_Out.close()
        print(f"Excited state decompositions \n {self.Excited_State_Decomp}")




    def get_MOs(self):
         ADF_Out = open(self.filename,'r')
         line = ADF_Out.readline()
         findline="E(eV)  Occ       MO           %     SFO (first member)   E(eV)  Occ   Fragment"
         while line:
             if findline in line:
                 for _ in range(3):
                     line = ADF_Out.readline()
                     orbital=None
                 while line !="\n": 
                     sp_line = line.split()
                     if len(sp_line) == 11:
                         s_char = 0.0
                         p_char = 0.0
                         d_char = 0.0
                         metal_char = 0.0
                         mol_char = 0.0
                         sp_line = line.split()
                         #print(sp_line," ",len(sp_line))
                         self.fill_orbital_occ(sp_line)
                         self.fill_orbital_energy(sp_line)
                         char_type = self.character_type_sum(sp_line)  
                         loc_type = self.localized_orbital_sum(sp_line)
#                         print(char_type)
                         metal_char += loc_type[0]
                         mol_char += loc_type[1]
                         s_char += char_type[0]
                         p_char += char_type[1]
                         d_char += char_type[2]
                         #print(f"{s_char}, {p_char}, {d_char}") 
                         orbital = str(sp_line[2])+str(sp_line[3])
                         orbital = orbital.split(":")
                         orbital = orbital[0]
                         #orbital= str(sp_line[2])+str(sp_line[3])
                     if len(sp_line) == 7:                         
                         sp_line = line.split()
                         #print(sp_line)
                         char_type = self.character_type_sum(sp_line) 
                         loc_type = self.localized_orbital_sum(sp_line)
#                         print(char_type)
                         s_char += char_type[0]
                         p_char += char_type[1]
                         d_char += char_type[2]
                         metal_char += loc_type[0]
                         mol_char += loc_type[1] 
                         #print(f"{s_char}, {p_char}, {d_char}") 
                         self.fill_orbital_characters(orbital,s_char,p_char,d_char)
                         self.fill_localized_characters(orbital,metal_char,mol_char)                    
                     line = ADF_Out.readline()
             line = ADF_Out.readline()
          
         ADF_Out.close()
         print(f"Orbital Energies \n{self.Orbital_Energy}")
         print(f"Orbital Occ\n{self.Orbital_Occ}") 
         print(f"Orbital Character \n {self.Orbital_Character}")
         print(f"Orbital Localized Character \n {self.Orbital_Localized_Character}")



#test = Molecule("ag8exc.out","Ag")
#test.get_MOs()
#test2 = Molecule("ag19+1exc.out")
#test2.get_MOs()
#test2.get_exc_states("A1")
#test2.get_exc_decomp("A1")
#
#print(test2.Excited_State_Decomp["A1"]["1"])
#test2.print_by_transition_dipole_moment("A1", 2, 2.2,"test.tex")





