#**************************************************************************
# Molecule class Reads ADF output, organizes MOs, organizes excited state 
# info by symmetry 
# Alva Dillon, Hongyi Wu, Hanwen Liu Spring 2022
#**************************************************************************

import traceback
import sys

class Molecule:
    filename = None
    Orbital_Occ = {} # {orbital label: occupation}
    Orbital_Energy = {} # {orbital label: MO Energy}
    Orbital_Localized_Character_Norm = {} # {orbital label: [percent metal , percent molecule]}
    Orbital_Localized_Character = {} # {orbital label: [percent metal , percent molecule]}
    Orbital_Character = {} # {orbital label:[s character, p character, d character]}
    Orbital_Character_Norm = {} # {orbital label:[s character, p character, d character]}
    metal_type= None
    Excited_States= {} # {Symm:{exc_num: [energy, osc_str, tdm_x, tdm_y, tdm_z]} } 
    Excited_State_Decomp = {} # {Symm: {Exc_num:{transition label: [weight, tdm_x, tdm_y, tdm_z]} } }
    CT_Osc_Product ={} #{Symm: {Exc_num: [metal_CT_char*osc_str, mol_CT_char*osc_str]}}
    CT_Excited_State={} #{Symm: {Exc_num: [metal_CT_char, mol_CT_char]}}

    def __init__(self, ADF_Outname, metal="Ag"):
        self.Orbital_Occ = {} # {orbital label: occupation}
        self.Orbital_Energy = {} # {orbital label: MO Energy}
        self.Orbital_Localized_Character_Norm = {} # {orbital label: [percent metal , percent molecule]}
        self.Orbital_Localized_Character = {} # {orbital label: [percent metal , percent molecule]}
        self.Orbital_Character = {} # {orbital label:[s character, p character, d character]}
        self.Orbital_Character_Norm = {} # {orbital label:[s character, p character, d character]}
        self.metal_type= None
        self.Excited_States= {} # {Symm:{exc_num: [energy, osc_str, tdm_x, tdm_y, tdm_z]} } 
        self.Excited_State_Decomp = {} # {Symm: {Exc_num:{transition label: [weight, tdm_x, tdm_y, tdm_z]} } }
        self.CT_Osc_Product ={} #{Symm: {Exc_num: [metal_CT_char*osc_str, mol_CT_char*osc_str]}}
        self.CT_Excited_State={} #{Symm: {Exc_num: [metal_CT_char, mol_CT_char]}}

        self.filename = str(ADF_Outname)
        self.metal_type = metal 
        self.get_MOs()
        self.orb_character_norm()
          
    def print_abs_spec(self,symmetry,energies, intensities):
        '''
        requires symmetry, array of energies, and array of intensities
        Function writes absorption intensities to file
        '''
        file1 = open(f'{self.filename}_{symmetry}_intensities.txt','w+')
        file1.write(f'Energy \t total \t metal \t mol \n')
        for point in range(len(intensities[0])):
            file1.write(f"{energies[point]} \t {intensities[0][point]} \t {intensities[1][point]} \t {intensities[2][point]} \n")
        file1.close()
 
    def pseudo_jablonski(self,plot_name):
        '''
        Function takes single input file and makes a jablonski like diagram 
        with linetypes that indicate the charge transfer character of an excited state
        '''
        import matplotlib.pyplot as plt
        Exc_state_prop = self.Excited_States
        CT_Char = self.CT_Excited_State
        energies1=[]
        char1 = []
        energies2=[]
        char2 = []
        energies3=[]
        char3 = []
        for sym in self.Excited_States:
            for num in self.Excited_States[sym]:
                print(f"state energies and ct character {self.Excited_States[sym][num]}  {self.CT_Excited_State[sym][num][0]}")
                if(self.CT_Excited_State[sym][num][0] < 0.25):
                    energies1.append(self.Excited_States[sym][num][0])
                    char1.append(self.CT_Excited_State[sym][num][0])
                elif(self.CT_Excited_State[sym][num][0] <= 0.5 and self.CT_Excited_State[sym][num][0] >= 0.25):
                    energies2.append(self.Excited_States[sym][num][0])
                    char2.append(self.CT_Excited_State[sym][num][0])
                elif(self.CT_Excited_State[sym][num][0] > 0.5):
                    energies3.append(self.Excited_States[sym][num][0])
                    char3.append(self.CT_Excited_State[sym][num][0])

        plt.hlines(energies1,xmin=0.25, xmax = 0.75, color='grey')
        plt.hlines(energies2,xmin=0.25, xmax = 0.75, color='blue')
        plt.hlines(energies3,xmin=0.25, xmax = 0.75, color='red')
        plt.xlim(0,1)
        plt.show()

    def multi_pseudo_jablonski(self, multiExcitedStates, multiCTStates,plt_title,labels):
        '''
        Function takes multiple dictionaries and makes a jablonski like diagram 
        with linetypes that indicate the charge transfer character of an excited state
        '''
        import matplotlib.pyplot as plt
        xs=0.25
        xe=0.75
        fig, ax = plt.subplots(1,1,figsize=(6,14))
        xe_list=[]
        max_e=0
        for group in range(len(multiExcitedStates)):
            energies1=[]
            char1 = []
            energies2=[]
            char2 = []
            energies3=[]
            char3 = []
            for sym in multiExcitedStates[group]:
                for num in multiExcitedStates[group][sym]:
                    #print(f"state energies and ct character {self.Excited_States[sym][num]}  {self.CT_Excited_State[sym][num][0]}")
                    if(multiExcitedStates[group][sym][num][0] > max_e):
                        max_e = multiExcitedStates[group][sym][num][0]
                    if(multiCTStates[group][sym][num][0] < 0.25):
                        energies1.append(multiExcitedStates[group][sym][num][0])
                        char1.append(multiCTStates[group][sym][num][0])
                    elif(multiCTStates[group][sym][num][0] <= 0.5 and multiCTStates[group][sym][num][0] >= 0.25):
                        energies2.append(multiExcitedStates[group][sym][num][0])
                        char2.append(multiCTStates[group][sym][num][0])
                    elif(multiCTStates[group][sym][num][0] > 0.5):
                        energies3.append(multiExcitedStates[group][sym][num][0])
                        char3.append(multiCTStates[group][sym][num][0])
             
            ax.hlines(energies1,xmin=xs, xmax=xe, color='grey',label='Ag12')
            ax.hlines(energies2,xmin=xs, xmax=xe, color='blue')
            ax.hlines(energies3,xmin=xs, xmax=xe, color='red')
            xe_list.append(xe-.25)
            xe+=0.75
            xs+=0.75
        ax.set_xticks(xe_list)
        ax.set_xticklabels(labels, fontsize=18)
        plt.title(f"${plt_title}$",fontsize=18)
        plt.ylim(0,max_e)
        plt.xlim(0,xe)
        plt.show() 


    def make_lorentzian_plot(self, max_energy, plot_name,sticks=False,createFile=False, noCT=False):
        '''
        Function to make absorption spectra and the absorption spectra scaled by charge transfer for each excited state
        input required: symmetry of system, maximum energy of absorption spectra, name of .png file
        '''

        import math
        import matplotlib.pyplot as plt
        #lorentezian parameters
        gamma = 0.1088
        energy_step= 0.02
        steps = int(max_energy/energy_step) + 1
        intensity = []
        intensity2 = []
        intensity3 = []
        intensity4 = []
        en_ax = []
        for i in range(0,steps):
            en_ax.append(float(i*energy_step))
            intensity.append(0.0)
            intensity2.append(0.0)
            intensity3.append(0.0)
            intensity4.append(0.0)

        fig = plt.figure(figsize=(6.54,6.54),dpi=900)
        plt.rcParams.update({'font.size': 20, 'font.weight':'medium'})
        ax = fig.add_subplot(111)
        stick=[]
        stick_h=[]
        #CT_Excited_State={} #{Symm: {Exc_num: [metal_CT_char, mol_CT_char]}}
        #Excited_States= {} # {Symm:{exc_num: [energy, osc_str, tdm_x, tdm_y, tdm_z]} } 
        print(self.Excited_States)
        for sym in self.Excited_States:
            for exc_num in self.Excited_States[sym]:
                if float(self.Excited_States[sym][exc_num][0]) < max_energy:
                    stick.append(float(self.Excited_States[sym][exc_num][0]))
                    stick_h.append(float(self.Excited_States[sym][exc_num][1]))
                    for i in range(len(intensity)):
                        energy_point = en_ax[i]
                        phi = gamma /(((energy_point - float(self.Excited_States[sym][exc_num][0]))**2 + gamma**2)* math.pi) 
                        if "E" in sym:
                            intensity[i] += phi * float(self.Excited_States[sym][exc_num][1]) * 2     
                            intensity2[i] += phi * float(self.CT_Excited_State[sym][exc_num][0]) * float(self.Excited_States[sym][exc_num][1]) * 2
                            intensity3[i] += phi * float(self.CT_Excited_State[sym][exc_num][1]) * float(self.Excited_States[sym][exc_num][1]) * 2
                        if "T" in sym:
                            intensity[i] += phi * float(self.Excited_States[sym][exc_num][1]) * 3   
                            intensity2[i] += phi * float(self.CT_Excited_State[sym][exc_num][0]) * float(self.Excited_States[sym][exc_num][1]) * 3
                            intensity3[i] += phi * float(self.CT_Excited_State[sym][exc_num][1]) * float(self.Excited_States[sym][exc_num][1]) * 3
                            #intensity4[i] += phi * float(self.CT_Excited_State[symmetry][exc_num][1]) * float(self.Excited_States[symmetry][exc_num][1]) * 2
                        if "A" in sym or "S" in sym or "B" in sym:
                            intensity[i] += phi * float(self.Excited_States[sym][exc_num][1]) * 1     
                            intensity2[i] += phi * (float(self.CT_Excited_State[sym][exc_num][0])) * float(self.Excited_States[sym][exc_num][1]) * 1
                            intensity3[i] += phi * (float(self.CT_Excited_State[sym][exc_num][1])) * float(self.Excited_States[sym][exc_num][1]) * 1
                            #intensity4[i] += phi * abs(float(self.CT_Excited_State[symmetry][exc_num][1])) * float(self.Excited_States[symmetry][exc_num][1]) * 1
            #print(f"intensity {intensity}")
        if (sticks == True):
            ax2 = ax.twinx()
            ax2.vlines(x = stick, ymin = 0, ymax = stick_h,linewidth=0.5,color="black")  
            ax2.set_ylim(ymin=0,ymax=max(stick_h)*2)
            ax2.set_ylabel("Osc Str.")

        ax.plot(en_ax, intensity,'-b',label="total abs")
        ax.set_xlabel("Energy")

        if (noCT == False):
            ax.plot(en_ax, intensity2,'-r',label="metal ct")
            ax.plot(en_ax, intensity3,'-g', label="mol ct")

        ax.set_ylabel("Intensity")
        ax.set_ylim(bottom=0)
        ax.legend(prop={'size': 10})
        plt.tight_layout()
        #plt.show()
        plt.savefig(f"{plot_name}.png")

        if (createFile == True):
            self.print_abs_spec('sigma',en_ax, [intensity,intensity2,intensity3])





    def orb_character_norm(self):
        '''
        normalizes s,p,d character of molecular orbitals
        '''
        for orb in self.Orbital_Character:
            self.Orbital_Character_Norm[orb] = [float(i)/sum(self.Orbital_Character[orb]) for i in self.Orbital_Character[orb]]





    def make_ct_table(self, symmetry, component_cutoff, outfile):
        '''
        prints latex table of excited state properties:
        state number, state energy, oscillator strength, charge transfer character
            , transition , weight, orbital character
        component_cutoff is for charge transfer oscillator strength product
        '''
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
 \\begin{longtable}{c|c|c|c}
   \\caption{Excited States and Electron Configurations Table}
   \\hline
     State & Energy/transition & Osc. Str./weight & charge transfer\\tiny({mol->metal,metal->mol})\\normalsize/Orb Char (s,p,d)\\\\
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

        states = self.Excited_States
        state_decomp = self.Excited_State_Decomp
        ct_char = self.CT_Excited_State
        ct_osc_prod = self.CT_Osc_Product
        orbs = self.Orbital_Character_Norm
        output = open(outfile,'a')
        #print(self.CT_Osc_Product)
        #print(self.CT_Excited_State)
        for exc in states[symmetry]:
            if abs(float(ct_osc_prod[symmetry][exc][1])) > component_cutoff or abs(float(ct_osc_prod[symmetry][exc][0])) > component_cutoff:
                output.write(f"{exc} &\t {states[symmetry][exc][0]:.2f} &\t {states[symmetry][exc][1]:.2f} &\t {ct_char[symmetry][exc][0]:.2f},{ct_char[symmetry][exc][1]:.2f} \\\\ \n  ")
                for tran in state_decomp[symmetry][exc]:
                    t_split = tran.split("$\\rightarrow$")
                    occ = t_split[0].replace(" ","")
                    unocc = t_split[1].replace(" ","")
                    #print(orbs) 
                    #print(f"occupied orb {occ} unoccupied orb {unocc}")
                    try:
                        if float(state_decomp[symmetry][exc][tran][0]) > 0.01: 
                            output.write(f"\t& {occ}({self.Orbital_Localized_Character_Norm[occ][0]:.2f})$\\rightarrow${unocc}({self.Orbital_Localized_Character_Norm[unocc][0]:.2f}) &\t {state_decomp[symmetry][exc][tran][0]:.4f} &\t {float(orbs[occ][0]):.2f}, {float(orbs[occ][1]):.2f}, {float(orbs[occ][2]):.2f} $\\rightarrow$ {float(orbs[unocc][0]):.2f}, {float(orbs[unocc][1]):.2f}, {float(orbs[unocc][2]):.2f} \\\\ \n ")                     
                    except:
                        output.write(f"\t& {tran} &\t {state_decomp[symmetry][exc][tran][0]:.4f} &\t ? $\\rightarrow$ ?  \\\\ \n ") 
                output.write("\\hline")
        output.close()
        self.make_latex_table_end(outfile)




    def calc_ct_character(self, symm):
        ''' 
        calculates charge transfer character of excited state from orbitals localized on metal or molecule
        '''
        ct_osc={}
        ct_exc={}
        self.CT_Excited_State[symm]={}
        self.Orbital_Localized_Character_Norm = {}
        self.CT_Osc_Product[symm]={}     
        count=1000
        for state in self.Excited_State_Decomp[symm]:
            self.CT_Excited_State[symm][state]={}     
            self.CT_Osc_Product[symm][state]={}
            #occ_vir_diff = [0,0]#[part metal , part molecule]
            metal_ct = 0.0
            mol_ct = 0.0
            #print("")
            for transition in self.Excited_State_Decomp[symm][state]:
                weight = self.Excited_State_Decomp[symm][state][transition][0]
                transition = transition.split("$\\rightarrow$")
                occ_orb = transition[0].replace(" ", "")
                vir_orb = transition[1].replace(" ", "")
                #print(f"occ_orb:{occ_orb}")
                #print(f"vir_orb:{vir_orb}")
                self.Orbital_Localized_Character_Norm[occ_orb] = []
                self.Orbital_Localized_Character_Norm[vir_orb] = [] 
                #occ_char = self.Orbital_Localized_Character[occ_orb]
                if occ_orb in self.Orbital_Localized_Character:
                    occ_norm = [float(i)/sum(self.Orbital_Localized_Character[occ_orb]) for i in self.Orbital_Localized_Character[occ_orb]]
                    self.Orbital_Localized_Character_Norm[occ_orb]=occ_norm
                    #print(f"occupied_norm {occ_norm}")
                    count +=1 
                if count == 1000:
                    exit()
                #vir_char = self.Orbital_Localized_Character[vir_orb]
                try:
                    vir_norm = [float(i)/sum(self.Orbital_Localized_Character[vir_orb]) for i in self.Orbital_Localized_Character[vir_orb]]
                    self.Orbital_Localized_Character_Norm[vir_orb]=vir_norm
                    #print(f"virtual_norm {vir_norm}")
                except:
                    #print(f"orbital {vir_orb} is not in your list of orbitals") 
                    pass
                    #vir_norm=[0,0]
                #print(f"occ_norm {occ_norm} vir_norm {vir_norm}")
                if vir_norm[0] - occ_norm[0] > 0 :
                    metal_ct += (vir_norm[0] - occ_norm[0]) * weight
                else:
                    mol_ct += (vir_norm[1] - occ_norm[1]) * weight
                #print(f"{metal_ct} { mol_ct}  inloop")
            try:
                #print(f"{metal_ct} {mol_ct} outloop")
                occ_vir_diff=[metal_ct, mol_ct] #occ_vir_diff = [part metal , part molecule]
                #print(f"ct {state} : {occ_vir_diff}") 
                #print(f"ct*osc_str : {occ_vir_diff[0]*float(self.Excited_States[symm][state][1])} ")                 
                self.CT_Excited_State[symm][state]=occ_vir_diff
                self.CT_Osc_Product[symm][state]=[occ_vir_diff[0]*float(self.Excited_States[symm][state][1]), occ_vir_diff[1]*float(self.Excited_States[symm][state][1]) ]
            except:
                pass 
                


    
    def make_xyz(self, outfile):
        '''
        makes file in .xyz format from ADF output file 
        and was created to view the progression of a 
        geometry optimization. 
        python ADF_output_file ouput_xyz_filename
        '''    
        filename=self.filename
        fi = open(filename,'r')
        outfi=open(outfile+".xyz","w")
        li = fi.readline()
        numAtoms=0
        while li:
            li= fi.readline()
            if "ATOMS" in li:
                #print(li)
                li=fi.readline()
                li=fi.readline()
                li=fi.readline()
                while "FRAGMENTS" not in li:
                    li=fi.readline()
                    if len(li.split()) >2:
                        numAtoms+=1
                    else:
                        break
                    #print(li)
                #print("number of Atoms %i "%(numAtoms)) 
            elif "Coordinates in Geometry Cycle" in li and numAtoms > 0: 
                outfi=open(outfile+".xyz","a")
                outfi.write(str(numAtoms))
                outfi.write("\n \n")
                coords=[]
                li=fi.readline()
                for i in range(numAtoms):
                    li =fi.readline()
                    line=li.split()
                    tempStr = line[0].split(".")
                    line[0] = tempStr[1]
                    for j in range(len(line)):
                        outfi.write(line[j] + "   ")
                    outfi.write("\n")
                    coords.append(line)
                    #print(line)
                #outfi.write(str(coords))
                #outfi.write("")
                
                #print(coords)
        outfi.close()
        fi.close()

    

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
 \\begin{longtable}{c|c|c|c|c|c|c}
   \\caption{Excited States and Electron Configurations Table}
   \\hline
     State & Energy & Osc. Str. & Transition & Weight & Transition Dipole Moment(x,y,z) & Orb Char (s,p,d)\\\\
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
        occ_orb = transition[0].replace(" ","")
        vir_orb = transition[2].replace(" ","")
        transition = f"{str(occ_orb)} $\\rightarrow$ {str(vir_orb)}"
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
                output.write(f"{exc} &\t {states[symmetry][exc][0]:.2f} &\t {states[symmetry][exc][1]:.2f} &\t &\t &\t {states[symmetry][exc][2]:.2f}, {states[symmetry][exc][3]:.2f}, {states[symmetry][exc][4]:.2f} &\t   \\\\ \n  ")
                for tran in state_decomp[symmetry][exc]:
                    t_split = tran.split("$\\rightarrow$")
                    occ = t_split[0].replace(" ","")
                    unocc = t_split[1].replace(" ","")
                    #print(orbs) 
                    #print(f"occupied orb {occ} unoccupied orb {unocc}")
                    try:
                        output.write(f"\t& \t& \t& {tran} &\t {state_decomp[symmetry][exc][tran][0]:.4f} &\t {state_decomp[symmetry][exc][tran][1]:.3f}, {state_decomp[symmetry][exc][tran][2]:.3f}, {state_decomp[symmetry][exc][tran][3]:.3f}  &\t {float(orbs[occ][0]):.2f}, {float(orbs[occ][1]):.2f}, {float(orbs[occ][2]):.2f} $\\rightarrow$ {float(orbs[unocc][0]):.2f}, {float(orbs[unocc][1]):.2f}, {float(orbs[unocc][2]):.2f} \\\\ \n ") 
                    
                    except:
                        output.write(f"\t& \t& \t& {tran} &\t {state_decomp[symmetry][exc][tran][0]:.4f} &\t {state_decomp[symmetry][exc][tran][1]:.3f}, {state_decomp[symmetry][exc][tran][2]:.3f}, {state_decomp[symmetry][exc][tran][3]:.3f}  &\t ? $\\rightarrow$ ?  \\\\ \n ") 
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
        #print(f"Excited States \n {self.Excited_States}")


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
        #print(f"Excited state decompositions \n {self.Excited_State_Decomp}")




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
                         #print(char_type)
                         metal_char += loc_type[0]
                         mol_char += loc_type[1]
                         s_char += char_type[0]
                         p_char += char_type[1]
                         d_char += char_type[2]
                         #print(f"{s_char}, {p_char}, {d_char}") 
                         orbital = str(sp_line[2])+str(sp_line[3])
                         orbital = orbital.split(":")
                         orbital = orbital[0].replace(" ","")
                         if "E1.u" in orbital:
                             orbital=orbital.replace("E1.u","pi.u")
                             #print(f"changed orbital: {orbital}")
                         if "E1.g" in orbital:
                             orbital=orbital.replace("E1.g","pi.g")
                             #print(f"changed orbital: {orbital}")
                         if "E1" in orbital:
                             orbital=orbital.replace("E1","e")
                             #print(f"changed orbital: {orbital}")
                         if "A1.g" in orbital:
                             orbital=orbital.replace("A1.g","s+.g")
                             #print(f"changed orbital: {orbital}")
                         if "A2.u" in orbital:
                             orbital=orbital.replace("A2.u","s+.u")
                             #print(f"changed orbital: {orbital}")
                         
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
         #print(f"Orbital Energies \n{self.Orbital_Energy}")
         #print(f"Orbital Occ\n{self.Orbital_Occ}") 
         #print(f"Orbital Character \n {self.Orbital_Character}")
         #print(f"Orbital Localized Character \n {self.Orbital_Localized_Character}")



#test = Molecule("ag8exc.out","Ag")
#test.get_MOs()
#test.get_exc_states("S+.u")
#test2 = Molecule("pyr_edge_exc.out")
#test2.get_MOs()
#test2.get_exc_states("A")
#test2.get_exc_decomp("A")
#test2.print_by_transition_dipole_moment("A1", 2, 2.2,"test.tex")
#test.calc_ct_character("S+.u")
#test2.make_ct_table("A", 0.001, "test.tex")
#test.make_lorentzian_plot("S+.u", 6, "mol_test")
#
#print(test2.Excited_State_Decomp["A1"]["1"])





class PlotFiles:
    def multi_spec(self,filenames,labels,plot_name,noCT=False):
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(6.54,6.54),dpi=600)
        plt.rcParams.update({'font.size': 20, 'font.weight':'medium'})
        ax = fig.add_subplot(111)
        for i in range(len(filenames)):
            fi = open(filenames[i],'r')
            intensity=[]
            intensity2=[]
            intensity3 = []
            en_ax = []
            line = fi.readline()
            line = fi.readline()
            while line:
                sp_line = line.split()
                en_ax.append(float(sp_line[0]))
                intensity.append(float(sp_line[1]))
                intensity2.append(float(sp_line[2]))
                intensity3.append(float(sp_line[3]))
                line = fi.readline()
            fi.close()

            ax.plot(en_ax, intensity,'-',label=f"{labels[i]}")
            if(noCT == False):
                ax.plot(en_ax, intensity2,'-',label=f"{labels[i]}"+"$_{metal ct}$")
                ax.plot(en_ax, intensity3,'-', label=f"{labels[i]}"+"$_{mol ct}$")

        ax.set_xlabel("Energy")
        ax.set_ylabel("Intensity")
        ax.set_ylim(bottom=0)
        ax.legend(prop={'size': 10})
        plt.tight_layout()
        #plt.show()
        plt.savefig(f"{plot_name}.png")





    def single_spec(self,filename,plot_name):
        '''
        Function takes input from file to make single absorption spectra
        '''

        import matplotlib.pyplot as plt
        fi = open(filename,'r')
        intensity=[]
        intensity2=[]
        intensity3 = []
        en_ax = []
        line = fi.readline()
        line = fi.readline()
        while line:
            sp_line = line.split()
            en_ax.append(float(sp_line[0]))
            intensity.append(float(sp_line[1]))
            intensity2.append(float(sp_line[2]))
            intensity3.append(float(sp_line[3]))
            line = fi.readline()
        fi.close()

        fig = plt.figure(figsize=(6.54,6.54),dpi=600)
        plt.rcParams.update({'font.size': 20, 'font.weight':'medium'})
        ax = fig.add_subplot(111)
        ax.plot(en_ax, intensity,'-b',label="total abs")
        ax.set_xlabel("Energy")
       # ax.set_ylabel("total")
        ax.plot(en_ax, intensity2,'-r',label="metal ct")
       # ax.set_ylabel("metal")
        ax.plot(en_ax, intensity3,'-g', label="mol ct")
        #ax.plot(en_ax, intensity4,'-m', label="abs ct")
        ax.set_ylabel("Intensity")
        ax.set_ylim(bottom=0)
        ax.legend(prop={'size': 10})
        plt.tight_layout()
        #plt.show()
        plt.savefig(f"{plot_name}.png")
            











## TESTS
class TestClass:
    def run_pyrEdge_test():
        pytest = Molecule("pyr_edge_exc.out","Ag")
        #pytest.get_MOs()
        print(pytest.Orbital_Energy)
        pytest.get_exc_states("A")
        print("finished get_exc_states")
        pytest.get_exc_decomp("A")
        print("finished get_exc_decomp")
        pytest.calc_ct_character("A")
        print("finished get_calc_ct_character")
        pytest.make_lorentzian_plot(6, "pyr_edge_test",createFile=True)
        del pytest




    def run_Ag8_test():
        agtest = Molecule("ag8exc.out","Ag")
        #agtest.get_MOs()
        print(agtest.Orbital_Energy)
        agtest.get_exc_states("S+.u")
        agtest.get_exc_states("Pi+.u")
        print("finished get_exc_states")
        agtest.get_exc_decomp("S+.u")
        agtest.get_exc_decomp("Pi+.u")
        print("finished get_exc_decomp")
        agtest.calc_ct_character("S+.u")
        agtest.calc_ct_character("Pi+.u")
        print("finished get_calc_ct_character")
        agtest.make_lorentzian_plot(6, "Ag8_test", noCT=True, createFile=True)
        del agtest

    def run_Ag19_test():
        ag19test = Molecule("Ag19_minus_exc.out","Ag")
        #ag19test.get_MOs()
        print(ag19test.Orbital_Energy)
        ag19test.get_exc_states("A1")
        ag19test.get_exc_states("E")
        print("finished get_exc_states")
        ag19test.get_exc_decomp("A1")
        ag19test.get_exc_decomp("E")
        print("finished get_exc_decomp")
        ag19test.calc_ct_character("A1")
        ag19test.calc_ct_character("E")
        print("finished get_calc_ct_character")
        ag19test.make_lorentzian_plot(6, "Ag19_minus_test",noCT=True,sticks=True, createFile=True)
        del ag19test

    def run_spec_test():
        test = PlotFiles()
        test.single_spec('pyr_edge_exc.out_sigma_intensities.txt',"test_single")
        print("single_spec run successfully! ")
        test.multi_spec(['pyr_edge_exc.out_sigma_intensities.txt','ag8exc.out_sigma_intensities.txt','Ag19_minus_exc.out_sigma_intensities.txt'],["pyr edge",'Ag 8', 'Ag 19'],"test_multi",noCT=True)
        print("multi_spec run successfully! ")
        
 

def run_test():
    try:
        TestClass.run_pyrEdge_test()
        print("pyr_edge_exc.out test successful !!")
    except Exception as e:
        print("Error in pyr_edge_exc.out test")
        print(e)
        traceback.print_exc()
    try:
        TestClass.run_Ag8_test()
        print("Ag8 test successful !!")
    except Exception as e:
        print("Error in Ag8 test")
        print(e)
        traceback.print_exc()
    
    try:
        TestClass.run_Ag19_test()
        print("Ag19 test successful !!")
    except Exception as e:
        print("Error in Ag19 test")
        print(e)
        traceback.print_exc()

    try:
        TestClass.run_spec_test()
        print("Spec test done")
    except Exception as e:
        print("Error in run_spec_test")
        print(e)
        traceback.print_exc()









