from ct_analysis import Molecule

multict=[]
multiex=[]
t2 = Molecule("pyr_edge_exc.out")
t2.get_MOs()
t2.get_exc_states("A")
t2.get_exc_decomp("A")
t2.calc_ct_character("A")
multiex.append(t2.Excited_States)
multict.append(t2.CT_Excited_State)
t2 = Molecule("Ag19_minus_exc.out")
t2.get_MOs()
t2.get_exc_states("A1")
t2.get_exc_decomp("A1")
t2.calc_ct_character("A1")
t2.get_exc_states("E")
t2.get_exc_decomp("E")
t2.calc_ct_character("E")
multiex.append(t2.Excited_States)
multict.append(t2.CT_Excited_State)
print(range(len(multiex)))
t2.multi_pseudo_jablonski(multiex, multict,"Plot Title",["pyr edge", "Ag19$^-$"])
