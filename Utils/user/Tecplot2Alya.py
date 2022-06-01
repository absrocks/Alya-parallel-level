'''

This program converts Tecplot's ascii point format to Alya's format.

TECPLOT: *.dat --> ALYA: *.dom.geo, *.fix.dat, *dims.dat

Author: Edgar Olivares

'''

import math
import operator

# Counting the number of physical entities in the problem
def sum_phy(namef):
    sum_p=0
    g=open(namef+'.dat','r')
    for l in g:
        fi = l.split()
        if fi[0].startswith('ZONE') and fi[1].startswith('T='):
            sum_p = sum_p+1
    g.close()
    return sum_p

# Write the name of each physical entity and turn on booleans "r" and "bound"
def physics(sum_physics,r,bound,phy_ent,fi): 
    if fi[1].startswith('T="Surf'):
        r = 'F'; bound = 'F'
        sum_physics = sum_physics + 1
        name = (l.split(':')[-1]).split('"')[0]
        phy_ent[sum_physics] = name
        zone = sum_physics
    if fi[1].startswith('T="Vol'):
        r = 'F'; bound = 'T'
        name = (l.split(':')[0]).split('"')[1]
        phy_ent[sum_physics+total] = name
        zone = 0
    return zone,phy_ent,r,bound,sum_physics

# Looking up boundaries nodes and saving them in a dictionary
def boundaries(sum_bound, sum_tot, bound_num, zone, fi, bound_dict_z):
    sum_bound = sum_bound + 1
    sum_tot = sum_tot + 1
    dic1 = {}
    bound_num[sum_bound] = zone
    for i in range(1,total+1):
        if zone == i:
            for j in range(0,len(fi)):
                nodeb = int(float(fi[j]))
                if not nodeb in dic1:
                    dic1[nodeb] = 1
                    bz = (sum_bound,zone)
                    if bz in bound_dict_z:
                        bound_dict_z[bz].append(nodeb)
                    else:
                        bound_dict_z[bz] = [nodeb]
    return sum_bound, sum_tot, bound_num, bound_dict_z

# Looking up boundaries coordinates and saving them in a dictionary
def bcoord(fi,total,sum_z,bcoord_dict,zone):
    coord = (fi[0],fi[1],fi[2])
    for i in range(1,total+1):
        if zone == i:
            sum_z[i] = sum_z[i]+1
            bbb = (sum_z[i],i)
            bcoord_dict[bbb] = coord
    return sum_z, bcoord_dict

# Lookipng up elemnts and savim them in a dictionary
def elements(sum_elm,sum_tot,fi,elm_node_dict,nodes_per_elm):
    sum_elm = sum_elm + 1
    sum_tot = sum_tot + 1
    dic2={}
    for i in range(0,len(fi)):
        node = int(float(fi[i]))
        if not node in dic2:
            dic2[node] = 1
            if sum_elm in elm_node_dict:
                elm_node_dict[sum_elm].append(node)
            else:
                elm_node_dict[sum_elm] = [node]
    nodes_per_elm[sum_elm] = len(dic2)
    return sum_elm,sum_tot,elm_node_dict,nodes_per_elm
    

if __name__ == "__main__":
    # Reading the file
    namef = raw_input('What is the name of the file? -- whthout any extension --'+'\n'+'EX: file_ex.dat, write only file_ex'+'\n')

    total = sum_phy(namef)
    print '### --> Computed', total+1, 'different physics entities ###'

    # Initialize variables
    bound_dict_z={}; bcoord_dict={}; bound_num={}; nodal_p={}; nodal_inv={}; elm_node_dict={}; nodes_per_elm={}; phy_ent={} #Dictioaries
    sum_physics = 0; sum_bound = 0; sum_tot = 0; sum_nodes = 0; sum_elm = 0; #Summations
    r = 'N'; bound = 'N' #Booleans
    sum_z=[] # Ssumations in an array for each physical boundary
    for i in range(1,total+1):
        sum_z.append(0)

    # Start reading the input file
    f=open(namef +'.dat','r')
    for l in f:
        fi = l.split()
        if fi[0].startswith('ZONE'):        
            (zone,phy_ent,r,bound,sum_physics)=physics(sum_physics,r,bound,phy_ent,fi)
        if fi[0].startswith('DT=') and r == 'F':
            r = "T"
        if r == "T" and not fi[0].startswith('DT='):
            if math.modf(float(fi[0]))[0] == 0.0 and bound == "F":
                (sum_bound, sum_tot, bound_num, bound_dict_z) = boundaries(sum_bound, sum_tot, bound_num, zone, fi, bound_dict_z)   
            if math.modf(float(fi[0]))[0] != 0.0 and bound == "F":
                (sum_z, bcoord_dict) = bcoord(fi,total,sum_z,bcoord_dict,zone)
            if math.modf(float(fi[0]))[0] != 0.0 and bound == "T": # Looking up nodes and saving them in a dictionary
                sum_nodes = sum_nodes +1
                coord = (fi[0],fi[1],fi[2])
                nodal_p[sum_nodes] = coord
                nodal_inv[coord] = sum_nodes
            if math.modf(float(fi[0]))[0] == 0.0 and bound == "T":
                (sum_elm,sum_tot,elm_node_dict,nodes_per_elm) = elements(sum_elm,sum_tot,fi,elm_node_dict,nodes_per_elm) 
    f.close()

    print '### Physical domains'
    for k,v in phy_ent.iteritems():
        print k,' --> ', v
    print '###'

    # Writing outputs
    o3 = open(namef+'.dims.dat','w')
    o3.write('NODAL_PONINTS'+'\t'+str(sum_nodes)+'\n'+'ELEMENTS'+'\t'+str(sum_elm)+'\n'+'BOUNDARIES'+'\t'+ str(sum_bound))
    o3.close()
    
    o2 = open(namef+'.geo.dat', 'w')

    o2.write('TYPES'+'\n')
    sorted_nodes_per_elm = sorted(nodes_per_elm.items(), key=operator.itemgetter(0))
    for kv in sorted_nodes_per_elm:
        o2.write(str(kv[0])+' ')
        if kv[1] == 4:
            o2.write('30'+'\n')
        elif kv[1] == 5:
            o2.write('32'+'\n')
        elif kv[1] == 6:
            o2.write('34'+'\n')
        else:
            print '### --> UNKNOWN ELEMENT WITH', kv[1], 'NODES ###'
    o2.write('END_TYPES' + '\n')

    o2.write('ELEMENTS' + '\n')
    for k,v in elm_node_dict.iteritems():
        o2.write(str(k) + ' ')
        for i in range(len(v)):
            o2.write(str(v[i]) + ' ')
        o2.write('\n')
    o2.write('END_ELEMENTS' + '\n')

    o2.write('COORDINATES' + '\n')
    for k,v in nodal_p.iteritems():
        o2.write(str(k)+' '+str(v[0])+' '+str(v[1])+' '+str(v[2])+'\n')
    o2.write('END_COORDINATES'+'\n')

    o2.write('BOUNDARIES' + '\n')
    sorted_b=sorted(bound_dict_z.items(), key=operator.itemgetter(0))
    for kv in sorted_b:
        o2.write(str(kv[0][0])+' ')
        for i in range(len(kv[1])):
            coord1 = bcoord_dict[kv[1][i],kv[0][1]]
            node_abs = nodal_inv[coord1]
            o2.write(str(node_abs)+' ')
        o2.write('\n')
    o2.write('END_BOUNDARIES' + '\n')
    o2.write('SKEW_SYSTEMS'+'\n')
    o2.write('END_SKEW_SYSTEMS'+'\n')
    o2.close()

    o = open(namef+'.fix.bou','w')
    o.write('ON_BOUNDAIRES' + '\n')
    for k,v in bound_num.iteritems():
        o.write(str(k)+' '+str(v)+'\n')
    o.write('END_ON_BOUNDARIES')
    o.close()
