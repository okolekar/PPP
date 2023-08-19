# Mesh List Generator

def gen_Meshlist(N,M,s,length):
    #N = No of elements in the liquid length
    nodes  = []
    values = []
    generated_list = {}
    for i in range(N+1):
        nodes.append(i)
        nodes.append(i+N+M+1)
        values.append(i*s/(N))
        values.append(i*s/N+0.5*s/N)
    nodes.pop(-1)
    values.pop(-1)
    for i in range(N):
        generated_list[i] = {}
        generated_list[i][nodes[2*i]] = values[2*i]
        generated_list[i][nodes[2*i+1]] = values[2*i+1]
        generated_list[i][nodes[2*i+2]] = values[2*i+2]
    yield generated_list

#===========================================================================================================================#
                                                #Solid Mesh list generator                                                  
#===========================================================================================================================#
    nodes  = []
    values = []
    generated_list = {}
    index = 0
    for i in range(N,N+M+1):
        nodes.append(i)
        nodes.append(i+N+M+1)
        values.append(s + index*(length-s)/(M))
        values.append(s + index*(length-s)/M+0.5*(length-s)/M)
        index += 1
    nodes.pop(-1)
    values.pop(-1)
    index = 0
    for i in range(N,M+N):
        generated_list[i] = {}
        generated_list[i][nodes[2*index]] = values[2*index]
        generated_list[i][nodes[2*index+1]] = values[2*index+1]
        generated_list[i][nodes[2*index+2]] = values[2*index+2]
        index += 1
    yield generated_list
