'''
  [Cancer growth model]
 The code using the Object functions and structs to define the ID
 and mutation type of daughter cell, and the overall code
 is shorter and several times faster, in addition, the test
 code could detect the whole cancer cell VAF.
'''
import numpy as np
import random
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors   # define the colour by myself
import pandas as pd     # Applied to cell pushing
import time     # calculate the spend time of code
from collections import Counter     # count the frequency of mutation



start_time = time.time()

### Basic parameter setting
divide_rate = 1
background_gene = 1   # The background gene that cancer originally carried
mutation_rate = 1
death_rate = 0.01
die_divide_times = 50   # The max times of cell divides,  arrived this point the cell will die
Poisson_lambda = 10
Max_ROW = 120
Max_COL = 120
Max_generation = 45

mesh_shape = (Max_ROW, Max_COL)     # growth grid
cancer_matrix = np.zeros(mesh_shape)    # create grid space of the cancer growth


def cancer_mesh(row,col):    # define cancer growth mesh size
    all_mesh = [(r,c)for r in np.arange(0,row) for c in np.arange(0, col)]
    central_pos = {cp:neighbor_mesh(cp,row,col) for cp in all_mesh}  # define the center position
    return central_pos

def neighbor_mesh(cp,row,col):  # define the neighbor of the center point.
    r, c = cp     # cp: central position   r: cp.row  c: cp.column
    neighbor_pos = [(r+i, c+j)
        for i in [-1, 0, 1]
        for j in [-1, 0, 1]
        if 0 <= r + i < row     # define the point in the grid
        if 0 <= c + j < col
        if not (j == 0 and i == 0)]
    return neighbor_pos

# Function simplification
binomial = np.random.binomial
shuffle = np.random.shuffle
randint = np.random.randint
random_choice = random.choice
poisson = np.random.poisson

def divide_r(): # divide rate
    divide = binomial(1, divide_rate)
    return divide
def mutation_r():   # mutation rate
    mutation = binomial(1, mutation_rate)
    return mutation
def death_r():  #death rate
    death = binomial(1, death_rate)
    return death
def Poisson():   # mutation numbers from Poisson distribution, where the lambda=10
    add_mutation =poisson(Poisson_lambda)     # cancer_num: the total number of cancer cells in the times
    # change the size to 1

    return add_mutation


class Cancercell ():

    def __init__(self, cp, neighbor_mesh):
        """
               Initialize Cancer Cell
               :param pos: position of cancer cell; tuple
               :param dictionary_of_neighbor_mesh: used to tell the cell the positions of its neighbors

        """
        self.cp = cp  # cp means central position
        self.die_divide_times = die_divide_times  # define the times of cell divide
        self.neighbor_pos= neighbor_mesh[self.cp]  # define the neighbor central point
        self.ID = 1  # the cell ID type, 1 main back_cell
        self.mutation = list(range(background_gene))  # mutation data add the back gene inform
        self.mu_times = 0  # the initial divide times of mutation divide
        self.pu_times = Max_generation
        self.all_mutation_ID = []
        global times
        times = 0

    def empty_neighbor(self, agent):
        """
                Search for empty positions in More neighborhood. If there is more than one free position,
                randomly select one and return it
                :param agent: dictionary of agents, key=position, value = cell.ID; dict
                :return: Randomly selected empty position, or None if no empty positions
        """
        empty_neighbor_list = [cp for cp in self.neighbor_pos if cp not in agent ]
        if empty_neighbor_list:
            empty_pos = random_choice(empty_neighbor_list)
            return empty_pos
        else:
            return None

    def act(self, agent, neighbor_mesh):
        """
                Cell carries out its actions, which are division and death. Cell will divide if it is lucky and
                there is an empty position in its neighborhood. Cell dies either spontaneously or if it exceeds its
                maximum number of divisions.

                :param agent: dictionary of agents, key=position, value = cell.ID; dict
                :return: None
        """
        ## Call definition function
        division = divide_r()
        variant = mutation_r()
        mu = Max_generation
        nu = 0
        self.mu_times += 1
        self.pu_times += 1
        empty_t = self.mu_times
        push_t = self.pu_times
        if division == 1:

            empty_pos = self.empty_neighbor(agent)



            poisson = Poisson()

            if empty_pos is not None:

                daughter_cell = Cancercell(empty_pos, neighbor_mesh)
                  # Creat new daughter cell and it to the cell dictionary
                daughter_cell.ID = Cancercell(empty_pos, neighbor_mesh).ID  # define the daughter ID same as parent
                  # define the daughter mutation type same as parent

                for i in range(max(self.mutation) + 1, max(self.mutation) + poisson + 1):
                    self.mutation.append(i)

                    agent[empty_pos] = daughter_cell
                    daughter_cell.mutation = self.mutation+list(range(max(daughter_cell.mutation) + 1, max(daughter_cell.mutation) + poisson + 1))
                    daughter_cell.ID = self.ID + random_choice(range(10))


                # run the max generation create new mutation
                # when the variant value 1 means cells create new mutation


                  # give new ID to new cells ()yi'chu
                # self.all_mutation_ID = Cancercell(empty_pos, neighbor_mesh).mutation + agent[empty_pos].mutation

            elif empty_pos is None:
                print('aaa',self.mutation)
                new_r, new_c = (randint(3)-1, randint(3)-1 )
                if new_r == 0 and new_c == 0:
                    return new_r, new_c

                # center position
                cp_r, cp_c = self.cp
                newcell_pos = (cp_r + new_r, cp_c + new_c)
            #
            #     # push_direction, 0-8 describe the 3*3 matrix，cp_pos = 4 point
            #     if new_r == -1 and new_c == -1:  # 0 point
            #         push_matrix = np.diag(cancer_matrix[:cp_r, :cp_c])
            #         push_matrix = np.delete(push_matrix, [0])
            #         push_matrix = np.append(push_matrix, 0)
            #         if cp_r >= cp_c:
            #             for i in range(cp_c):
            #                 cancer_matrix[cp_r - i, cp_c - i] = push_matrix[i]
            #         else:
            #             for i in range(cp_r):
            #                 cancer_matrix[cp_r - i, cp_c - i] = push_matrix[i]
            #
            #     elif new_r == -1 and new_c == 0:  # 1 point
            #         push_matrix = cancer_matrix[:cp_r, cp_c:(cp_c+1)]
            #         push_matrix = np.delete(push_matrix, [0])
            #         push_matrix = np.append(push_matrix, 0)
            #         for i in range(cp_r):
            #             cancer_matrix[i, cp_c] = push_matrix[i]
            #
            #     # elif new_r == -1 and new_c == 1:  # 2 point
            #     #     push_matrix = []
            #     #     if cp_r >= cp_c:
            #     #         for i in range((Max_COL - cp_c-1)):
            #     #             push_matrix = np.append(push_matrix, cancer_matrix[cp_r - i, cp_c + i])
            #     #             push_matrix = np.delete(push_matrix, [0])
            #     #             push_matrix = np.append(push_matrix, 0)
            #     #             cancer_matrix[cp_r - i, cp_c + i] = push_matrix[i]
            #     #     else:
            #     #         for i in range(Max_ROW - cp_r-2):
            #     #             push_matrix = np.delete(push_matrix, [0])
            #     #             push_matrix = np.append(push_matrix, 0)
            #     #             cancer_matrix[cp_r - i, cp_c + i] = push_matrix[i]
            #
            #
            #     elif new_r == 0 and new_c == -1:  # 3 point
            #         push_matrix = cancer_matrix[cp_r:(cp_r + 1), :cp_c]
            #         push_matrix_pd = pd.DataFrame(push_matrix)
            #         push_matrix_pd.shift(-1, axis='columns')
            #         cancer_matrix[cp_r:(cp_r + 1), :cp_c] = push_matrix_pd.as_matrix()
            #
            #     elif new_r == 0 and new_c == 1:  # 5 point
            #         push_matrix = cancer_matrix[cp_r:(cp_r + 1), cp_c:]
            #         push_matrix_pd = pd.DataFrame(push_matrix)
            #         push_matrix_pd.shift(1, axis='columns')
            #         cancer_matrix[cp_r:(cp_r + 1), cp_c:] = push_matrix_pd.as_matrix()
            #
            #     elif new_r == 1 and new_c == -1:  # 6 point
            #         push_matrix = np.diag(cancer_matrix[cp_r:, :cp_c])
            #         push_matrix_pd = pd.DataFrame(push_matrix)
            #         push_matrix_pd.shift(-1)
            #         cancer_matrix[cp_r:, :cp_c] == push_matrix_pd.as_matrix()
            #
            #     elif new_r == 1 and new_c == 0:  # 7 point
            #         push_matrix = cancer_matrix[cp_r:, cp_c:(cp_c + 1)]
            #         push_matrix_pd = pd.DataFrame(push_matrix)
            #         push_matrix_pd.shift(1)
            #         cancer_matrix[cp_r:, cp_c:(cp_c + 1)] = push_matrix_pd.as_matrix()
            #
            #     elif new_r == 1 and new_c == 1:  # 8 point
            #         push_matrix = np.diag(cancer_matrix[cp_r:, cp_c:])
            #         push_matrix_pd = pd.DataFrame(push_matrix)
            #         push_matrix_pd.shift(1)
            #         cancer_matrix[cp_r:, cp_c:] == push_matrix_pd.as_matrix()
            #
                # daughter_cell = Cancercell(newcell_pos, neighbor_mesh)
                # Creat new daughter cell and it to the cell dictionary
                # daughter_cell.ID = Cancercell(newcell_pos, neighbor_mesh).ID
                # daughter_cell.ID = self.ID + random_choice(range(10))
            #     for i in range(Max_generation):
            #         for i in range(max(self.mutation) + 2*poisson, max(self.mutation) + 3*poisson ):
            #             Cancercell(newcell_pos, neighbor_mesh).mutation = self.mutation.append(i)
            #     daughter_cell.mutation = Cancercell(newcell_pos, neighbor_mesh).mutation+list(range(max(self.mutation) + 1, max(self.mutation) + poisson+1 ))
            #     # mutation process the background_gene hall_mutation_ID lost , so remove it
            #     # agent[newcell_pos].mutation = self.mutation
            #     agent[newcell_pos] = daughter_cell
            #
            #
            #
            #
            # self.die_divide_times -= 1
        '''the cancer cell death process'''
        spontaneous_death = death_r()
        # (add divide death code: 'or self.die_divide_times <= 0')
        if spontaneous_death == 1:
            del agent[self.cp]


if __name__ == "__main__":
    '''# Slicing parameters  
    clipx_star = 20
    clipx_end = 60
    clipy_star = 40
    clipy_end = 60
    N_clip = (clipx_end-clipx_star)*(clipy_end-clipy_star)
    '''
    cancer_grid = cancer_mesh(Max_ROW, Max_COL)
    central_r = int(round(Max_ROW/2))
    central_c = int(round(Max_COL/2))
    central_pos = (central_r, central_c)
    initial_cell = Cancercell(central_pos, cancer_grid)
    cell_dictionary = {central_pos: initial_cell}

    for rep in range(Max_generation):   # replicate cell in the generation
        cell_list = list(cell_dictionary.values())      # get the total cell in list
        shuffle(cell_list)      # random growth order
        for cell in cell_list:  # Iterate through the list to run act
            cell.act(cell_dictionary, cancer_grid)

    cell_r = []
    cell_c = []
    for cell in cell_dictionary.values():   # Iterate through cell_dictionary to find the center point
        c_r, c_c = cell.cp
        c_r = int(c_r)
        c_c = int(c_c)
        cancer_matrix[cell.cp] = cell.ID    # Assign the ID to each center point
        cell_r.append(c_r)
        cell_c.append(c_c)
    cancer_num = Max_ROW * Max_COL - np.sum(cancer_matrix == 0)
    mutation_number = Max_ROW * Max_COL - np.sum(cancer_matrix == 0) - np.sum(cancer_matrix == 1)

    all_mutation_ID = []
    for cut in range(0, len(cell_r)):
        x1 = cell_r[cut]
        y1 = cell_c[cut]
    # for x1 in range(clipx_star,clipx_end ):   # For clipping area information
    #     for y1 in range(clipy_star,clipy_end):
        mutation_ID = cell_dictionary[x1, y1].mutation  # get all the mutation ID in the dictionary
        print('mmm',x1,y1,mutation_ID)
        all_mutation_ID += mutation_ID      # Collect every mutation.ID in the all_mutation_ID

    print(type(cell_dictionary),len(cell_dictionary),'len',len(all_mutation_ID))
    #print('what ID for the pos',type(all_mutation_ID),all_mutation_ID)
    # counter the mutation in the clip area
    result = Counter(all_mutation_ID)       # Count frequency of each mutation ID
    count_mu = []       # empty list: record the mutation ID in the VAF
    count_times = []    # empty list: record the frequency of mutation ID in the VAF
    for i in result:
        count_mu.append(i)
    for j in count_mu:
        count_times.append(result[j])

    # counter the mutation proberbility in the clip area
    VAF = list(map(lambda n:n/(2*cancer_num), count_times))     # Object function change to list, and calculate each element of count_times
    result_VAF = Counter(VAF)   # Count each element in VAF, feedback the numbers
    VAF_mu = []     # empty list: record the mutation ID in the VAF
    VAF_times = []      # empty list: record the each ID frequency in the VAF
    for i1 in result_VAF:
        VAF_mu.append(i1)
    for j1 in VAF_mu:
        VAF_times.append(result_VAF[j1])

    end_time = time.time()
    run_time = end_time - start_time

    #VAF = [(count_times_element - 1) for count_times_element in range(count_times)] #variant allele frequence
    print(result, '\n', 'counter mutation: ', count_mu, '\n', len(result), result[1], 'counter times: ', count_times)
    print(VAF, result_VAF, 'VAF mutaion', VAF_mu, 'VAF_times', VAF_times)
    print( 'Growth Run time is: {}'.format(run_time), '\n','Max Mutation',max(all_mutation_ID),'\n','mutation ID number per cell', len(all_mutation_ID)/cancer_num, '\n','generation times:{}'.format(Max_generation),'\n','Cancer Number is: {}'.format(cancer_num),'\n','Mutation nember: ',mutation_number )
    # print('cell', cell_dictionary[100, 101].mutation, '\n', cell_dictionary[99, 101].mutation, '\n',  '[99, 99]',cell_dictionary[99, 99].mutation, '\n', '[111, 111]',cell_dictionary[111, 111].mutation,'\n',cell_dictionary[110, 111].mutation, '\n',cell_dictionary[112, 111].mutation)
    '''
    # mutant_total = []
    # for m_t in range(0, (100)):
    #     mutant_total.append(cell.mutation)
    # print('mutant_total is', mutant_total)
    # p = 2
    # b = len(mutant_total)
    # plt.subplot(131)
    # plt.bar(count_mu, count_times, width=5, color='seagreen')
    '''
    plt.subplot(121)
    plt.bar(VAF_mu, VAF_times, width=0.005, color='seagreen')
    plt.xlabel("Variant Allele Frequency")
    plt.ylabel("Number of Mutation")

    '''
    poisson_dis = poisson()
    for a in poisson_dis:
        b = sum(poisson_dis == a)
        plt.subplot(131)
        plt.bar(a, b, color='red')
        plt.title('cancer number:{}'.format(cell_num))
        plt.xlabel("mutation class")
        plt.ylabel("Number of mutation")
    '''
    print('mutant_total is',type(count_mu), type(count_times),)
    plt.subplot(122)
    map = matplotlib.colors.ListedColormap(['#FFFFFF', '#0000FF', '#6495ED', '#FF7F50','#ff8d00','#006400',
    '#8FBC8F','#9400D3','#00BFFF','#FFD700','#CD5C5C', '#556B2F', '#9932CC', '#8FBC8F','#2F4F4F',
    '#CD853F', '#FFC0CB', '#4169E1', '#FFFFE0', '#ADD8E6', '#008000','#9400D3',
    '#9ACD32', '#D2B48C', '#008B8B','#9400D3', '#00BFFF',  '#CD5C5C',])
    plt.title('Cancer Number is: {}'.format(cancer_num))
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    plt.imshow(cancer_matrix, cmap = map)
    plt.axis('off')
    plt.show()
    plt.ion()





