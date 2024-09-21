#!/usr/bin/env python
# coding: utf-8

# # N-Queens problem using genetic algorithm
# # Maryam Hosseinali 610398209

# ## general explanation:
# ### In this implementation, the fitness score of each chromosome which is a solution candidate, is calculated based on the number of conflicts between queens in the same row, same main diagonal, and same off-diagonal. The lower the score, the better the solution.
# 

# ## import required libraries

# In[21]:


import numpy as np
import random
import copy


# ## Gene class:
# ### represents a single gene in the chromosome, its value representing the position of a queen on the chessboard.

# In[22]:


class Gene:
    def __init__(self, _value):
        self.value = _value
    
    def __str__(self):
        return str(self.value)


# ## Chromosome class:
# 
# ### Score and Fitness:
# ### represents a solution candidate, consists of N genes with methods for calculating the fitness score of a chromosome, performing crossover and mutation operations, and getting a string representation of the chromosome.
# ### in general we have 4 types of conflicts:
#     
# ### 1: two queens in the same row
# ### 2: two queens in the same column
# ### 3: two queens with main diagonal conflict
# ### 4: two queens with off-diagonal conflict
# ### - the first case does not occure because of the choromosome definition (in each column exactly on queen)
# ### - the second case does not occure either cause of definition of cross-over and mutation (but we consider in case of functions later changes)
# 
# ### calculate_score function will calculate the score for a chromosome using score_func function.
# ### for example for the first case counts the repeated numbers in a gene or for the third case counts the condtions that the difference between the value and its index is equal.
# 
# ### the total number of conflicts calculated as fallows:
# ### for each possible case we discussed, we calculate the sum of C(n, 2) for each n.( C(n,2) is n choose 2) in other words it calculates the number of pairs of queens that are attacking each other.
# ### for clarity if in the third case we find three number 3 it means we have 3 queens who threaten each other on the main diagonal and we should add C(3,2) to our number of conflicts.
# 
# 
# ### get_fitness_score function calculates and returns the total fitness score of the chromosome based on its row, main diagonal and off-diagonal scores.
# 

# ### Crossover: unified_crossover takes two chromosomes and a crossover probability(p_c) and produces two new chromosomes by exchanging some of their genes.
# ### it's remarkable that in a usual and simple crossover, to produce new off_springs from parent A and B with probability p_c, some part of parents swap with another but in unified_crossover we iterate from first bit to last bit (using for loops) and with crossover probability p_c their i_th bit may change with each other(i_th bit of parent A swap with i_th bit of parent B). in fact in this mode of crossover(unified) we may have crossover for each bit. this makes us not have repeated number in a chromosme and the process becomes faster. 

# ### Mutation: swap_mutate takes a chromosome and a mutation probability (p_m) or mutation rate and randomly swaps two genes within the chromosome.

# In[23]:


class Chromosome:
    
# this function creates an array of genes of length _length, with values ranging from 1 to _length, 
# and then shuffles the genes randomly.
    def __init__(self, _length):
        self.length = _length
        self.genes = np.array([Gene(_value=i) for i in range(1, _length+1)])
        np.random.shuffle(self.genes)
        
#row_score_func takes a gene value and index and returns the gene value. 
    @staticmethod
    def row_score_func(gene_value, gene_ind):
        return gene_value
    
#main_diag_score_func takes a gene value and index 
#and returns the difference between the gene value and its index.  
    @staticmethod
    def main_diag_score_func(gene_value, gene_ind):
        return gene_value - gene_ind
    
    
#off_diag_score_func takes a gene value and index 
#and returns the sum of the gene value and its index.    
    @staticmethod
    def off_diag_score_func(gene_value, gene_ind):
        return gene_value + gene_ind

    
#calculate the score for a chromosome using a scoring function.
#(the process is described above)
    def calculate_score(self, score_func):
        score = 0
        all_genes_values_cnt = {value: 0 for value in range(1-self.length, 2 * self.length + 1)}
        for gene_ind in range(len(self.genes)):
            gene = self.genes[gene_ind]
            all_genes_values_cnt[score_func(gene.value, gene_ind)] += 1
        for gene_value in all_genes_values_cnt.keys():
            gene_value_cnt = all_genes_values_cnt[gene_value]
            score += (gene_value_cnt * (gene_value_cnt - 1)) // 2 
        return score
    
 #calculate and return the total fitness score of choromosome
#by its row, main diagonal and off-diagonal scores
    def get_fitness_score(self):
        row_score = self.calculate_score(self.row_score_func)
        main_diag_score = self.calculate_score(self.main_diag_score_func)
        off_diag_score = self.calculate_score(self.off_diag_score_func)

        return row_score + main_diag_score + off_diag_score
    
    def unified_crossover(self, other_chromosome, p_c):
        #iterates over the genes of the chromosome. 
        for gene_ind in range(self.length):
            #checks whether a randomly generated float value between 0 and 1
            #is less than or equal to the crossover probability p_c.
            if random.random() <= p_c:
            #to keep track of the index of the first and second genes with
            #the same value as the gene being crossed over.
                first_ind = 0
                second_ind = 0
                
                #iterates over the genes of the other chromosome.
                for ind in range(self.length):
                    
                    #checks whether the value of the gene in the other chromosome at
                    #index ind is equal to the value of the gene in the current chromosome at index gene_ind.
                    if self.genes[ind].value == other_chromosome.genes[gene_ind].value:
                    #sets the variable first_ind to the current value of ind.
                    #updates first_ind to the index of the first gene with the same value as the gene being crossed over.
                        first_ind = ind
                    if other_chromosome.genes[ind].value == self.genes[gene_ind].value:
                        second_ind = ind
                
                tmp = self.genes[gene_ind].value
                self.genes[gene_ind].value = other_chromosome.genes[gene_ind].value
                other_chromosome.genes[gene_ind].value = tmp
                
                tmp = self.genes[first_ind].value
                self.genes[first_ind].value = other_chromosome.genes[second_ind].value
                other_chromosome.genes[second_ind].value = tmp
  
 #randomly swaps two genes within the chromosome.
    def swap_mutate(self, p_m):
        if random.random() <= p_m:
            first_gene_ind = random.randint(0, self.length-1)
            second_gene_ind = random.randint(0, self.length-1)
            
            tmp = self.genes[first_gene_ind].value
            self.genes[first_gene_ind].value = self.genes[second_gene_ind].value
            self.genes[second_gene_ind].value = tmp
            
            
    def __str__(self):
        result = "([" + " ".join(map(str, self.genes.tolist())) + "], Score: " +                                 str(self.get_fitness_score()) + ")"
        return result


# ## Population class:
# ### Population class is for creating, shuffling, sorting, and performing crossover and mutation operations on a population of chromosomes.
# ### initialize a population by creating a numpy array of Chromosome objects and then randomly generating chromosomes for the initial population.(by shuffle and numpy)
# ### - crossover_operator method performs the crossover operation on the population with a probability p_c.by selecting pairs of adjacent chromosomes in the population and applying the unified crossover operator on them.
# ### - mutatation_operator method performs the mutation operation on the population with a probability p_m. by applying the swap mutation operator on each chromosome in the population.
# ### - sort_by_fitness method sorts the chromosomes in the population in descending order based on their fitness score, which is calculated by calling the get_fitness_score on each chromosome.

# In[24]:


class Population:
    
#initial a population.
#length of the array is the population size and
#the length of each chromosome is chromosome_length parameter

    def __init__(self, population_length, chromosome_length):
        self.length = population_length
        self.chromosomes = np.array([Chromosome(_length=chromosome_length) for _ in range(population_length)])
        

#shuffles the order of the chromosomes in the population randomly 
    def shuffle(self):
        np.random.shuffle(self.chromosomes)
        
        
        
#crossover operation on the population with p_c
    def crossover_operator(self, p_c):
        for chromosome_ind in range(3, self.length - 1):
            self.chromosomes[chromosome_ind].unified_crossover(self.chromosomes[chromosome_ind + 1], p_c)

# mutation operation on the population with          
    def mutatation_operator(self, p_m):
        for chromosome_ind in range(self.length - 1):
            self.chromosomes[chromosome_ind].swap_mutate(p_m)
#orts the chromosomes in descending order based on their fitness score    
    def sort_by_fitness(self):
        self.chromosomes = np.array(sorted(self.chromosomes, key=lambda x: x.get_fitness_score()))
    

#returns a string representation of the population
    def __str__(self):
        result = ", ".join(map(str, self.chromosomes.tolist()))
        return result


# ## Solver class:
# ### solving the N-Queens problem using a genetic algorithm. It initializes a population of chromosomes with random values, where each chromosome represents a possible solution to the problem. It then repeatedly generates new populations by selecting the best individuals from the previous population, applying crossover and mutation operators to create new individuals, and sorting the resulting population by fitness score.
# 
# ### - rank_selection method performs rank-based selection, where each individual in the population is assigned a rank based on its fitness score, and a random individual is selected with a probability proportional to its rank. it helps to ensure that the fittest individuals are more likely to be selected for the next generation.
# ### - is_goal_state method checks whether the best individual in the population has a fitness score of 0, which is a valid solution to the N-Queens problem.
# ### - The generate_new_population method applies the genetic operators to create a new population, including crossover to create new individuals by combining the genetic material of two parents, and mutation to introduce random changes to an individual's genetic material.
# 
# ### and Finally the solve method repeatedly generates new populations until a valid solution is found, and prints out the best individual in the final population.

# In[34]:


class Solver:
    def __init__(self, n_value):
        self.POP_LENGTH = 300
        self.P_C = 0.005
        self.P_M = 0.2
        
        self.population = Population(population_length = self.POP_LENGTH, chromosome_length=n_value)
    
    def rank_selection(self):
        total = (self.population.length * (self.population.length + 1)) // 2
        population_backup = self.population.chromosomes.copy()
        for chromosome_ind in range(self.population.length):
            select_value = random.randint(1 , total)
            cur_chromosome_rank = self.population.length
            select_ind = 0
            while cur_chromosome_rank != 0:
                if select_value <= cur_chromosome_rank:
                    self.population.chromosomes[chromosome_ind] = copy.deepcopy(population_backup[select_ind])
                    break
                else:
                    select_value -= cur_chromosome_rank
                    select_ind += 1
                    cur_chromosome_rank -= 1
            
    
    def is_goal_state(self):
        best_chromosome = self.population.chromosomes[0]
        return best_chromosome.get_fitness_score() == 0
    
    def generate_new_population(self):
        self.rank_selection()
        self.population.shuffle()
        self.population.crossover_operator(self.P_C)
        self.population.mutatation_operator(self.P_M)
        self.population.sort_by_fitness()
    
    def solve(self):
        self.population.sort_by_fitness()
        while not self.is_goal_state():
            self.generate_new_population()
            
            print("total conflicts so far:")
            print(self.population.chromosomes[0].get_fitness_score())
        print("solution found: ")
        print(self.population.chromosomes[0])


# In[43]:


solver = Solver(100)


# In[44]:


solver.solve()

