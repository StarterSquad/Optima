import add_optima_paths
import program


# Let's say modality 1 can reach 100% of the population, but modalities 2 and 3 can reach 50%. Then
# (1) - 1
# (2) - 0.5
# (3) - 0.5
# (1,2) = 1*0.5
# (1,3) = 1*0.5
# (2,3) = 0.5*0.5 = 0.25
# (1,2,3) = 1*0.5*0.5 = 0.25

# This includes double counting
# Now, we remove the double counting
# The unique coverage we want is
# (1) - 0.25
# (2) - 0
# (3) - 0
# (1,2) = 0.25
# (1,3) = 0.25
# (2,3) = 0
# (1,2,3) = 0.25

# Algorithmically, we go
# 0.25 people are in (1,2,3), so they get subtracted from any population containing them


# Then we have 50% of the population reached only by modality 1
# 


r =  program.Program('Increased testing')
m = []
m.append(r.add_modality('TV ads',1))
m.append(r.add_modality('Mobile clinic',0.5))
m.append(r.add_modality('GP flyers',0.5))
m.append(r.add_modality('GP flyers2',0.5))

# Clearly if we have 4 programs, we have
# (1,2,3,4) = 0.0625
# (1,2,3) = 0.12 -> 0.0625
# (2,3,4) = 0.12 -> 0.0625
# (1,3,4) = 0.12 -> 0.0625
# Now consider (1,3)

print r.modalities
print r.metamodalities
#p = programs.Program('Test program')



