import copy
import numpy
from collections import Counter

spectrum = [57, 57, 71, 99, 129, 137, 170, 186, 194, 208, 228, 265, 285, 299, 307, 323, 356, 364, 394, 422, 493]
M = 20
N = 60
# file = '/content/dataset_104_7.txt'
# with open(file) as fh:
#   data = fh.readlines()
#   M = int(data[0])
#   N = int(data[1])
#   spectrum_str = data[2].split()
#   spectrum = [int(x) for x in spectrum_str]
spectrum.sort()

def extended_mass_table():
  """Creates a dictionary using a unicode as the key and the mass as the value"""
  extended_list = {}
  for i in range(57,201):
      extended_list[chr(i)] = int(i)
  return(extended_list)

def spectral_convolution(spectrum):
  """Creates a dictionary by using the convolution differences as the key and
     using the number of occurances as the value"""
  convolution = []
  for i in range(len(spectrum)-1):
    for j in range(i+1, len(spectrum)):
      diff = spectrum[j] - spectrum[i]
      if diff == 0 or diff > 200 or diff < 57:
        pass
      else:
        convolution.append(diff)
  convolution_dict = Counter(convolution)
  return convolution_dict

def best_convolutions(convolution_dict, M):
  """Takes the convolutions with the M most occurances and returns them as a
     list"""
  best_mass = numpy.array(list(convolution_dict.keys()))
  scores = numpy.array(list(convolution_dict.values()))
  indx = scores.argsort()[::-1] # sorts the scores indices in descending order
  best_mass = list(best_mass[indx])
  scores = list(scores[indx])
  for j in range(M, len(best_mass)):
    if scores[j] < scores[M-1]:
      del best_mass[j:]
      return best_mass
  return best_mass

def mass_to_aa(top_mass, mass_initial):
  """Converts a list of masses to a dictionary where the key is the respective
     unicode and the mass is the value"""
  masses = {}
  for mass in top_mass:
    for aa in mass_initial:
      if mass_initial[aa] == mass:
        masses[aa] = mass
        break
  return masses
    
def leaderboard_sequencing(spectrum, N):
  """Takes in a list spectrum and an int N. A leaderboard is made up of spectrums
     and at the last iteration the highest scoring spectrum compared to the given
     spectrum is returned"""
  leaderboard = {""}
  leader_p = ''
  convolution = spectral_convolution(spectrum)
  top_mass = best_convolutions(convolution, M)
  masses = mass_to_aa(top_mass, mass_initial)
  while len(leaderboard) != 0:
    leaderboard = expand(leaderboard, masses)
    temp_board = copy.deepcopy(leaderboard)
    for peptide in temp_board:
      if total_mass(peptide, masses) == parent_mass(spectrum):
        if score(peptide, masses) > score(leader_p, masses):
          leader_p = peptide
      elif total_mass(peptide, masses) > parent_mass(spectrum):
        leaderboard.remove(peptide)
    leaderboard = trim(leaderboard, spectrum, masses, N)
  leader_p = convert_str(leader_p, masses)
  return leader_p

def trim(leaderboard:set, spectrum:list, masses:dict, N:int):
  """Takes in a set, leaderboard, list, spectrum, dictionary, masses, int, N
     and trims the leaderboard based on the N highest scoring spectrums unless
     there are ties in scores with the N highest score, once there is a spectrum
     that has a score lower than the N highest score the new leaderboard is
     returned"""
  leaderboard = list(leaderboard)
  scores = []
  for j in range(len(leaderboard)):
    peptide = leaderboard[j]
    scores.append(lin_score(peptide, masses))
  leaderboard = numpy.array(leaderboard)
  scores = numpy.array(scores)
  indx = scores.argsort()[::-1] # sorts the scores indices in descending order
  leaderboard = list(leaderboard[indx])
  scores = list(scores[indx])
  for j in range(N, len(leaderboard)):
    if scores[j] < scores[N-1]:
      del leaderboard[j:]
      return leaderboard
  return leaderboard

def score(peptide, masses):
  """Compares a cyclopeptide spectrum to the given spectrum and returns the
     score"""
  theoretical = cyc_spectrum(peptide, masses)
  count = 0
  for num in spectrum:
    if num in theoretical:
      count += 1
      theoretical.remove(num)
  return count

def lin_score(peptide, masses):
  """Compares a linear peptide spectrum to the given spectrum and returns the
     score"""
  consistent = lin_spectrum(peptide, masses)
  count = 0
  for num in spectrum:
    if num in consistent:
      count += 1
      consistent.remove(num)
  return count

def expand(candidate_peptides, masses):
  """Iterates through each peptide adding each letter in masses and returns the
     new peptides"""
  peptides = set()  
  for peptide in candidate_peptides:  
    for letter in masses:
      new_peptide = peptide + letter
      peptides.add(new_peptide)
  return peptides

def total_mass(peptide, masses):
  """Returns the total mass of the peptide passed in"""
  total = 0
  for aa in peptide:
    total += masses[aa]
  return total

def parent_mass(spectrum):
  """Returns the largest mass in the spectrum"""
  largest = spectrum[-1]
  return largest

def convert_str(peptide, masses):
  """Converts peptide to a representation of masses"""
  nums = ''
  for i,letter in enumerate(peptide):
    num = str(masses[letter])
    if i == 0:
      nums += num
    else:
      nums += '-' + num
  return nums

def lin_spectrum(peptide, masses):
  """Reutrns the linear spectrum for a peptide"""
  prefix_mass = [0]
  for i in range(len(peptide)):
    letter = peptide[i]
    prefix_current = prefix_mass[i] + masses[letter]
    prefix_mass.append(prefix_current)
  linear_spectrum = [0]
  for i in range(len(prefix_mass)-1):
    for j in range(i+1, len(prefix_mass)):
      current_mass = prefix_mass[j] - prefix_mass[i]
      linear_spectrum.append(current_mass)
  linear_spectrum.sort()
  return linear_spectrum

def cyc_spectrum(peptide, masses):
  """Returns the cyclo spectrum for a peptide"""
  prefix_mass = [0]
  for i in range(len(peptide)):
    letter = peptide[i]
    prefix_current = prefix_mass[i] + masses[letter]
    prefix_mass.append(prefix_current)
  peptide_mass = prefix_mass[len(peptide)]
  cyclic_spectrum = [0]
  for i in range(len(prefix_mass)-1):
    for j in range(i+1, len(prefix_mass)):
      current_mass = prefix_mass[j] - prefix_mass[i]
      cyclic_spectrum.append(current_mass)
      if i > 0 and j < len(prefix_mass)-1:
        current_mass = peptide_mass - (prefix_mass[j] - prefix_mass[i])
        cyclic_spectrum.append(current_mass)  
  cyclic_spectrum.sort()
  return cyclic_spectrum

mass_initial = extended_mass_table()
leaderboard_sequencing(spectrum, N)