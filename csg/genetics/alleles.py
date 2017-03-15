#!/usr/bin/env python
import pandas as pd
import numpy as np

def flip_strand(x):
  if x == "A":
    return "T"
  elif x == "T":
    return "A"
  elif x == "G":
    return "C"
  elif x == "C":
    return "G"
  else:
    raise Exception("No matching allele for: %s" % str(x))

def flip_effects(a1,a2,o1,o2,effect):
  """
  This function flips effect sizes (effect) specified relative to a given set of alleles (o1, o2) towards
  a target set of alleles (a1, a2). We assume the effects are relative to o1. 
  
  Args: 
    a1 (Series or ndarray): Allele 1 
    a2 (Series or ndarray): Allele 2
    o1 (Series or ndarray): Other allele 1 
    o2 (Series or ndarray): Other allele 2 
    effect (Series[float] or ndarray[float]): Effect sizes
  
  Returns:
    Nothing, arrays are modified in place (o1,o2,effect are modified)

    If alleles cannot be matched at all, the effect will be nan'd
  """
  
  assert len(a1) == len(a2) == len(o1) == len(o2) == len(effect)
  nan = pd.np.nan
    
  # Use the underlying numpy array if we can (it is orders of magnitude faster.)
  # This is a total hack, though. 
  try:
    a1 = a1.values
    a2 = a2.values
    o1 = o1.values
    o2 = o2.values
    effect = effect.values
  except:
    pass

  #for c in "a1 a2 o1 o2 effect".split():
    #if isinstance(eval(c),pd.Series):
      #exec "{v} = {v}.values".format(v=c)
    #elif isinstance(eval(c),np.ndarray):
      #pass
    #else:
      #raise ValueError, "Argument {} is not a pd.Series or np.ndarray".format(c)
  
  # Quick check that alleles are upper case. 
  if a1[0] != a1[0].upper() or o1[0] != o1[0].upper():
    raise ValueError("Alleles should be upper case")
  
  for i in range(len(a1)):
    i_a1 = a1[i]
    i_a2 = a2[i]
    i_o1 = o1[i]
    i_o2 = o2[i]
    
    # Do the alleles match? 
    match = set([i_a1,i_a2]) == set([i_o1,i_o2])
    if not match:
      # Hmm. How about strand flip? 
      flip_o1 = flip_strand(i_o1)
      flip_o2 = flip_strand(i_o2)
      
      flip_match = set([i_a1,i_a2]) == set([flip_o1,flip_o2])
      if not flip_match:
        # We couldn't even match the alleles with a strand flip. 
        # This one is bogus!
        effect[i] = nan
        continue
      else:
        # Strand flipping resulted in matching alleles
        i_o1 = flip_o1
        o1[i] = flip_o1
        
        i_o2 = flip_o2
        o2[i] = flip_o2
        
    # If we're here, alleles are now matched. But do we need to swap them? 
    if i_a1 != i_o1:
      # Swap is needed. 
      i_o1, i_o2 = i_o2, i_o1
      o1[i] = i_o1
      o2[i] = i_o2
      effect[i] = -effect[i]

def _tests():
  test_allele_flip = pd.DataFrame(
    [
      [1,"A","G","G","A",-2,"Swap alleles, invert effect"],
      [2,"C","T","C","T",4.3,"No action required"],
      [3,"G","A","C","T",5.5,"Strand flip, no invert effect"],
      [4,"G","A","T","C",-1.3,"Strand flip and swap, invert effect"],
      [5,"G","T","T","A",-4,"Cannot resolve, alleles do not match even with strand flip, effect should be set to nan"]
    ],
    columns = [
      "test_id",
      "allele1",
      "allele2",
      "other1",
      "other2",
      "effect",
      "comment"
    ]
  )

  test_allele_flip.insert(5,"original_other1",test_allele_flip["other1"])
  test_allele_flip.insert(6,"original_other2",test_allele_flip["other2"])
  test_allele_flip.insert(7,"original_effect",test_allele_flip["effect"])
        
  large_flip_speed_test = pd.DataFrame(
    {
      "allele1": np.random.choice(["A","C","G","T"],1000),
      "allele2": np.random.choice(["A","C","G","T"],1000),
      "other1": np.random.choice(["A","C","G","T"],1000),
      "other2": np.random.choice(["A","C","G","T"],1000),
      "effect": np.random.uniform(-2,2,1000)
    }
  )

  flip_effects(
    test_allele_flip.allele1,
    test_allele_flip.allele2,
    test_allele_flip.other1,
    test_allele_flip.other2,
    test_allele_flip.effect
  )

