#!/usr/bin/env python3
import ast, patsy

def get_variables(formula):
  """
  Given a patsy formula, extract out the names of all variables (even ones within function calls.)
  """

  names = set()
  funcs = set()

  def traverse(node):
    for n in ast.iter_child_nodes(node):
      if isinstance(n,ast.Call):
        funcs.add(n.func.id)
      elif isinstance(n,ast.Name):
        names.add(n.id)
      elif isinstance(n,ast.Str):
        names.add(n.s)

      traverse(n)

  desc = patsy.ModelDesc.from_formula(formula)

  for term in desc.rhs_termlist:
    for factor in term.factors:
      traverse(ast.parse(factor.code))
  
  for term in desc.lhs_termlist:
    for factor in term.factors:
      traverse(ast.parse(factor.code))

  variables = names.difference(funcs)
  return variables

