import unittest
from IntervalTree import *

class TestIntervalTreeNode(unittest.TestCase):

   def test_node_construction(self):
      left_child = IntervalTreeNode(0, 1, 1, 'value', 0, IntervalTreeNode.RED, None, None)
      self.assertEqual(left_child.start, 0)
      self.assertEqual(left_child.end, 1)
      self.assertEqual(left_child.max_end, 1)
      self.assertEqual(len(left_child.values), 1)
      self.assertEqual(left_child.values[0], 'value')
      self.assertEqual(left_child.deviation, 0)
      self.assertEqual(left_child.color, IntervalTreeNode.RED)
      self.assertIsNone(left_child.left)
      self.assertIsNone(left_child.right)

      right_child = IntervalTreeNode(0, 1, 1, 'value', 0, IntervalTreeNode.RED, None, None)

      parent = IntervalTreeNode(0, 1, 1, 'value', 0, IntervalTreeNode.RED, left_child, right_child)
      self.assertEqual(id(left_child.parent), id(parent))
      self.assertEqual(id(right_child.parent), id(parent))
      self.assertEqual(id(parent.left), id(left_child))
      self.assertEqual(id(parent.right), id(right_child))
      self.assertEqual(id(left_child), id(right_child.get_sibling()))
      self.assertEqual(id(right_child), id(left_child.get_sibling()))

      uncle = IntervalTreeNode(0, 1, 1, 'value', 0, IntervalTreeNode.RED, None, None)
      grandparent = IntervalTreeNode(0, 1, 1, 'value', 0, IntervalTreeNode.RED, parent, uncle)
      self.assertEqual(id(left_child.get_uncle()), id(uncle))
      self.assertEqual(id(right_child.get_uncle()), id(uncle))
      self.assertEqual(id(left_child.get_grandparent()), id(grandparent))
      self.assertEqual(id(right_child.get_grandparent()), id(grandparent))      


   def test_node_compare(self):
      node1 = IntervalTreeNode(-2, 2, 2, 'value', 0, IntervalTreeNode.RED, None, None)
      
      node2 = IntervalTreeNode(-4, -3, -3, 'value', 0, IntervalTreeNode.RED, None, None)
      self.assertEqual(node2.compare(node1), -1)

      node2 = IntervalTreeNode(-3, -2, -2, 'value', 0, IntervalTreeNode.RED, None, None)
      self.assertEqual(node2.compare(node1), -1)

      node2 = IntervalTreeNode(-2, -1, -1, 'value', 0, IntervalTreeNode.RED, None, None)
      self.assertEqual(node2.compare(node1), -1)

      node2 = IntervalTreeNode(-2, 2, 2, 'value', 0, IntervalTreeNode.RED, None, None)
      self.assertEqual(node2.compare(node1), 0)

      node2 = IntervalTreeNode(-2, 3, 3, 'value', 0, IntervalTreeNode.RED, None, None)
      self.assertEqual(node2.compare(node1), 1)

      node2 = IntervalTreeNode(-1, 1, 1, 'value', 0, IntervalTreeNode.RED, None, None)
      self.assertEqual(node2.compare(node1), 1)

      node2 = IntervalTreeNode(-1, 2, 2, 'value', 0, IntervalTreeNode.RED, None, None)
      self.assertEqual(node2.compare(node1), 1)

      node2 = IntervalTreeNode(-1, 3, 3, 'value', 0, IntervalTreeNode.RED, None, None)
      self.assertEqual(node2.compare(node1), 1)

      node2 = IntervalTreeNode(2, 3, 3, 'value', 0, IntervalTreeNode.RED, None, None)
      self.assertEqual(node2.compare(node1), 1)

      node2 = IntervalTreeNode(3, 4, 4, 'value', 0, IntervalTreeNode.RED, None, None)
      self.assertEqual(node2.compare(node1), 1)

      node2 = IntervalTreeNode(-3, 3, 1, 'value', 0, IntervalTreeNode.RED, None, None)
      self.assertEqual(node2.compare(node1), -1)



class TestIntervalTree(unittest.TestCase):
          
   def test_tree_construction(self):   
      tree = IntervalTree()
      self.assertIsNone(tree.root)
      self.assertEqual(tree.get_intervals_count(), 0)

   def test_tree_balance_case1(self):
      tree = IntervalTree()
      
      tree.add(-2, 2, '1')
      self.assertIsNotNone(tree.root)
      self.assertEqual(tree.root.color, IntervalTreeNode.BLACK) 
      self.assertEqual(tree.root.max_end, 2)
      self.assertEqual(tree.get_intervals_count(), 1)
      self.assertEqual(tree.get_height(), 1)
      
      tree.add(-3, 1, '3')
      self.assertEqual(tree.root.color, IntervalTreeNode.BLACK) 
      self.assertEqual(tree.root.max_end, 2) 
      self.assertEqual(tree.get_intervals_count(), 2)
      self.assertEqual(tree.get_height(), 2)

      tree.add(-4, 0, '4')
      self.assertEqual(tree.root.color, IntervalTreeNode.BLACK) 
      self.assertEqual(tree.root.start, -3) 
      self.assertEqual(tree.root.end, 1) 
      self.assertEqual(tree.root.max_end, 2) 
      self.assertEqual(tree.get_intervals_count(), 3)
      self.assertEqual(tree.get_height(), 2)


   def test_tree_balance_case2(self):
      tree = IntervalTree()
      
      tree.add(-2, 2, '1')
      self.assertIsNotNone(tree.root)
      self.assertEqual(tree.root.color, IntervalTreeNode.BLACK) 
      self.assertEqual(tree.root.max_end, 2)
      self.assertEqual(tree.get_intervals_count(), 1)
      self.assertEqual(tree.get_height(), 1)
      
      tree.add(-4, 0, '4')
      self.assertEqual(tree.root.color, IntervalTreeNode.BLACK) 
      self.assertEqual(tree.root.max_end, 2) 
      self.assertEqual(tree.get_intervals_count(), 2)
      self.assertEqual(tree.get_height(), 2)

      tree.add(-3, 1, '3')
      self.assertEqual(tree.root.color, IntervalTreeNode.BLACK) 
      self.assertEqual(tree.root.start, -3) 
      self.assertEqual(tree.root.end, 1) 
      self.assertEqual(tree.root.max_end, 2) 
      self.assertEqual(tree.get_intervals_count(), 3)
      self.assertEqual(tree.get_height(), 2)


   def test_tree_balance_case3(self):
      tree = IntervalTree()
      
      tree.add(-2, 2, '1')
      self.assertIsNotNone(tree.root)
      self.assertEqual(tree.root.color, IntervalTreeNode.BLACK) 
      self.assertEqual(tree.root.max_end, 2)
      self.assertEqual(tree.get_intervals_count(), 1)
      self.assertEqual(tree.get_height(), 1)
      
      tree.add(1, 3, '4')
      self.assertEqual(tree.root.color, IntervalTreeNode.BLACK) 
      self.assertEqual(tree.root.max_end, 3) 
      self.assertEqual(tree.get_intervals_count(), 2)
      self.assertEqual(tree.get_height(), 2)

      tree.add(2, 4, '3')
      self.assertEqual(tree.root.color, IntervalTreeNode.BLACK) 
      self.assertEqual(tree.root.start, 1) 
      self.assertEqual(tree.root.end, 3) 
      self.assertEqual(tree.root.max_end, 4) 
      self.assertEqual(tree.get_intervals_count(), 3)
      self.assertEqual(tree.get_height(), 2)


   def test_tree_balance_case4(self):
      tree = IntervalTree()
      
      tree.add(-2, 2, '1')
      self.assertIsNotNone(tree.root)
      self.assertEqual(tree.root.color, IntervalTreeNode.BLACK) 
      self.assertEqual(tree.root.max_end, 2)
      self.assertEqual(tree.get_intervals_count(), 1)
      self.assertEqual(tree.get_height(), 1)
      
      tree.add(2, 4, '3')
      self.assertEqual(tree.root.color, IntervalTreeNode.BLACK) 
      self.assertEqual(tree.root.max_end, 4) 
      self.assertEqual(tree.get_intervals_count(), 2)
      self.assertEqual(tree.get_height(), 2)

      tree.add(1, 3, '4')
      self.assertEqual(tree.root.color, IntervalTreeNode.BLACK) 
      self.assertEqual(tree.root.start, 1) 
      self.assertEqual(tree.root.end, 3)  
      self.assertEqual(tree.root.max_end, 4) 
      self.assertEqual(tree.get_intervals_count(), 3)
      self.assertEqual(tree.get_height(), 2)


   def test_tree_balance_case5(self):
      tree = IntervalTree()
      
      tree.add(-2, 2, '1')
      self.assertIsNotNone(tree.root)
      self.assertEqual(tree.root.color, IntervalTreeNode.BLACK) 
      self.assertEqual(tree.root.max_end, 2)
      self.assertEqual(tree.get_intervals_count(), 1)
      self.assertEqual(tree.get_height(), 1)
      
      tree.add(-3, 1, '3')
      self.assertEqual(tree.root.color, IntervalTreeNode.BLACK) 
      self.assertEqual(tree.root.max_end, 2) 
      self.assertEqual(tree.get_intervals_count(), 2)
      self.assertEqual(tree.get_height(), 2)

      tree.add(0, 4, '4')
      self.assertEqual(tree.root.color, IntervalTreeNode.BLACK) 
      self.assertEqual(tree.root.start, -2) 
      self.assertEqual(tree.root.end, 2)  
      self.assertEqual(tree.root.max_end, 4) 
      self.assertEqual(tree.get_intervals_count(), 3)
      self.assertEqual(tree.get_height(), 2)

      tree.add(1, 5, '5')
      self.assertEqual(tree.root.color, IntervalTreeNode.BLACK)
      self.assertEqual(tree.root.start, -2) 
      self.assertEqual(tree.root.end, 2)  
      self.assertEqual(tree.root.max_end, 5) 
      self.assertEqual(tree.get_intervals_count(), 4)
      self.assertEqual(tree.get_height(), 3)
     
      tree.add(-1, 6, '6')
      self.assertEqual(tree.root.color, IntervalTreeNode.BLACK)
      self.assertEqual(tree.root.start, -2) 
      self.assertEqual(tree.root.end, 2)   
      self.assertEqual(tree.root.max_end, 6)
      self.assertEqual(tree.get_intervals_count(), 5)
      self.assertEqual(tree.get_height(), 3)
     
      tree.add(2, 7, '6')
      self.assertEqual(tree.root.color, IntervalTreeNode.BLACK)
      self.assertEqual(tree.root.start, -2) 
      self.assertEqual(tree.root.end, 2)   
      self.assertEqual(tree.root.max_end, 7)
      self.assertEqual(tree.get_intervals_count(), 6)
      self.assertEqual(tree.get_height(), 4)



   def test_duplicated_intervals(self):
      tree = IntervalTree()
      
      tree.add(-2, 2, '1')
      self.assertIsNotNone(tree.root)
      self.assertEqual(tree.root.color, IntervalTreeNode.BLACK) 
      self.assertEqual(tree.root.max_end, 2)
      self.assertEqual(tree.get_intervals_count(), 1)
      self.assertEqual(tree.get_values_count(), 1)
      self.assertEqual(tree.get_height(), 1)

      tree.add(-2, 2, '2')
      self.assertEqual(tree.root.color, IntervalTreeNode.BLACK) 
      self.assertEqual(tree.root.max_end, 2)
      self.assertEqual(tree.get_intervals_count(), 1)
      self.assertEqual(tree.get_values_count(), 2)
      self.assertEqual(tree.get_height(), 1)
      
      tree.add(-3, 1, '3')
      self.assertEqual(tree.root.color, IntervalTreeNode.BLACK) 
      self.assertEqual(tree.root.max_end, 2) 
      self.assertEqual(tree.get_intervals_count(), 2)
      self.assertEqual(tree.get_values_count(), 3)
      self.assertEqual(tree.get_height(), 2)

      tree.add(-3, 1, '4')
      self.assertEqual(tree.root.color, IntervalTreeNode.BLACK) 
      self.assertEqual(tree.root.start, -2) 
      self.assertEqual(tree.root.end, 2) 
      self.assertEqual(tree.root.max_end, 2) 
      self.assertEqual(tree.get_intervals_count(), 2)
      self.assertEqual(tree.get_values_count(), 4)
      self.assertEqual(tree.get_height(), 2)


   def test_point_intersect(self):
      tree = IntervalTree()
      intervals = [(0, 5, '1'), (-2, 1, '2'), (-1, 20, '3'), (4, 6, '4'), (1, 3, '5'), (3, 10, '6')]
      expected_result = set(['1', '2', '3', '5'])
      for x in intervals:
         tree.add(x[0], x[1], x[2])
      result = set([i.values[0] for i in tree.point_intersect(1)])
      self.assertSetEqual(result, expected_result)
      

   def test_interval_overlap(self):
      tree = IntervalTree()
      intervals = [(0, 5, '1'), (-2, 1, '2'), (-1, 20, '3'), (4, 6, '4'), (1, 3, '5'), (3, 10, '6')]
      expected_result = set(['1', '3', '4', '6'])
      for x in intervals:
         tree.add(x[0], x[1], x[2])
      result = set([i.values[0] for i in tree.interval_overlap(5, 6)])
      self.assertSetEqual(result, expected_result)

 
   def test_nearest_left(self):
      tree = IntervalTree()
      intervals = [(0, 5, '1'), (-2, 1, '2'), (-1, 20, '3'), (4, 6, '4'), (1, 3, '5'), (3, 10, '6')]
      for x in intervals:
         tree.add(x[0], x[1], x[2])
      self.assertIsNone(tree.nearest_left(-10))
      self.assertEqual(tree.nearest_left(0).values[0], '3')
      self.assertEqual(tree.nearest_left(1).values[0], '1')
      self.assertEqual(tree.nearest_left(10).values[0], '4')


   def test_nearest_right(self):
      tree = IntervalTree()
      intervals = [(0, 5, '1'), (-2, 1, '2'), (-1, 20, '3'), (4, 6, '4'), (1, 3, '5'), (3, 10, '6')]
      for x in intervals:
         tree.add(x[0], x[1], x[2])
      self.assertIsNone(tree.nearest_right(10))
      self.assertEqual(tree.nearest_right(0).values[0], '1')
      self.assertEqual(tree.nearest_right(2).values[0], '6')
      self.assertEqual(tree.nearest_right(-3).values[0], '2')


   def test_k_first(self):
      tree = IntervalTree()
      intervals = [(0, 5, '1'), (-2, 1, '2'), (-1, 20, '3'), (4, 6, '4'), (1, 3, '5'), (3, 10, '6')]
      for x in intervals:
         tree.add(x[0], x[1], x[2])

      result = [i.values[0] for i in tree.k_first(0)]
      self.assertListEqual(result, [])
      
      result = [i.values[0] for i in tree.k_first(1)]
      self.assertListEqual(result, ['2'])

      result = [i.values[0] for i in tree.k_first(3)]
      self.assertListEqual(result, ['2', '3', '1'])
      
      result = [i.values[0] for i in tree.k_first(20)]
      self.assertListEqual(result, ['2', '3', '1', '5', '6', '4'])
 

   def test_k_last(self):
      tree = IntervalTree()
      intervals = [(0, 5, '1'), (-2, 1, '2'), (-1, 20, '3'), (4, 6, '4'), (1, 3, '5'), (3, 10, '6')]
      for x in intervals:
         tree.add(x[0], x[1], x[2])

      result = [i.values[0] for i in tree.k_last(0)]
      self.assertListEqual(result, [])
      
      result = [i.values[0] for i in tree.k_last(1)]
      self.assertListEqual(result, ['4'])

      result = [i.values[0] for i in tree.k_last(3)]
      self.assertListEqual(result, ['4', '6', '5'])
      
      result = [i.values[0] for i in tree.k_last(20)]
      self.assertListEqual(result, ['4', '6', '5', '1', '3', '2'])
 

   def test_k_nearest_left(self):
      tree = IntervalTree()
      intervals = [(0, 5, '1'), (-2, 1, '2'), (-1, 20, '3'), (4, 6, '4'), (1, 3, '5'), (3, 10, '6')]
      for x in intervals:
         tree.add(x[0], x[1], x[2])

      result = [i.values[0] for i in tree.k_nearest_left(-20, 20)]
      self.assertListEqual(result, [])

      result = [i.values[0] for i in tree.k_nearest_left(3, 1)]
      self.assertListEqual(result, ['5'])

      result = [i.values[0] for i in tree.k_nearest_left(3, 2)]
      self.assertListEqual(result, ['5', '1'])

      result = [i.values[0] for i in tree.k_nearest_left(3, 3)]
      self.assertListEqual(result, ['5', '1', '3'])

      result = [i.values[0] for i in tree.k_nearest_left(20, 20)]
      self.assertListEqual(result, ['4', '6', '5', '1', '3', '2'])


   def test_k_nearest_right(self):
      tree = IntervalTree()
      intervals = [(0, 5, '1'), (-2, 1, '2'), (-1, 20, '3'), (4, 6, '4'), (1, 3, '5'), (3, 10, '6')]
      for x in intervals:
         tree.add(x[0], x[1], x[2])

      result = [i.values[0] for i in tree.k_nearest_right(20, 20)]
      self.assertListEqual(result, [])

      result = [i.values[0] for i in tree.k_nearest_right(0, 1)]
      self.assertListEqual(result, ['1'])

      result = [i.values[0] for i in tree.k_nearest_right(0, 2)]
      self.assertListEqual(result, ['1', '5'])

      result = [i.values[0] for i in tree.k_nearest_right(0, 3)]
      self.assertListEqual(result, ['1', '5', '6'])

      result = [i.values[0] for i in tree.k_nearest_right(-20, 20)]
      self.assertListEqual(result, ['2', '3', '1', '5', '6', '4'])



if __name__ == '__main__':
   suite = unittest.TestSuite()
   test_loader = unittest.TestLoader()
   suite.addTest(test_loader.loadTestsFromTestCase(TestIntervalTreeNode))
   suite.addTest(test_loader.loadTestsFromTestCase(TestIntervalTree))
   unittest.TextTestRunner(verbosity = 2).run(suite)
