"""A tree structure for efficient interval-based queries.

This module implements the interval tree datastucture for time effecient (log-complexity) querying of interval-based data.
The implementation is based on red-black binary tree.

.. moduleauthor:: Daniel Taliun <dtaliun@umich.edu>

"""

class IntervalTreeNode:
   """This class represents an interval in the interval tree.

   Every interval has start and end positions, and (optional) associated value of any type.

   Attributes:
      start (long): start position.
      end (long): end position.
      values (list): list of associated values (of any type).

   """

   _BLACK = 0
   _RED = 1

   def __init__(self, start, end, max_end, value, deviation, color, left, right):
      self.start = start
      self.end = end
      self.max_end = max_end
      if value is not None:
         self.values = [value]
      else:
         self.values = []
      self.deviation = deviation
      self.color = color
      self.parent = None

      self.left = left
      if self.left is not None:
         left.parent = self

      self.right = right
      if self.right is not None:
          right.parent = self


   def compare(self, interval):
      if self.start == interval.start:
         if self.end == interval.end:
            return 0
         elif self.end < interval.end:
            return -1
      elif self.start < interval.start:
         return -1
      return 1


   def get_sibling(self):
      if self.parent is not None:
         if self is self.parent.left:
            return self.parent.right
         else:
            return self.parent.left
      return None


   def get_uncle(self):
      if self.parent is not None:
         return self.parent.get_sibling()
      return None


   def get_grandparent(self):
      if self.parent is not None:
         return self.parent.parent
      return None


   def get_height(self):
      heights = [1]
      for node in [self.left, self.right]:
         if node is not None:
            heights.append(1 + node.get_height())
      return max(heights)

   def get_values_count(self):
      count = len(self.values)
      for node in [self.left, self.right]:
         if node is not None:
            count += node.get_values_count()
      return count

   def get_intervals_count(self):
      count = 1
      for node in [self.left, self.right]:
         if node is not None:
            count += node.get_intervals_count()
      return count


class IntervalTree:
   """This class implements the interval tree based on red-black binary tree.
   """
   def __init__(self):
      self.root = None


   def add(self, start, end, value = None):
      """Adds interval to the interval tree.

      Note:
         No check if start < end is done.

      Args:
         start (long): interval's start position.
         end (long): interval's end position.
         value: (optional) associated value of any type.

      """
      self.__add(start, end, value, 0)


   def __add(self, start, end, value, deviation):
      new_node = IntervalTreeNode(start, end, end, value, deviation, IntervalTreeNode._RED, None, None)

      if self.root is None:
         self.root = new_node
      else:
         parent = self.root
         comparison = 0
         while  True:
            comparison = new_node.compare(parent)
            if comparison == 0:
               if value is not None:
                  parent.values.append(value)
               return
            elif comparison < 0:
               parent.max_end = max(parent.max_end, new_node.max_end)
               if parent.left is None:
                  parent.left = new_node
                  break
               else:
                  parent = parent.left
            else:
               parent.max_end = max(parent.max_end, new_node.max_end)
               if parent.right is None:
                  parent.right = new_node
                  break
               else:
                  parent = parent.right
         new_node.parent = parent

      self.__balance(new_node)


   def __balance(self, node):
      if node.parent is None:
         node.color = IntervalTreeNode._BLACK  # Root is always black
      elif node.parent.color != IntervalTreeNode._BLACK: # Parent of RED node must be black
         uncle = node.get_uncle()
         grandparent = node.get_grandparent() # at this point, grandparent exists always
         if uncle is not None and uncle.color == IntervalTreeNode._RED: # All children of red node must be black
            node.parent.color = IntervalTreeNode._BLACK
            uncle.color = IntervalTreeNode._BLACK
            grandparent.color = IntervalTreeNode._RED
            self.__balance(grandparent)
         else:
            # Unify subtree for the last step with appropriate rotations: parent and its child must be on the same side relatively to grandparent (either left or right).
            if node is node.parent.right and node.parent is grandparent.left:
               self.__rotate_left(node.parent)
               node = node.left
            elif node is node.parent.left and node.parent is grandparent.right:
               self.__rotate_right(node.parent)
               node = node.right
            # The last step:
            node.parent.color = IntervalTreeNode._BLACK
            grandparent.color = IntervalTreeNode._RED
            if node is node.parent.left and node.parent is grandparent.left:
               self.__rotate_right(grandparent)
            else:
               self.__rotate_left(grandparent)


   def __rotate_left(self, node):
      pivot_node = node.right
      node.right = pivot_node.left
      if pivot_node.left is not None:
         pivot_node.left.parent = node
      pivot_node.left = node
      pivot_node.parent = node.parent
      if node is self.root:
         self.root = pivot_node
      elif node is node.parent.left:
         node.parent.left = pivot_node
      else:
         node.parent.right = pivot_node
      node.parent = pivot_node
      self.__update_max_end(node)
      self.__update_max_end(pivot_node)


   def __rotate_right(self, node):
      pivot_node = node.left
      node.left = pivot_node.right
      if pivot_node.right is not None:
         pivot_node.right.parent = node
      pivot_node.right = node
      pivot_node.parent = node.parent
      if node is self.root:
         self.root = pivot_node
      elif node is node.parent.left:
         node.parent.left = pivot_node
      else:
         node.parent.right = pivot_node
      node.parent = pivot_node
      self.__update_max_end(node)
      self.__update_max_end(pivot_node)


   def __update_max_end(self, node):
      if node.left is None:
         if node.right is None:
            node.max_end = node.end
         else:
            node.max_end = max(node.end, node.right.max_end)
      else:
         if node.right is None:
            node.max_end = max(node.end, node.left.max_end)
         else:
            node.max_end = max(node.end, node.left.max_end, node.right.max_end)


   def point_intersect(self, position):
      """Finds all intervals that intersect the specified chromosomal position.

      Args:
         position (long): point position.

      Yields:
         IntervalTreeNode: intersecting interval.

      """
      if self.root is not None:
         nodes = [self.root]
         while nodes:
            node = nodes.pop()
            if node.left is not None and node.left.max_end >= position:
               nodes.append(node.left)
            if node.right is not None and node.right.max_end >= position and node.start <= position:
               nodes.append(node.right)
            if node.start <= position and position <= node.end:
               yield node

   def interval_overlap(self, start, end):
      """Finds all intervals that overlap specified interval.

      Args:
         start (long): start position of the interval.
         end (long): end position of the interval.

      Yields:
         IntervalTreeNode: overlapping interval.

      """
      if self.root is not None:
         nodes = [self.root]
         while nodes:
            node = nodes.pop()
            if node.left is not None and node.left.max_end >= start:
               nodes.append(node.left)
            if node.right is not None and node.right.max_end >= start and end >= node.start:
               nodes.append(node.right)
            if start <= node.end and end >= node.start:
               yield node


   def nearest_left(self, position):
      """Find the closest interval to the left from the specified position.

      The distance to the closest interval is measured with respect to its start position.

      Args:
         position (long): point position.

      Returns:
         IntervalTreeNode: interval.

      """
      if self.root is None:
         return None
      nodes = [self.root]
      nearest_node = None
      while nodes:
         node = nodes.pop()
         if node.start < position:
            nearest_node = node
            if node.right is not None:
               nodes.append(node.right)
         elif node.left is not None:
            nodes.append(node.left)
      return nearest_node


   def nearest_right(self, position):
      """Find the closest interval to the right from the specified position.

      The distance to the closest interval is measured with respect to its start position.

      Args:
         position (long): point position.

      Returns:
         IntervalTree: interval.

      """
      if self.root is None:
         return None
      nodes = [self.root]
      nearest_node = None
      while nodes:
         node = nodes.pop()
         if node.start >= position:
            nearest_node = node
            if node.left is not None:
               nodes.append(node.left)
         elif node.right is not None:
            nodes.append(node.right)
      return nearest_node


   def k_first(self, k):
      """Iterates over K leftmost intervals.

      Args:
         k (int): maximal number of leftmost intervals.

      Yields:
         IntervalTreeNode: interval.

      """
      nodes = []
      node = self.root
      while k > 0:
         while node is not None:
            nodes.append(node)
            node = node.left
         if nodes:
            node = nodes.pop()
            k -= 1
            yield node
            node = node.right
         else:
            break


   def __k_first(self, node, k):
      nodes = []
      while k > 0:
         while node is not None:
            nodes.append(node)
            node = node.left
         if nodes:
            node = nodes.pop()
            k -= 1
            yield node
            node = node.right
         else:
            break

   def ascending(self):
      """Iterates over all intervals in ascending order (i.e. from the leftmost to rightmost interval).

      Yields:
         IntervalTreeNode: interval.

      """
      nodes = []
      node = self.root
      while True:
         while node is not None:
            nodes.append(node)
            node = node.left
         if nodes:
            node = nodes.pop()
            yield node
            node = node.right
         else:
            break

   def k_last(self, k):
      """Iterates over K rightmost intervals.

      Args:
         k (int): maximal number of rightmost intervals.

      Yields:
         IntervalTreeNode: interval.

      """
      nodes = []
      node = self.root
      while k > 0:
         while node is not None:
            nodes.append(node)
            node = node.right
         if nodes:
            node = nodes.pop()
            k -= 1
            yield node
            node = node.left
         else:
            break

   def __k_last(self, node, k):
      nodes = []
      while k > 0:
         while node is not None:
            nodes.append(node)
            node = node.right
         if nodes:
            node = nodes.pop()
            k -= 1
            yield node
            node = node.left
         else:
            break

   def descending(self):
      """Iterates over all intervals in descending order (i.e. from the rightmost to the leftmost interval).

      Yields:
         IntervalTreeNode: interval.

      """
      nodes = []
      node = self.root
      while True:
         while node is not None:
            nodes.append(node)
            node = node.right
         if nodes:
            node = nodes.pop()
            yield node
            node = node.left
         else:
            break

   def k_nearest_left(self, position, k):
      """Iterates over K closest interval to the left from the specified position.

      The distance to the closest interval is measured with respect to its start position.

      Args:
         position (long): point position.
         k (int): maximal number of closest intervals.

      Yields:
         IntervalTreeNode: interval.

      """
      node = self.nearest_left(position)
      if node is None:
         return
      nodes = [node]
      while nodes:
         if k <= 0:
            break

         node = nodes.pop()
         k -= 1
         yield node

         if k > 0 and node.left is not None:
            k1 = k
            for x in self.__k_last(node.left, k1):
               k -= 1
               yield x

         if k > 0:
            while node.parent is not None and node.parent.right is not node:
               node = node.parent
            if node.parent is not None:
               nodes.append(node.parent)


   def k_nearest_right(self, position, k):
      """Iterates over K closest interval to the right from the specified position.

      The distance to the closest interval is measured with respect to its start position.

      Args:
         position (long): point position.
         k (int): maximal number of closest intervals.

      Yields:
         IntervalTreeNode: interval.
      """
      node = self.nearest_right(position)
      if node is None:
         return
      nodes = [node]
      while nodes:
         if k <= 0:
            break

         node = nodes.pop()
         k -= 1
         yield node

         if k > 0 and node.right is not None:
            k1 = k
            for x in self.__k_first(node.right, k1):
               k -= 1
               yield x

         if k > 0:
            while node.parent is not None and node.parent.left is not node:
               node = node.parent
            if node.parent is not None:
               nodes.append(node.parent)


   def merge(self):
      """Returns new interval tree with no overlapping intervals.

      Overlapping intervals are merged into a new single interval.
      The value of the new interval is the list of all values from the merged overlapping intervals.

      Returns:
         IntervalTree: new interval tree.

      """
      new_tree = IntervalTree()
      ascending = self.ascending()
      merged_start = None
      merged_end = None
      merged_values = None
      x = next(ascending, None)
      if x is not None:
         merged_start = x.start
         merged_end = x.end
         merged_values = x.values[:]
      for x in ascending:
         if merged_end < x.start:
            if len(merged_values) > 0:
               for v in merged_values:
                  new_tree.add(merged_start, merged_end, v)
            else:
               new_tree.add(merged_start, merged_end)
            merged_start = x.start
            merged_end = x.end
            merged_values = x.values[:]
         else:
            if merged_end < x.end:
               merged_end = x.end
            merged_values.extend(x.values)
      if merged_start is not None and merged_end is not None and merged_values is not None:
         if len(merged_values) > 0:
            for v in merged_values:
               new_tree.add(merged_start, merged_end, v)
         else:
             new_tree.add(merged_start, merged_end)
      return new_tree


   def complementary(self):
      """Returns new interval tree constructed from the current interval tree by turning gaps that are not covered by any interval into intervals in the new tree.

      Note:
         Intervals in the new interval tree do not store any values.

      Returns:
         IntervalTree: new interval tree.

      """
      new_tree = IntervalTree()
      start = None
      for x in self.ascending():
         if start is not None and start < x.start:
            new_tree.add(start, x.start - 1)
         start = x.end + 1
      return new_tree


   def get_values_count(self):
      """int: Total number of stored values within intervals in the interval tree."""
      return 0 if self.root is None else self.root.get_values_count()
   __len__ = get_values_count


   def get_intervals_count(self):
      """int: Total number of intervals in the interval tree."""
      return 0 if self.root is None else self.root.get_intervals_count()


   def get_height(self):
      """int: Number of levels in the interval tree based on the longest branch."""
      return 0 if self.root is None else self.root.get_height()


   def __str__(self):
      rv = []
      height = self.get_height()
      def print_node(node, depth):
         if node is not None:
            print_node(node.left, depth+1)
            rv.append('{:<{}}'.format(
               ' '*3*depth + ('blk' if node.color == IntervalTreeNode._BLACK else 'red'),
               height*3) +
                      ' [{start:2},{end:2}]({max_end:2}) : {values}'.format(**node.__dict__))
            print_node(node.right, depth+1)
      print_node(self.root, 0)
      return '\n'.join(rv)
