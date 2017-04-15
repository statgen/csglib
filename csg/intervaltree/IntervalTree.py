
class IntervalTreeNode:
   BLACK = 0
   RED = 1

   def __init__(self, start, end, max_end, value, deviation, color, left, right):
      self.start = start
      self.end = end
      self.max_end = max_end
      self.values = [value]
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
         if id(self) == id(self.parent.left):
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



class IntervalTree:

   def __init__(self):
      self.root = None


   def add(self, start, end, value):
      self.__add(start, end, value, 0)


   def __add(self, start, end, value, deviation):
      new_node = IntervalTreeNode(start, end, end, value, deviation, IntervalTreeNode.RED, None, None)

      if self.root is None:
         self.root = new_node
      else:
         parent = self.root
         comparison = 0
         while  True:
            comparison = new_node.compare(parent)
            if comparison == 0:
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
         node.color = IntervalTreeNode.BLACK  # Root is always black
      elif node.parent.color != IntervalTreeNode.BLACK: # Parent of RED node must be black
         uncle = node.get_uncle()
         grandparent = node.get_grandparent() # at this point, grandparent exists always
         if uncle is not None and uncle.color == IntervalTreeNode.RED: # All children of red node must be black
            node.parent.color = IntervalTreeNode.BLACK
            uncle.color = IntervalTreeNode.BLACK
            grandparent.color = IntervalTreeNode.RED
            self.__balance(grandparent)
         else:
            # Unify subtree for the last step with appropriate rotations: parent and its child must be on the same side relatively to grandparent (either left or right).
            if id(node) == id(node.parent.right) and id(node.parent) == id(grandparent.left):
               self.__rotate_left(node.parent)
               node = node.left
            elif id(node) == id(node.parent.left) and id(node.parent) == id(grandparent.right):
               self.__rotate_right(node.parent)
               node = node.right
            # The last step:
            node.parent.color = IntervalTreeNode.BLACK
            grandparent.color = IntervalTreeNode.RED
            if id(node) == id(node.parent.left) and id(node.parent) == id(grandparent.left):
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
      if id(node) == id(self.root):
         self.root = pivot_node
      elif id(node) == id(node.parent.left):
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
      if id(node) == id(self.root):
         self.root = pivot_node
      elif id(node) == id(node.parent.left):
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
            while node.parent is not None and id(node.parent.right) != id(node):
               node = node.parent
            if node.parent is not None:
               nodes.append(node.parent)


   def k_nearest_right(self, position, k):
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
            while node.parent is not None and id(node.parent.left) != id(node):
               node = node.parent
            if node.parent is not None:
               nodes.append(node.parent)



   def get_values_count(self):
      return self.__get_values_count(self.root, 0)

   def __get_values_count(self, node, count):
      if node is not None:
         count += len(node.values)
      else:
         return count
      count = self.__get_values_count(node.left, count)
      count = self.__get_values_count(node.right, count)
      return count


   def get_intervals_count(self):
      return self.__get_intervals_count(self.root, 0)

   def __get_intervals_count(self, node, count):
      if node is not None:
         count += 1
      else:
         return count
      count = self.__get_intervals_count(node.left, count)
      count = self.__get_intervals_count(node.right, count)
      return count


   def get_height(self):
      return self.__get_height(self.root, 0, 0)

   def __get_height(self, node, path_height, max_height):
      if node is not None:
         path_height += 1
      else:
         return max(path_height, max_height)
      max_height = self.__get_height(node.left, path_height, max_height)
      max_height = self.__get_height(node.right, path_height, max_height)
      return max_height
