""" author : shakthi visagan (shak360), adam weiner (adamcweiner)
some edits: Farnaz Mohammadi
description: a file to hold the cell class
"""
import math
import numpy as np
from scipy import optimize, stats as sp


class CellNode:
    """Each cell in our tree will consist of a node containing these traits.

    This class includes many functions that assign the cell its properties, i.e.,
    the cell's generation, lifetime, true_state, etc."""

    def __init__(self, gen=1, linID=0, startT=0, endT=float('nan'), fate=None, left=None, right=None, parent=None, true_state=None, g1=float('nan'), g2=float('nan')):
        """
        Args:
        -----
        gen (int): the generation of the cell, root cells are of generation 1,
            each division adds 1 to the previous generation.

        linID (int): the lineage identity of the cell, keeps track of what
        lineage a cell belongs to.

        startT (float): the starting time of the cell, the point at which
        it spawned into existence.

        endT (float): the end time of the cell, the point at which it either
        divided or died, can be NaN.

        tau (float): [avoiding self.t, since that is a common function
        (i.e. transposing matrices)] tau is how long the cell lived.

        fate (0/1): the fate at the endT of the cell, 0 is death, 1 is division.

        left (obj): the left daughter of the cell, either returns a CellNode
        or NoneType object.

        right (obj): the right daughter of the cell, either returns a CellNode
        or NoneType object.

        parent (obj): the parent of the cell, returns a CellNode object
        (except at the root node)


        true_state (0/1): indicates whether cell is PC9 (0) or H1299 (1)

        fateObserved (T/F): marks whether the cell reached the true end
        of its lifetime (has truely died or divided)

        """
        self.gen = gen
        self.linID = linID
        self.startT = startT
        self.endT = endT
        self.tau = self.endT - self.startT
        self.fate = fate
        self.left = left
        self.right = right
        self.parent = parent

        self.true_state = true_state
        self.fateObserved = False
        self.g1 = g1
        self.g2 = g2

    def isParent(self):
        """ Returns true if the cell has at least one daughter; i.e., if either of the left or right daughter cells exist, it returns True. """
        return self.left or self.right

    def isChild(self):
        """ Returns true if this cell has a known parent. """
        return self.parent.isParent()

    def isRootParent(self):
        """ Returns whether this is a starting cell with no parent. """
        bool_parent = False
        if self.parent is None:
            assert self.gen == 1
            bool_parent = True
        return bool_parent

    def isLeaf(self):
        """ Returns True when a cell is a leaf with no children. """
        return self.left is None and self.right is None

    def calcTau(self):
        """
        Find the cell's lifetime by subtracting its endTime from startTime
        it makes sure the cell's lifetime is not NaN.
        """
        self.tau = self.endT - self.startT   # calculate tau here
        assert np.isfinite(self.tau), "Warning: your cell lifetime, {}, is a nan".format(self.tau)

    def isUnfinished(self):
        """ See if the cell is living or has already died/divided. """
        return math.isnan(self.endT) and self.fate is None   # returns true when cell is still alive

    def setUnfinished(self):
        """ Set a finished cell back to being unfinished. """
        self.endT = float('nan')
        self.fate = None
        self.tau = float('nan')

    def start_G2(self):
        """
        Returns the start time point of cell's G2 phase by adding the start time and the duration of G1.
        """
        return self.startT + self.g1

    def die(self, endT):
        """
        Cell dies without dividing.

        If the cell dies, the endTime is reached so we calculate the lifetime (tau)
        and the death is observed.
        Args:
        -----
            endT (float): end time of the cell.

        This function doesn't return
        """
        self.fate = False   # no division
        self.endT = endT    # mark endT
        self.calcTau()      # calculate Tau when cell dies
        self.fateObserved = True  # this cell has truly died

    def divide(self, endT):
        """
        Cell life ends through division. The two optional trackID arguments
        represents the trackIDs given to the two daughter cells.

        Args:
        -----
        endT (float): end time of the cell

        Returns:
        --------
        self.left (obj): left new born daughter cell
        self.right (obj): right new born daughter cell

        """
        self.endT = endT
        self.fate = True    # division
        self.calcTau()      # calculate Tau when cell dies
        self.fateObserved = True  # this cell has truly divided

        if self.isRootParent():

            self.left = CellNode(gen=self.gen + 1, linID=self.linID, startT=endT, parent=self, true_state=self.true_state)
            self.right = CellNode(gen=self.gen + 1, linID=self.linID, startT=endT, parent=self, true_state=self.true_state)
        else:
            self.left = CellNode(gen=self.gen + 1, linID=self.linID, startT=endT, parent=self, true_state=self.true_state)
            self.right = CellNode(gen=self.gen + 1, linID=self.linID, startT=endT, parent=self, true_state=self.true_state)

        return (self.left, self.right)

    def get_root_cell(self):
        """
        Gets the root cell associated with the cell.

        it keeps going up to the root by jumping from daughter cells to their
        mother cells by keeping track of their lineage ID (linID) until it reaches
        the first generation and then returns the last hold cell which is the
        ancestor of the lineage.

        Returns:
        --------
        curr_cell (obj): the ancestor cell if a lineage

        """
        cell_linID = self.linID
        curr_cell = self
        while curr_cell.gen > 1:
            curr_cell = curr_cell.parent
            assert cell_linID == curr_cell.linID
        assert cell_linID == curr_cell.linID
        return curr_cell
