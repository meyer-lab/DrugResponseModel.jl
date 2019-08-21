"""This file inclues functions to:
    1. Generating the Lineages of cells,
    2. Number of cells in G1/G2/total in every time point,
    3. Separates cells which have the same root parent in list of lists."""

from .cell import CellNode
import math
import numpy as np
from scipy import optimize, stats as sp

##----------------------------------- Create Lineage ------------------------------------##


def generateLineageWithTime(initCells, experimentTime, locBern, g1_a=None, g1_b=None, g2_a=None, g2_b=None):
    """
    generates a list of objects (cells) in a lineage.

    Given an experimental end time, a Bernoulli distribution for
    dividing/dying and two Gamma parameters for cell lifetime,
    it creates objects as cells and puts them in a list.


    Bernoulli distribution:
        It has one parameter (locBern)
        If locBern = 0.80 then 80% of the times the cells will divide and 20%
        of the times they die; meaning, for every cell that is generated,
        a random number (either 1:divide, or 0:die) is picked from the
        distribution with this parameter and the fate of the cell is assigned.


    Gamma distribution:
        It has two parameters(shape , scale) {here are a and b}
        Used as alifetime generator for cells. Here to generate the cells we
        specify the two parameters and it will return a number that we assign to cell's
        lifetime, as G1 phase and G2 phase of the cell cycle.

    Args:
    -----
        initCells (int): the number of initial cells to initiate the tree with
        experimentTime (int) [hours]: the time that the experiment will be running
        to allow for the cells to grow
        locBern (float): the Bernoulli distribution parameter
        (p = success) for fate assignment (either the cell dies or divides)
        range = [0, 1]
        g1_a, g2_a: shape parameters of Gamma for G1 and G2 phase of the cell cycle.
        g1_b, g2_b: scale parameters of Gamma for G1 and G2 phase of the cell cycle.

    Returns:
    --------
        lineage (list): A list of objects (cells) that creates the tree.
    """

    # create an empty lineage
    lineage = []

    # initialize the list with cells
    for ii in range(initCells):
        lineage.append(CellNode(startT=0, linID=ii))

    # have cell divide/die according to distribution
    for cell in lineage:
        cell.g1 = sp.gamma.rvs(g1_a, scale=g1_b)
        cell.g2 = sp.gamma.rvs(g2_a, scale=g2_b)

        if cell.isUnfinished():
            cell.tau = cell.g1 + cell.g2
            cell.endT = cell.startT + cell.tau

            if cell.endT < experimentTime:   # determine fate only if endT is within range
                # assign cell fate
                cell.fate = sp.bernoulli.rvs(locBern)  # assign fate

                # divide or die based on fate
                if cell.fate:
                    temp1, temp2 = cell.divide(cell.endT)  # cell divides
                    # append the children to the list
                    lineage.append(temp1)
                    lineage.append(temp2)
                else:
                    cell.die(cell.endT)

    return lineage

##------------------------- How many cells are in G1 or G2? --------------------------------##


def inG1_or_G2(X, time):
    """
    This function determines whether the cell is in G1 phase or in G2 phase.

    Args:
    -----
        X (list): is the lineage, a list of objects representing cells.
        time (list): a list -- could be np.linspace() -- including time points of
        duration of the time experiment is being conducted.

    Returns:
    --------
        num_G1 (list):  a list of # of cells in G1 at each time point
        num_G2 (list):  a list of # of cells in G2 at each time point
        num_cell (list):  a list of total # of cells at each time point
    """

    num_G1 = []
    num_G2 = []
    num_cell = []

    for t in time:
        count_G1 = 0
        count_G2 = 0
        count_numCell = 0

        for cell in X:
            g2 = cell.start_G2()

            # if the time point is between the cell's start time and the end of G1 phase, then count it as being in G1.
            if cell.startT <= t <= g2:
                count_G1 += 1

            # if the time point is between the start of the cell's G1 phase and the end time, then count it as being in G2.
            if g2 <= t <= cell.endT:
                count_G2 += 1

            # if the time point is within the cell's lifetime, count it as being alive.
            if cell.startT <= t <= cell.endT:
                count_numCell += 1

        num_G1.append(count_G1)
        num_G2.append(count_G2)
        num_cell.append(count_numCell)

    return num_G1, num_G2, num_cell


##--------------------- Separate different lineages by their root parent -------------------------##

def separate_pop(numLineages, X):
    """
    This function separates each lineage by their root parent, using their linID.

    Args:
    -----
        numLineages (int): the number of lineages, which here basically is the number of initial cells.
        X (list): is the lineage, a list of objects representing cells.

    Returns:
    --------
        population (list of lists): a list that holds lists of cells that belong tothe same parent.
    """

    population = []
    for i in range(numLineages):
        list_cell = []

        for cell in X:
            if cell.linID == i:
                list_cell.append(cell)
        population.append(list_cell)

    return population
