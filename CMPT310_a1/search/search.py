# search.py
# ---------
# Licensing Information:  You are free to use or extend these projects for
# educational purposes provided that (1) you do not distribute or publish
# solutions, (2) you retain this notice, and (3) you provide clear
# attribution to UC Berkeley, including a link to http://ai.berkeley.edu.
# 
# Attribution Information: The Pacman AI projects were developed at UC Berkeley.
# The core projects and autograders were primarily created by John DeNero
# (denero@cs.berkeley.edu) and Dan Klein (klein@cs.berkeley.edu).
# Student side autograding was added by Brad Miller, Nick Hay, and
# Pieter Abbeel (pabbeel@cs.berkeley.edu).

# I MUST ACKNOWLEDGE I USED SOME AI TOOLS TO DEBUG AND CLEANUP MY CODE SNIPPETS - ARSHIA HAMEDI
"""
In search.py, you will implement generic search algorithms which are called by
Pacman agents (in searchAgents.py).
"""

import util
from game import Directions
from typing import List

class SearchProblem:
    """
    This class outlines the structure of a search problem, but doesn't implement
    any of the methods (in object-oriented terminology: an abstract class).

    You do not need to change anything in this class, ever.
    """

    def getStartState(self):
        """
        Returns the start state for the search problem.
        """
        util.raiseNotDefined()

    def isGoalState(self, state):
        """
          state: Search state

        Returns True if and only if the state is a valid goal state.
        """
        util.raiseNotDefined()

    def getSuccessors(self, state):
        """
          state: Search state

        For a given state, this should return a list of triples, (successor,
        action, stepCost), where 'successor' is a successor to the current
        state, 'action' is the action required to get there, and 'stepCost' is
        the incremental cost of expanding to that successor.
        """
        util.raiseNotDefined()

    def getCostOfActions(self, actions):
        """
         actions: A list of actions to take

        This method returns the total cost of a particular sequence of actions.
        The sequence must be composed of legal moves.
        """
        util.raiseNotDefined()




def tinyMazeSearch(problem:SearchProblem) -> List[Directions]:
    """
    Returns a sequence of moves that solves tinyMaze.  For any other maze, the
    sequence of moves will be incorrect, so only use this for tinyMaze.
    """
    s = Directions.SOUTH
    w = Directions.WEST
    return  [s, s, w, s, w, w, s, w]

def depthFirstSearch(problem):
    """
    Search the deepest nodes in the search tree first.

    Your search algorithm needs to return a list of actions that reaches the
    goal. Make sure to implement a graph search algorithm.

    To get started, you might want to try some of these simple commands to
    understand the search problem that is being passed in:

    print("Start:", problem.getStartState())
    print("Is the start a goal?", problem.isGoalState(problem.getStartState()))
    print("Start's successors:", problem.getSuccessors(problem.getStartState()))
    """
    "*** YOUR CODE HERE ***"
    # util.raiseNotDefined()
    frontier = util.Stack()
    explored = set()

    start_state = problem.getStartState()
    frontier.push((start_state, []))

    while not frontier.isEmpty():
        current_state, path = frontier.pop()
        if current_state not in explored:
            explored.add(current_state)
            if problem.isGoalState(current_state):
                return path
            for next_state, action, _ in problem.getSuccessors(current_state):
                if next_state not in explored:
                    new_path = path + [action]
                    frontier.push((next_state, new_path))
    return []

def breadthFirstSearch(problem: SearchProblem) -> List[Directions]:
    """Search the shallowest nodes in the search tree first."""
    "*** YOUR CODE HERE ***"

    frontier = util.Queue()
    explored = set()
    start_state = problem.getStartState()
    frontier.push((start_state, []))

    while not frontier.isEmpty():
        state, actions = frontier.pop()
        if state in explored:
            continue
        explored.add(state)
        if problem.isGoalState(state):
            return actions
        for next_state, action, _ in problem.getSuccessors(state):
            if next_state not in explored:
                new_actions = actions + [action]
                frontier.push((next_state, new_actions))
    return []
    # util.raiseNotDefined()

def uniformCostSearch(problem: SearchProblem) -> List[Directions]:
    """Search the node of least total cost first."""
    "*** YOUR CODE HERE ***"
    frontier = util.PriorityQueue()
    initial_state = problem.getStartState()
    frontier.push((initial_state, []), 0)  
    cost_so_far = {initial_state: 0}

    while not frontier.isEmpty():
        current_state, path = frontier.pop()
        if problem.isGoalState(current_state):
            return path
        for successor, action, step_cost in problem.getSuccessors(current_state):
            new_cost = cost_so_far[current_state] + step_cost
            if successor not in cost_so_far or new_cost < cost_so_far[successor]:
                cost_so_far[successor] = new_cost
                new_path = path + [action]
                frontier.push((successor, new_path), new_cost)
    return []
    # util.raiseNotDefined()

def nullHeuristic(state, problem=None) -> float:
    """
    A heuristic function estimates the cost from the current state to the nearest
    goal in the provided SearchProblem.  This heuristic is trivial.
    """
    return 0

def aStarSearch(problem: SearchProblem, heuristic=nullHeuristic) -> List[Directions]:
    """Search the node that has the lowest combined cost and heuristic first."""
    "*** YOUR CODE HERE ***"
    from util import PriorityQueue

    frontier = PriorityQueue()
    startState = problem.getStartState()
    frontier.push((startState, [], 0), 0)  
    visited = {}

    while not frontier.isEmpty():
        currentState, actions, g = frontier.pop()

        if currentState in visited and visited[currentState] <= g:
            continue
        visited[currentState] = g

        if problem.isGoalState(currentState):
            return actions

        for successor, action, stepCost in problem.getSuccessors(currentState):
            newActions = actions + [action]
            newCost = g + stepCost
            f = newCost + heuristic(successor, problem)  
            frontier.push((successor, newActions, newCost), f)
    return []
    # util.raiseNotDefined()

# Abbreviations
bfs = breadthFirstSearch
dfs = depthFirstSearch
astar = aStarSearch
ucs = uniformCostSearch