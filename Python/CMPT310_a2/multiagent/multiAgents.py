# multiAgents.py
# --------------
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

# I MUST ACKNOWLEDGE I USED SOME AI TOOLS AND ONLINE WEBSITES TO DEBUG AND CLEANUP MY CODE SNIPPETS - ARSHIA HAMEDI

from util import manhattanDistance
from game import Directions
import random, util

from game import Agent
from pacman import GameState

class ReflexAgent(Agent):
    """
    A reflex agent chooses an action at each choice point by examining
    its alternatives via a state evaluation function.

    The code below is provided as a guide.  You are welcome to change
    it in any way you see fit, so long as you don't touch our method
    headers.
    """


    def getAction(self, gameState: GameState):
        """
        You do not need to change this method, but you're welcome to.

        getAction chooses among the best options according to the evaluation function.

        Just like in the previous project, getAction takes a GameState and returns
        some Directions.X for some X in the set {NORTH, SOUTH, WEST, EAST, STOP}
        """
        # Collect legal moves and successor states
        legalMoves = gameState.getLegalActions()

        # Choose one of the best actions
        scores = [self.evaluationFunction(gameState, action) for action in legalMoves]
        bestScore = max(scores)
        bestIndices = [index for index in range(len(scores)) if scores[index] == bestScore]
        chosenIndex = random.choice(bestIndices) # Pick randomly among the best

        "Add more of your code here if you want to"

        return legalMoves[chosenIndex]

    def evaluationFunction(self, currentGameState: GameState, action):
        """
        Design a better evaluation function here.

        The evaluation function takes in the current and proposed successor
        GameStates (pacman.py) and returns a number, where higher numbers are better.

        The code below extracts some useful information from the state, like the
        remaining food (newFood) and Pacman position after moving (newPos).
        newScaredTimes holds the number of moves that each ghost will remain
        scared because of Pacman having eaten a power pellet.

        Print out these variables to see what you're getting, then combine them
        to create a masterful evaluation function.
        """
        # Useful information you can extract from a GameState (pacman.py)
        successorGameState = currentGameState.generatePacmanSuccessor(action)
        newPos = successorGameState.getPacmanPosition()
        newFood = successorGameState.getFood()
        newGhostStates = successorGameState.getGhostStates()
        newScaredTimes = [ghostState.scaredTimer for ghostState in newGhostStates]

        "*** YOUR CODE HERE ***"
        # return successorGameState.getScore()

        foodDist = [manhattanDistance(newPos, food) for food in newFood.asList()]
        closestFoodDist = min(foodDist) if foodDist else 1
        ghostDist = [manhattanDistance(newPos, ghost.getPosition()) for ghost in newGhostStates]
        ghostPenalty = sum([-200 / (dist + 1) for dist, scared in zip(ghostDist, newScaredTimes) if dist < 2 and scared == 0])
        scaredGhostBonus = sum([200 / (dist + 1) for dist, scared in zip(ghostDist, newScaredTimes) if scared > 0])
        foodPenaltyLeft = -len(newFood.asList()) * 10
        score = successorGameState.getScore()
        score += 10 / closestFoodDist  # Reward for being closer to food
        score += ghostPenalty
        score += scaredGhostBonus
        score += foodPenaltyLeft

        return score
    
def scoreEvaluationFunction(currentGameState: GameState):
    """
    This default evaluation function just returns the score of the state.
    The score is the same one displayed in the Pacman GUI.

    This evaluation function is meant for use with adversarial search agents
    (not reflex agents).
    """
    return currentGameState.getScore()

class MultiAgentSearchAgent(Agent):
    """
    This class provides some common elements to all of your
    multi-agent searchers.  Any methods defined here will be available
    to the MinimaxPacmanAgent, AlphaBetaPacmanAgent & ExpectimaxPacmanAgent.

    You *do not* need to make any changes here, but you can if you want to
    add functionality to all your adversarial search agents.  Please do not
    remove anything, however.

    Note: this is an abstract class: one that should not be instantiated.  It's
    only partially specified, and designed to be extended.  Agent (game.py)
    is another abstract class.
    """

    def __init__(self, evalFn = 'scoreEvaluationFunction', depth = '2'):
        self.index = 0 # Pacman is always agent index 0
        self.evaluationFunction = util.lookup(evalFn, globals())
        self.depth = int(depth)

class MinimaxAgent(MultiAgentSearchAgent):
    """
    Your minimax agent (question 2)
    """

    def getAction(self, gameState: GameState):
        """
        Returns the minimax action from the current gameState using self.depth
        and self.evaluationFunction.

        Here are some method calls that might be useful when implementing minimax.

        gameState.getLegalActions(agentIndex):
        Returns a list of legal actions for an agent
        agentIndex=0 means Pacman, ghosts are >= 1

        gameState.generateSuccessor(agentIndex, action):
        Returns the successor game state after an agent takes an action

        gameState.getNumAgents():
        Returns the total number of agents in the game

        gameState.isWin():
        Returns whether or not the game state is a winning state

        gameState.isLose():
        Returns whether or not the game state is a losing state
        """
        "*** YOUR CODE HERE ***"
        # util.raiseNotDefined()
        # Start the minimax process and determine the best action for Pacman
        action, _ = self.minimax(gameState, 0, 0)
        return action

    def minimax(self, state, currentDepth, agentIndex):
        """
        Recursively calculate the minimax value and return the best action and its value.
        """
        numAgents = state.getNumAgents()
        
        # Check for terminal states or if the maximum depth is reached
        if state.isWin() or state.isLose() or currentDepth == self.depth * numAgents:
            return None, self.evaluationFunction(state)

        # Determine the next agent and depth
        nextAgent = (agentIndex + 1) % numAgents
        nextDepth = currentDepth + 1
        
        # Generate possible actions
        actions = state.getLegalActions(agentIndex)
        if not actions:
            return None, self.evaluationFunction(state)

        # Evaluate all possible successor states and calculate their minimax values
        results = [(action, self.minimax(state.generateSuccessor(agentIndex, action), nextDepth, nextAgent)[1])
                   for action in actions]

        # If the current agent is Pacman, maximize; if it is a ghost, minimize
        if agentIndex == 0:
            # Pacman's turn (maximizing player)
            return max(results, key=lambda x: x[1])
        else:
            # Ghosts' turn (minimizing player)
            return min(results, key=lambda x: x[1])



class AlphaBetaAgent(MultiAgentSearchAgent):
    """
    Your minimax agent with alpha-beta pruning (question 3)
    """

    def getAction(self, gameState: GameState):
        """
        Returns the minimax action using self.depth and self.evaluationFunction
        """
        "*** YOUR CODE HERE ***"
        # util.raiseNotDefined()
        def alphaBeta(agentIndex, depth, state, alpha, beta):
            if depth == self.depth or state.isWin() or state.isLose():
                return self.evaluationFunction(state)

            if agentIndex == 0:  # Pacman (Maximizer)
                value = float('-inf')
                for action in state.getLegalActions(agentIndex):
                    successor = state.generateSuccessor(agentIndex, action)
                    value = max(value, alphaBeta(1, depth, successor, alpha, beta))
                    if value > beta:
                        return value
                    alpha = max(alpha, value)
                return value
            else:  # Ghosts (Minimizer)
                value = float('inf')
                nextAgent = (agentIndex + 1) % state.getNumAgents()
                nextDepth = depth + 1 if nextAgent == 0 else depth
                for action in state.getLegalActions(agentIndex):
                    successor = state.generateSuccessor(agentIndex, action)
                    value = min(value, alphaBeta(nextAgent, nextDepth, successor, alpha, beta))
                    if value < alpha:
                        return value
                    beta = min(beta, value)
                return value

        actions = gameState.getLegalActions(0)
        bestAction = None
        bestValue = float('-inf')
        alpha = float('-inf')
        beta = float('inf')

        for action in actions:
            successor = gameState.generateSuccessor(0, action)
            value = alphaBeta(1, 0, successor, alpha, beta)
            if value > bestValue:
                bestValue = value
                bestAction = action
            alpha = max(alpha, bestValue)

        return bestAction


class ExpectimaxAgent(MultiAgentSearchAgent):
    """
      Your expectimax agent (question 4)
    """

    def getAction(self, gameState: GameState):
        """
        Returns the expectimax action using self.depth and self.evaluationFunction

        All ghosts should be modeled as choosing uniformly at random from their
        legal moves.
        """
        "*** YOUR CODE HERE ***"
        # util.raiseNotDefined()
        # Initiate the expectimax decision process to determine the best action for Pacman
        action = max(gameState.getLegalActions(0), key=lambda a: self.expectiMax(1, 0, gameState.generateSuccessor(0, a)))
        return action

    def expectiMax(self, agentIndex, depth, state):
        # Recursively calculate the expectimax value and return the expected value of the state.
        # Base case: check if the depth limit is reached or the game state is a terminal state
        if state.isWin() or state.isLose() or depth == self.depth:
            return self.evaluationFunction(state)
        actions = state.getLegalActions(agentIndex)
        numberOfAgents = state.getNumAgents()
        nextDepth = depth if (agentIndex + 1) % numberOfAgents else depth + 1
    
        if agentIndex == 0:  # Pacman's turn
            values = (self.expectiMax((agentIndex + 1) % numberOfAgents, nextDepth, state.generateSuccessor(agentIndex, a)) for a in actions)
            return max(values)
        else:  # Ghost's turn
            total = 0
            for action in actions:
                total += self.expectiMax((agentIndex + 1) % numberOfAgents, nextDepth, state.generateSuccessor(agentIndex, action))
            return total / len(actions) if actions else self.evaluationFunction(state)


def betterEvaluationFunction(currentGameState: GameState):
    """
    Your extreme ghost-hunting, pellet-nabbing, food-gobbling, unstoppable
    evaluation function (question 5).

    DESCRIPTION: <write something here so we know what you did>
    This evaluation function considers the following factors:
    - The reciprocal of the Manhattan distance to the nearest food pellet.
    - A penalty for proximity to active ghosts (not scared).
    - A bonus for proximity to scared ghosts (to eat them for points).
    - A penalty for the remaining number of food pellets to encourage faster wins.
    - A penalty for uneaten power pellets to encourage safer gameplay.
    - Incorporates the game's score to prioritize winning states.
    """
    "*** YOUR CODE HERE ***"
    # util.raiseNotDefined()
    pacmanPos = currentGameState.getPacmanPosition()
    foodList = currentGameState.getFood().asList()
    ghostStates = currentGameState.getGhostStates()
    scaredTimes = [ghostState.scaredTimer for ghostState in ghostStates]
    foodDistances = [manhattanDistance(pacmanPos, food) for food in foodList] # Distance to the nearest food
    closestFoodDistance = min(foodDistances) if foodDistances else 1
    foodScore = 10 / closestFoodDistance  
    ghostDistances = [manhattanDistance(pacmanPos, ghost.getPosition()) for ghost in ghostStates] # Ghost proximity penalty
    activeGhostPenalty = sum([-300 / (dist + 1) for dist, scared in zip(ghostDistances, scaredTimes) if scared == 0 and dist < 3])
    scaredGhostBonus = sum([250 / (dist + 1) for dist, scared in zip(ghostDistances, scaredTimes) if scared > 0]) # Scared ghost bonus
    remainingFoodPenalty = -len(foodList) * 15 # Remaining food penalty
    powerPellets = currentGameState.getCapsules() # Power pellet penalty (encourage eating power pellets)
    powerPelletPenalty = -len(powerPellets) * 20
    score = currentGameState.getScore()
    totalScore = score + foodScore + activeGhostPenalty + scaredGhostBonus + remainingFoodPenalty + powerPelletPenalty # Combine factors with the game score

    return totalScore

# Abbreviation
better = betterEvaluationFunction
