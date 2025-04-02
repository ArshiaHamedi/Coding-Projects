# valueIterationAgents.py
# -----------------------
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


# valueIterationAgents.py
# -----------------------
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


import mdp, util

from learningAgents import ValueEstimationAgent
import collections

class ValueIterationAgent(ValueEstimationAgent):
    """
        * Please read learningAgents.py before reading this.*

        A ValueIterationAgent takes a Markov decision process
        (see mdp.py) on initialization and runs value iteration
        for a given number of iterations using the supplied
        discount factor.
    """
    def __init__(self, mdp, discount = 0.9, iterations = 100):
        """
          Your value iteration agent should take an mdp on
          construction, run the indicated number of iterations
          and then act according to the resulting policy.

          Some useful mdp methods you will use:
              mdp.getStates()
              mdp.getPossibleActions(state)
              mdp.getTransitionStatesAndProbs(state, action)
              mdp.getReward(state, action, nextState)
              mdp.isTerminal(state)
        """
        self.mdp = mdp
        self.discount = discount
        self.iterations = iterations
        self.values = util.Counter() # A Counter is a dict with default 0
        self.runValueIteration()

    def runValueIteration(self):
        # Write value iteration code here
        "*** YOUR CODE HERE ***"
        for _ in range(self.iterations):
            new_values = self.values.copy()  # Copy current values for batch update

            for state in self.mdp.getStates():
                if self.mdp.isTerminal(state):
                    continue  # Keep terminal state values unchanged

                actions = self.mdp.getPossibleActions(state)
                if actions:
                    new_values[state] = max(
                        self.computeQValueFromValues(state, action) for action in actions
                    )

            self.values = new_values  # Update values after all states are processed
        
    def getValue(self, state):
        """
          Return the value of the state (computed in __init__).
        """
        return self.values[state]


    def computeQValueFromValues(self, state, action):
        """
          Compute the Q-value of action in state from the
          value function stored in self.values.
        """
        "*** YOUR CODE HERE ***"
        return sum(
            prob * (self.mdp.getReward(state, action, nextState) + self.discount * self.values[nextState])
            for nextState, prob in self.mdp.getTransitionStatesAndProbs(state, action)
        )
        # util.raiseNotDefined()

    def computeActionFromValues(self, state):
        """
          The policy is the best action in the given state
          according to the values currently stored in self.values.

          You may break ties any way you see fit.  Note that if
          there are no legal actions, which is the case at the
          terminal state, you should return None.
        """
        "*** YOUR CODE HERE ***"
        actions = self.mdp.getPossibleActions(state)
        if not actions:
            return None  # Return None for terminal states

        return max(actions, key=lambda action: self.computeQValueFromValues(state, action), default=None)
        # util.raiseNotDefined()

    def getPolicy(self, state):
        return self.computeActionFromValues(state)

    def getAction(self, state):
        "Returns the policy at the state (no exploration)."
        return self.computeActionFromValues(state)

    def getQValue(self, state, action):
        return self.computeQValueFromValues(state, action)

class AsynchronousValueIterationAgent(ValueIterationAgent):
    """
        * Please read learningAgents.py before reading this.*

        An AsynchronousValueIterationAgent takes a Markov decision process
        (see mdp.py) on initialization and runs cyclic value iteration
        for a given number of iterations using the supplied
        discount factor.
    """
    def __init__(self, mdp, discount = 0.9, iterations = 1000):
        """
          Your cyclic value iteration agent should take an mdp on
          construction, run the indicated number of iterations,
          and then act according to the resulting policy. Each iteration
          updates the value of only one state, which cycles through
          the states list. If the chosen state is terminal, nothing
          happens in that iteration.

          Some useful mdp methods you will use:
              mdp.getStates()
              mdp.getPossibleActions(state)
              mdp.getTransitionStatesAndProbs(state, action)
              mdp.getReward(state)
              mdp.isTerminal(state)
        """
        ValueIterationAgent.__init__(self, mdp, discount, iterations)

    def runValueIteration(self):
        "*** YOUR CODE HERE ***"
        states = self.mdp.getStates()
    
        for i in range(self.iterations):
            state = states[i % len(states)]  # Process states cyclically

            if self.mdp.isTerminal(state):
                continue  # Skip terminal states

            actions = self.mdp.getPossibleActions(state)
            if actions:
                self.values[state] = max(
                    self.computeQValueFromValues(state, action) for action in actions
            )

class PrioritizedSweepingValueIterationAgent(AsynchronousValueIterationAgent):
    """
        * Please read learningAgents.py before reading this.*

        A PrioritizedSweepingValueIterationAgent takes a Markov decision process
        (see mdp.py) on initialization and runs prioritized sweeping value iteration
        for a given number of iterations using the supplied parameters.
    """
    def __init__(self, mdp, discount = 0.9, iterations = 100, theta = 1e-5):
        """
          Your prioritized sweeping value iteration agent should take an mdp on
          construction, run the indicated number of iterations,
          and then act according to the resulting policy.
        """
        self.theta = theta
        ValueIterationAgent.__init__(self, mdp, discount, iterations)

    def runValueIteration(self):
        "*** YOUR CODE HERE ***"
        # Step 1: Compute predecessors for all states
        predecessors = {state: set() for state in self.mdp.getStates()}  # Use sets to avoid duplicates

        for state in self.mdp.getStates():
            for action in self.mdp.getPossibleActions(state):
                for next_state, prob in self.mdp.getTransitionStatesAndProbs(state, action):
                    if prob > 0:  # Only consider meaningful transitions
                        predecessors[next_state].add(state)

        # Step 2: Initialize the priority queue
        priority_queue = util.PriorityQueue()

        for state in self.mdp.getStates():
            if not self.mdp.isTerminal(state):  # Ignore terminal states
                max_q_value = max(
                    [self.computeQValueFromValues(state, action) for action in self.mdp.getPossibleActions(state)],
                    default=0
                )
                diff = abs(self.values[state] - max_q_value)
                priority_queue.push(state, -diff)  # Use negative diff since PriorityQueue is a min-heap

        # Step 3: Perform the prioritized sweeping value iteration
        for _ in range(self.iterations):
            if priority_queue.isEmpty():
                break  # Stop if no more updates are needed

            state = priority_queue.pop()

            if not self.mdp.isTerminal(state):  # Only update non-terminal states
                self.values[state] = max(
                    [self.computeQValueFromValues(state, action) for action in self.mdp.getPossibleActions(state)],
                    default=0
                )

            # Step 4: Update predecessors
            for p in predecessors[state]:
                if self.mdp.isTerminal(p):
                    continue  # Skip terminal predecessors

                max_q_value = max(
                    [self.computeQValueFromValues(p, action) for action in self.mdp.getPossibleActions(p)],
                    default=0
                )
                diff = abs(self.values[p] - max_q_value)

                if diff > self.theta:
                    priority_queue.update(p, -diff)  # Push/update only if significant change
