## Slide 1 - Bachelor's Project

- introduce myself + my project


## Slide 2 - Project Aim

- suggest TA as a solution to the JSP
- aim to solve other NP-complete problems like the TSP


## Slide 3 - The JSP

- introduce the JSP with some e.g.
    - talk through the constraints
        - **precedence rules within the same job**


## Slide 4 - The JSP

- show an optimal solution
    - emphasise that not unique
- solutions frequently offered from the field of Operations Research => want to offer another perspective on it: TA, which brings us to the next slide


## Slide 5 - TA for Tasks

- TA is similar to the automata seen in the Automata & Complexity course, with the addition of a state clock
- explain the components of TA on the e.g.:
    - TA for each tasks (name) and its corresponding machine
    - start
    - clocks named for the corresponding tasks,  show the time unit at transition time
    - state with bar = start state
    - state without bar = in execution state
    - final accepting state _f_


## Slide 6 - TA for Jobs

- combine the TA for tasks into TA for jobs
- explain the components of TA of jobs:
    - same as TA for tasks, except the automaton can have multiple machines
    - clocks named for the corresponding tasks


## Slide 7 - Global TA

- combine the TA for jobs into a global TA = the start states of each TA for jobs are combined into its start state
- supposed to emphasise the resource conflict - only one of the paths can be taken
- explain the components of the global TA:
    - start
    - "flow" through all the paths
    - goal is to reach _ff_ final state
- given the paths, we still need to compute the time needed to go through them - global automaton doesn't show the time


## Slide 8 - Immediate Runs Graph

- shows the time elapsed between states and the global time
- created from the transitions in the global TA
(- some state collapsed, because similar)
- explain the components of the graph:
    - perpendicular signs = not in use
    - 1st + 2nd entries <-> machine local time, 3rd entry = global time
    - walk through an e.g.
- shortest global time = optimal scheduling

- now that we've seen how it works for the JSP, time to show how it can be applied to the TSP


## Slide 9 - TSP

- as seen in courses such as Networks & Graphs, just a brief refresher:
    - number of cities to visit by a salesman
    - cost to go from one city to the other
    - modelled as:
        - complete graph/distance matrix
        - given a starting point, the goal is to tour all the cities and return to the start with the lowest cost possible
    - talk briefly about 1st 2 bullet points


## Slide 10 - TSP as JSP Modelling

- we now want to model the TSP as a JSP in order to solve it as one
- we have the recipe:
    - one salesman = one machine
    - _n_ cities = _n_ tasks
    - how to group the cities into tasks, given that there's no precedence rules?


## Slide 11 - Average-Link Clustering

- form clusters based on the shortest path to a group/city
- walk through e.g.:
    - B & D closest, so linked together
    - A closer to BD than C
- will use the same to cluster cities into jobs, random order, the grouping is important - cities with the same colours belong in the same job


## Slide 12 - TSP as JSP Modelling Example

- trivial e.g. for the simplicity of drawing automata
- from the image, A & B can be clustered together, C is separate
- as before, create the timed automata for everything

## Slide 13 - Solving the TSP using TA

- directly create global automaton, with the addition that the cost here depends sequentially on the previous city, so added on the transition in the global automaton
- clocks have the name of the cities
- graph also shows the currently visited city
- add the cost to return in the starting point at the end
- A - B - C - A shortest route

## Slide 14 - Conclusion

- show how we solved the TSP with TA, given that it can be modelled as a JSP
- this theoretical approach calls for an actual implementation and complexity analysis