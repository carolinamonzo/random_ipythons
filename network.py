#! usr/bin/env python2.7

# 2018-01-08; Carolina Monzo

import networkx as nx
import random
import matplotlib.pyplot as plt

def random_network(nodes):
    '''
    Function to create a random network with N nodes and average degree <k> = 4
    Input: Number of nodes
    Output: random graph, each node will have a degree of 4
    '''

    grafo = nx.random_regular_graph(4, nodes, seed = random.seed(a=long))

    return(grafo)

def efficiency(grafo, lambd):
    '''
    Function to calculate the efficiency of a network
    Input: network and value of lambda
    Output: efficiency of the network as variable E0
    '''

    # Calculate average shortest path in our network of interest
    avg_length = nx.average_shortest_path_length(grafo)

    # Get the number of edges in our network
    num_edges = grafo.number_of_edges()

    # Apply formula to calculate energy
    E0 = lambd * avg_length + ((1 - lambd) * num_edges)

    return(E0)

def modify_connection(grafo, a, b):
    '''
    Function to check if two nodes are connected, if they are, remove the connection
        if they are not, create a connection. Check efficiency of the graph agains the 
        initial one
    Input: initial graph, two nodes to evaluate (a, b)
    Output: modified graph
    '''

    # Check if the node exists
    connection = grafo.has_edge(a, b)

    # If both nodes are connected, remove the edge
    if connection is True:
        grafo.remove_edge(a, b)
        # Inform the user of the action taken
        print('Removed edge between {} and {}'.format(a, b))
    # If the nodes are not connected, add an edge
    else:
        grafo.add_edge(a, b)
        # Inform the user of the action taken
        print('Added edge between {} and {}'.format(a, b))

    return(connection, grafo)


def return_to_ori(connection, grafo, a, b):
    '''
    Function to remove the network modifications
    Input: network, boolean checking the original connection, nodes of interest
    Ouptut: network in the original state
    '''

    # Remember, connection refers to the original network
    if connection is True:
        grafo.add_edge(a, b)
    else:
        grafo.remove_edge(a, b)

    return(grafo)

def main():
    '''
    Function to implement the minimization procedure and create an optimized network
    Input:
    Ouptut: Optimized network with reduced efficiency
    '''

    # Create a random graph with 50 nodes
    nodes = 50
    grafo = random_network(nodes)

    # Calculate efficiency
    for lambd in [0.01, 0.08, 0.99]:
        E0 = efficiency(grafo, lambd)

        print('Original efficiency for lambda = {} : {}'.format(lambd, E0))
        
        # Draw original network
        plt.figure(figsize = (8, 6))
        plt.title('Original Random network with lambda = {} and efficiency = {}'.format(lambd, E0))
        nx.draw(grafo, with_labels = True, node_size = 20, node_color = 'blue')

        plt.savefig('Random_L{}_E{}.png'.format(lambd, E0))

        # Modify graph and check efficiency
        for i in range(100):
            a = random.randrange(nodes)
            b = random.randrange(nodes)

            connection, grafo = modify_connection(grafo, a, b)

            # Check energy and compare with initial energy
            E1 = efficiency(grafo, lambd)

            # If E1 > E0, undo modifications
            if E1 >= E0:

                grafo = return_to_ori(connection, grafo, a, b)
            else:
                # Keep modifications and the new network energy
                E0 = E1
        
        # Print final efficiency
        print('Final efficiency for lambda = {} : {}'.format(lambd, E0))

        # Draw final efficient network
        plt.figure(figsize = (8, 6))
        plt.title('Final optimized network with lambda = {} and efficiency = {}'.format(lambd, E0))
        nx.draw(grafo, with_labels = True, node_size = 20, node_color = 'green')

        plt.savefig('Optimized_L{}_E{}.png'.format(lambd, E0))



if __name__ == '__main__':
    main()
 
