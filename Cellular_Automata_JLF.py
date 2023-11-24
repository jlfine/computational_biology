# Cellular Automata for Conway's Game of Life
# Jacob L. Fine
# August 31st, 2023
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def game(board):  # applies the function to the current board at each iteration

    updated_board = np.zeros(board.shape)  # generates the board, empty at first
    for r, c in np.ndindex(board.shape): # goes through each cell in the board, based on iterating through the indices
        alive_neighbors_count = np.sum(board[r-1:r+2, c-1:c+2]) - board[r,c]  # for a given cell, get the count of its alive neighbors
        if board[r,c] == 1: # i.e., if the given cell is alive
            if (alive_neighbors_count < 2) or (alive_neighbors_count > 3): # if it has less than 2 or more than 3 neighbors in the current round, i.e., under/overpop
                updated_board[r,c] = 0
            elif (alive_neighbors_count == 2) or (alive_neighbors_count == 3):  # if it has 2 or 3 neighbors in the current round, then it stays alive
                updated_board[r,c] = 1
        elif board[r,c] == 0: # i.e., if the cell is dead in the current round
            if (alive_neighbors_count == 3): # if its currently dead but has 3 live neighbors, in the next round it will become alive
                updated_board[r,c] = 1
            else:
                updated_board[r,c] = 0  # otherwise, dead cells stay dead in the next round
    board = updated_board  # changes the whole board to the updated board
    return board
    
board = np.random.randint(2, size=(100, 100))  # initializes a board to have random integers between 0 and 1 everywhere
board = np.zeros((100,100))

# different initializations
def spaceship_1():
    for i in range(30,60,3):
        board[i,i] = 1  # updates one entry to 1
        board[i,i+1] = 1  # updates one entry to 1
        board[i,i+2] = 1  # updates one entry to 1
        board[i-1,i+2] = 1  # updates one entry to 1
        board[i-2,i+1] = 1  # updates one entry to 1


def spaceship_2():
    for i in range(30,60,2):
        board[i,i] = 1  # updates one entry to 1
        board[i,i+1] = 1  # updates one entry to 1
        board[i,i+2] = 1  # updates one entry to 1
        board[i-1,i+2] = 1  # updates one entry to 1
        board[i-2,i+1] = 1  # updates one entry to 1

spaceship_1()  # calls spaceship_1 initialization

steps = 100

fig, ax = plt.subplots()  # initializes the figure to display

im = ax.imshow(board, cmap='binary')  # makes an image out of the matrix

def update(i):
    global board  # references the global variable
    board = game(board)  # gets the board from the game function
    im.set_data(board)  # updates the image according to the current board

ani = animation.FuncAnimation(fig, update, frames=steps, interval=200)  # makes an animation out of the game

plt.show()  # shows the result of the animation


