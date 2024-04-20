# Jacob L. Fine
# April 20th, 2024

# Tic-tac-toe is an example of an M,n,k game (with M,n,k = 3) (see: https://en.wikipedia.org/wiki/M,n,k-game) with properties that may be explored computationally. 
# Connect 4 is another examlpe of an M,n,k game with (m=6,n=7,k=4)
# In M,n,k games, each player is a differnet symbol must obtain k symbols of their type in a row, column, or diagonal, on a matrix that is m x n. 
# This code will compute the board states of all possible tic-tac-toe games, and returns the number of unqiue games and the winners of each game (first mover, second mover). 
# It can be shown that there exist 255168 possible unique tic-tac-toe games, and that the first mover has the advantage. 
# Specifically, if the first mover is 'O', we will obtain the win_O : win_X : draw ratio to be 131184 : 77904 : 46080
# This code will therefore reproduce the results of Table 2 of Z.M.K. Zuhri (2022); ref: https://informatika.stei.itb.ac.id/~rinaldi.munir/Matdis/2021-2022/Makalah2021/Makalah-Matdis-2021%20(148).pdf

class TicTacToeGame:
    # initializes the board, sets the size and chooses the first player (set to O)
    def __init__(self, board_size=3, first_player='O'):
        self.BOARD_SIZE = board_size
        self.first_player = first_player
        self.game_states = []
        self.scores = {}

    # a function to check each possible place in the board where there could be a winner: each row, each col, and each diagonal
    def check_if_board_is_winning(self, board, player):
        # checks the ROWS for a winner, and sets check_winner to True
        for row in board:
            if all(cell == player for cell in row):
                return True

        # checks the COLS for a winner, and sets check_winner to True
        for col in range(self.BOARD_SIZE):
            if all(board[row][col] == player for row in range(self.BOARD_SIZE)):
                return True

        # checks the DIAGS for a winner, and sets check_winner to True
        if all(board[i][i] == player for i in range(self.BOARD_SIZE)) or \
           all(board[i][self.BOARD_SIZE - 1 - i] == player for i in range(self.BOARD_SIZE)):
            return True

        return False

    # checks if the current board is full, returns True if so
    def check_if_board_full(self, board):
        return all(cell != ' ' for row in board for cell in row)

    # generates all possible TTT game states, given the first initalized empty 3x3 board and the first player, using a helper function
    def make_all_game_states(self):
        # calls the helper function first at the current state being an empty 3x3 board
        self.generate_states_helper([[' ' for _ in range(self.BOARD_SIZE)] for _ in range(self.BOARD_SIZE)], self.first_player)
    
    # uses a helper function to check if the current state has a winner, after checking, it generates a copy of the previous state, and updates the current row, col value with the player letter (X/O)
    def generate_states_helper(self, current_state, player):

        # appends the final game state and the current winner to the list of final game states; converts the matrix representations of boards to strings for storing final board states
        if self.check_if_board_is_winning(current_state, 'O'):  
            self.game_states.append((''.join(cell for row in current_state for cell in row), 'O'))
            return
        elif self.check_if_board_is_winning(current_state, 'X'):
            self.game_states.append((''.join(cell for row in current_state for cell in row), 'X'))
            return
        
        # if the board is full but no winner has been attained, it will return a draw
        elif self.check_if_board_full(current_state):  
            self.game_states.append((''.join(cell for row in current_state for cell in row), 'Draw'))
            return
        # goes through all possible starting positions, and explores game states that stem from that starting position
        for row in range(self.BOARD_SIZE):
            for col in range(self.BOARD_SIZE):

                # checks if the current row,col value is blank in the TTT board
                if current_state[row][col] == ' ':
                    # generates a copy of the current state, and updates it with the player logo, if this board happens to be a winning board it will not overwrite the current_state
                    new_state = [list(row) for row in current_state]
                    new_state[row][col] = player
                    # will use recrusion to explore all new_state paths for a given starting position until the board is filled for that starting position
                    self.generate_states_helper(new_state, 'X' if player == 'O' else 'O')

    # computes the scores
    def compute_winners(self):
        self.scores = {player: sum(1 for _, winner in self.game_states if winner == player) for player in ['O', 'X', 'Draw']}
        print('Outcomes:')
        print(self.scores)

    # computes the number of unique games
    def print_winners(self):
        total_unique_games = sum(self.scores.values())
        print('Total unique games:')
        print(total_unique_games)

# create an instance of the tic tac toe game
tic_tac_toe = TicTacToeGame()

# generate unique game states
tic_tac_toe.make_all_game_states()

# calculate scores
tic_tac_toe.compute_winners()
tic_tac_toe.print_winners()
