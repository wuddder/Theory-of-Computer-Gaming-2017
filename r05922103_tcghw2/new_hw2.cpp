/*
 * This code is provied as a sample code of Hw 2 of "Theory of Computer Game".
 * The "genmove" function will randomly output one of the legal move.
 * This code can only be used within the class.
 *
 * 2015 Nov. Hung-Jui Chang
 * */
#include <iostream>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <cmath>
#include <vector>

#include <stdio.h>
#include <fstream>

#define BOARDSIZE        9
#define BOUNDARYSIZE    11
#define COMMANDLENGTH 1000
#define DEFAULTTIME     10
#define SIMULATETIME     1
#define DEFAULTKOMI      7

#define SAFEMOVE        10
#define UCBCONSTANT    1.4
#define TRIALS          10 

#define MAXGAMELENGTH 1000
#define MAXSTRING       50
#define MAXDIRECTION     4

#define NUMINTERSECTION 81
#define HISTORYLENGTH  200

#define EMPTY            0
#define BLACK            1
#define WHITE            2
#define BOUNDARY         3

#define SELF             1
#define OPPONENT         2

#define NUMGTPCOMMANDS      15

#define LOCALVERSION      1
#define GTPVERSION        2
 
using namespace std;
using namespace std::chrono;

int _board_size = BOARDSIZE;
int _board_boundary = BOUNDARYSIZE;
double _komi =  DEFAULTKOMI;
bool _newgame = true;
const int DirectionX[MAXDIRECTION] = {-1, 0, 1, 0};
const int DirectionY[MAXDIRECTION] = { 0, 1, 0,-1};
const char LabelX[]="0ABCDEFGHJ";
int move_counter = 0;

typedef struct mcsNode{
    struct mcsNode* parent;
    vector <struct mcsNode> child;
    double score;
    double ucb_value;
    double square_score;
    int visit; 
    int move;
    int turn;
    int Board[BOUNDARYSIZE][BOUNDARYSIZE];
} MCSNODE;

void reset(int Board[BOUNDARYSIZE][BOUNDARYSIZE]);
int find_liberty(int X, int Y, int label, int Board[BOUNDARYSIZE][BOUNDARYSIZE], int ConnectBoard[BOUNDARYSIZE][BOUNDARYSIZE]);
void count_liberty(int X, int Y, int Board[BOUNDARYSIZE][BOUNDARYSIZE], int Liberties[MAXDIRECTION]);
void count_neighboorhood_state(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int X, int Y, int turn, int* empt, int* self, int* oppo ,int* boun, int NeighboorhoodState[MAXDIRECTION]);
int remove_piece(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int X, int Y, int turn);
void update_board(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int X, int Y, int turn);
int update_board_check(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int X, int Y, int turn);
int gen_legal_move(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int turn, int game_length, int GameRecord[MAXGAMELENGTH][BOUNDARYSIZE][BOUNDARYSIZE], int MoveList[HISTORYLENGTH]);
int rand_pick_move(int num_legal_moves, int MoveList[HISTORYLENGTH]);
int rand_simulate(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int turn, int game_length, int GameRecord[MAXGAMELENGTH][BOUNDARYSIZE][BOUNDARYSIZE]);
int rand_simulate_score(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int turn, int game_length, int GameRecord[MAXGAMELENGTH][BOUNDARYSIZE][BOUNDARYSIZE]);
int MCS_pure_move(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int num_legal_moves, int MoveList[HISTORYLENGTH], clock_t end_t, int turn, int game_length, int GameRecord[MAXGAMELENGTH][BOUNDARYSIZE][BOUNDARYSIZE]);
float UCBvalue(mcsNode* node);
void do_move(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int turn, int move);
void record(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int GameRecord[MAXGAMELENGTH][BOUNDARYSIZE][BOUNDARYSIZE], int game_length);
double final_score(int Board[BOUNDARYSIZE][BOUNDARYSIZE]);
void gtp_showboard(int Board[BOUNDARYSIZE][BOUNDARYSIZE]);
/*
 * This function reset the board, the board intersections are labeled with 0,
 * the boundary intersections are labeled with 3.
 * */
void reset(int Board[BOUNDARYSIZE][BOUNDARYSIZE]) {
    for (int i = 1 ; i <= BOARDSIZE; ++i) {
    for (int j = 1 ; j <= BOARDSIZE; ++j) {
        Board[i][j] = EMPTY;
    }
    }
    for (int i = 0 ; i < BOUNDARYSIZE; ++i) {
    Board[0][i] = Board[BOUNDARYSIZE-1][i] = Board[i][0] = Board[i][BOUNDARYSIZE-1] = BOUNDARY;
    }
}

/*
 * This function return the total number of liberity of the string of (X, Y) and
 * the string will be label with 'label'.
 * */
int find_liberty(int X, int Y, int label, int Board[BOUNDARYSIZE][BOUNDARYSIZE], int ConnectBoard[BOUNDARYSIZE][BOUNDARYSIZE]) {
    // Label the current intersection
    ConnectBoard[X][Y] |= label;
    int total_liberty = 0;
    for (int d = 0 ; d < MAXDIRECTION; ++d) {
        // Check this intersection has been visited or not
        if ((ConnectBoard[X+DirectionX[d]][Y+DirectionY[d]] & (1<<label) )!= 0) continue;

        // Check this intersection is not visited yet
        ConnectBoard[X+DirectionX[d]][Y+DirectionY[d]] |=(1<<label);
        // This neighboorhood is empty
        if (Board[X+DirectionX[d]][Y+DirectionY[d]] == EMPTY)
            total_liberty++;
        // This neighboorhood is self stone
        else if (Board[X+DirectionX[d]][Y+DirectionY[d]] == Board[X][Y])
            total_liberty += find_liberty(X+DirectionX[d], Y+DirectionY[d], label, Board, ConnectBoard);
    }
    return total_liberty;
}

/*
 * This function count the liberties of the given intersection's neighboorhod
 * */
void count_liberty(int X, int Y, int Board[BOUNDARYSIZE][BOUNDARYSIZE], int Liberties[MAXDIRECTION]) {
    int ConnectBoard[BOUNDARYSIZE][BOUNDARYSIZE];
    // Initial the ConnectBoard
    for (int i = 0 ; i < BOUNDARYSIZE; ++i) {for (int j = 0 ; j < BOUNDARYSIZE; ++j) {ConnectBoard[i][j] = 0;}}
    // Find the same connect component and its liberity
    for (int d = 0 ; d < MAXDIRECTION; ++d) {
        Liberties[d] = 0;
        if (Board[X+DirectionX[d]][Y+DirectionY[d]] == BLACK || Board[X+DirectionX[d]][Y+DirectionY[d]] == WHITE ) {
            Liberties[d] = find_liberty(X+DirectionX[d], Y+DirectionY[d], d, Board, ConnectBoard);
        }
    }
}

/*
 * This function count the number of empty, self, opponent, and boundary intersections of the neighboorhod
 * and saves the type in NeighboorhoodState.
 * */
void count_neighboorhood_state(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int X, int Y, int turn, int* empt, int* self, int* oppo ,int* boun, int NeighboorhoodState[MAXDIRECTION]) {
    for (int d = 0 ; d < MAXDIRECTION; ++d) {
    // check the number of nonempty neighbor
    switch(Board[X+DirectionX[d]][Y+DirectionY[d]]) {
        case EMPTY:    (*empt)++; 
               NeighboorhoodState[d] = EMPTY;
               break;
        case BLACK:    if (turn == BLACK) {
                   (*self)++;
                   NeighboorhoodState[d] = SELF;
               }
               else {
                   (*oppo)++;
                   NeighboorhoodState[d] = OPPONENT;
               }
               break;
        case WHITE:    if (turn == WHITE) {
                   (*self)++;
                   NeighboorhoodState[d] = SELF;
               }
               else {
                   (*oppo)++;
                   NeighboorhoodState[d] = OPPONENT;
               }
               break;
        case BOUNDARY: (*boun)++;
               NeighboorhoodState[d] = BOUNDARY;
               break;
    }
    }
}

/*
 * This function remove the connect component contains (X, Y) with color "turn" 
 * And return the number of remove stones.
 * */
int remove_piece(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int X, int Y, int turn) {
    int remove_stones = (Board[X][Y]==EMPTY)?0:1;
    Board[X][Y] = EMPTY;
    for (int d = 0; d < MAXDIRECTION; ++d) {
    if (Board[X+DirectionX[d]][Y+DirectionY[d]] == turn) {
        remove_stones += remove_piece(Board, X+DirectionX[d], Y+DirectionY[d], turn);
    }
    }
    return remove_stones;
}
/*
 * This function update Board with place turn's piece at (X,Y).
 * Note that this function will not check if (X, Y) is a legal move or not.
 * */
void update_board(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int X, int Y, int turn) {
    int num_neighborhood_self = 0;
    int num_neighborhood_oppo = 0;
    int num_neighborhood_empt = 0;
    int num_neighborhood_boun = 0;
    int Liberties[4];
    int NeighboorhoodState[4];
    count_neighboorhood_state(Board, X, Y, turn,
        &num_neighborhood_empt,
        &num_neighborhood_self,
        &num_neighborhood_oppo,
        &num_neighborhood_boun, NeighboorhoodState);
    // check if there is opponent piece in the neighboorhood
    if (num_neighborhood_oppo != 0) {
    count_liberty(X, Y, Board, Liberties);
    for (int d = 0 ; d < MAXDIRECTION; ++d) {
        // check if there is opponent component only one liberty
        if (NeighboorhoodState[d] == OPPONENT && Liberties[d] == 1 && Board[X+DirectionX[d]][Y+DirectionY[d]]!=EMPTY) {
        remove_piece(Board, X+DirectionX[d], Y+DirectionY[d], Board[X+DirectionX[d]][Y+DirectionY[d]]);
        }
    }
    }
    Board[X][Y] = turn;
}
/*
 * This function update Board with place turn's piece at (X,Y).
 * Note that this function will check if (X, Y) is a legal move or not.
 * */
int update_board_check(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int X, int Y, int turn) {
    // Check the given coordination is legal or not
    if ( X < 1 || X > BOARDSIZE || Y < 1 || Y > BOARDSIZE || Board[X][Y]!=EMPTY)
    return 0;
    int num_neighborhood_self = 0;
    int num_neighborhood_oppo = 0;
    int num_neighborhood_empt = 0;
    int num_neighborhood_boun = 0;
    int Liberties[4];
    int NeighboorhoodState[4];
    count_neighboorhood_state(Board, X, Y, turn,
        &num_neighborhood_empt,
        &num_neighborhood_self,
        &num_neighborhood_oppo,
        &num_neighborhood_boun, NeighboorhoodState);
    // Check if the move is a legal move
    // Condition 1: there is a empty intersection in the neighboorhood
    int legal_flag = 0;
    count_liberty(X, Y, Board, Liberties);
    if (num_neighborhood_empt != 0) {
    legal_flag = 1;
    }
    else {
        // Condition 2: there is a self string has more than one liberty
        for (int d = 0; d < MAXDIRECTION; ++d) {
            if (NeighboorhoodState[d] == SELF && Liberties[d] > 1) {
            legal_flag = 1;
            }
        }
        if (legal_flag == 0) {
        // Condition 3: there is a opponent string has exactly one liberty
            for (int d = 0; d < MAXDIRECTION; ++d) {
            if (NeighboorhoodState[d] == OPPONENT && Liberties[d] == 1) {
                legal_flag = 1;
            }
            }
        }
    }

    if (legal_flag == 1) {
        // check if there is opponent piece in the neighboorhood
        if (num_neighborhood_oppo != 0) {
            for (int d = 0 ; d < MAXDIRECTION; ++d) {
                // check if there is opponent component only one liberty
                if (NeighboorhoodState[d] == OPPONENT && Liberties[d] == 1 && Board[X+DirectionX[d]][Y+DirectionY[d]]!=EMPTY){
                    remove_piece(Board, X+DirectionX[d], Y+DirectionY[d], Board[X+DirectionX[d]][Y+DirectionY[d]]);
                }
            }
        }
        Board[X][Y] = turn;
    }

    return (legal_flag==1)?1:0;
}

/*
 * This function return the number of legal moves with clor "turn" and
 * saves all legal moves in MoveList
 * */
int gen_legal_move(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int turn, int game_length, int GameRecord[MAXGAMELENGTH][BOUNDARYSIZE][BOUNDARYSIZE], int MoveList[HISTORYLENGTH]) {
    int NextBoard[BOUNDARYSIZE][BOUNDARYSIZE];
    int num_neighborhood_self = 0;
    int num_neighborhood_oppo = 0;
    int num_neighborhood_empt = 0;
    int num_neighborhood_boun = 0;
    int legal_moves = 0;
    int next_x, next_y;
    int Liberties[4];
    int NeighboorhoodState[4];
    bool eat_move = 0;
    int num_eat = 0;;
    for (int x = 1 ; x <= BOARDSIZE; ++x) {
    for (int y = 1 ; y <= BOARDSIZE; ++y) {
        if(game_length < SAFEMOVE &&( x == 1 || y == 1 || x == BOARDSIZE || y == BOARDSIZE))continue;
        // check if current 
        if (Board[x][y] == 0) {
            // check the liberty of the neighborhood intersections
            num_neighborhood_self = 0; num_neighborhood_oppo = 0;
            num_neighborhood_empt = 0; num_neighborhood_boun = 0;
            // count the number of empy, self, opponent, and boundary neighboorhood
            count_neighboorhood_state(Board, x, y, turn, &num_neighborhood_empt, &num_neighborhood_self, &num_neighborhood_oppo, &num_neighborhood_boun, NeighboorhoodState);
            // check if the emtpy intersection is a legal move
            next_x = next_y = 0;
            eat_move = 0;
            count_liberty(x, y, Board, Liberties);
            // Case 1: exist empty intersection in the neighborhood
            if (num_neighborhood_empt > 0) {
             next_x = x;
             next_y = y;
             // check if it is a capture move
            for (int d = 0 ; d < MAXDIRECTION; ++d) {
                if (NeighboorhoodState[d] == OPPONENT && Liberties[d] == 1) {
                    eat_move = 1;
                }   
            }

            }
            // Case 2: no empty intersection in the neighborhood
            else {
                // Case 2.1: Surround by the self piece
                if (num_neighborhood_self + num_neighborhood_boun == MAXDIRECTION) {
                    int check_flag = 0, check_eye_flag = num_neighborhood_boun;
                    for (int d = 0 ; d < MAXDIRECTION; ++d) {
                        // Avoid fill self eye
                        if (NeighboorhoodState[d]==SELF && Liberties[d] > 1)
                            check_eye_flag++;
                        // Check if there is one self component which has more than one liberty
                        if (NeighboorhoodState[d]==SELF && Liberties[d] > 1)
                            check_flag = 1;
                    }
                    if (check_flag == 1 && check_eye_flag!=4) {
                        next_x = x;
                        next_y = y;
                    }
                }   
                // Case 2.2: Surround by opponent or both side's pieces.
                else if (num_neighborhood_oppo > 0) {
                    int check_flag = 0;
                    int eat_flag = 0;
                    for (int d = 0 ; d < MAXDIRECTION; ++d) {
                        // Check if there is one self component which has more than one liberty
                        if (NeighboorhoodState[d]==SELF && Liberties[d] > 1)
                            check_flag = 1;
                        // Check if there is one opponent's component which has exact one liberty
                        if (NeighboorhoodState[d]==OPPONENT && Liberties[d] == 1) 
                            eat_flag = 1;
                    }
                    if (check_flag == 1) {
                        next_x = x;
                        next_y = y;
                        if (eat_flag == 1)
                            eat_move = 1;
                    }
                    else { // check_flag == 0
                        if (eat_flag == 1) {
                            next_x = x;
                            next_y = y;
                            eat_move = 1;
                        }
                    }
                }   
            }
            if (next_x !=0 && next_y !=0) {
                // copy the current board to next board
                for (int i = 0 ; i < BOUNDARYSIZE; ++i)for (int j = 0 ; j < BOUNDARYSIZE; ++j)NextBoard[i][j] = Board[i][j];
                // do the move
                // The move is a capture move and the board needs to be updated.
                if (eat_move == 1) 
                    update_board(NextBoard, next_x, next_y, turn);
                else 
                    NextBoard[x][y] = turn;
                // Check the history to avoid the repeat board
                bool repeat_move = 0;
                for (int t = 0 ; t < game_length; ++t) {
                    bool repeat_flag = 1;
                    for (int i = 1; i <=BOARDSIZE; ++i) {
                    for (int j = 1; j <=BOARDSIZE; ++j) {
                        if (NextBoard[i][j] != GameRecord[t][i][j])
                            repeat_flag = 0;
                    }
                    }
                    if (repeat_flag == 1) {
                        repeat_move = 1;
                        break;
                    }
                }
                if (repeat_move == 0) {
                    // 3 digit zxy, z means eat or not, and put at (x, y)
                    MoveList[legal_moves] = eat_move * 100 + next_x * 10 + y ;
                    legal_moves++;
                    if(eat_move == 1) num_eat++;
                }
            }
        }
    }
    }
   return legal_moves;
}
/*
 * This function randomly selects one move from the MoveList.
 * */
int rand_pick_move(int move_count, int MoveList[HISTORYLENGTH]) {
    if (move_count == 0)
    return 0;
    else {
    int move_id = rand()%move_count;
    return MoveList[move_id];
    }
}
/*
 * This function simlates the rand_move a single game.
 * */
int rand_simulate(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int turn, int game_length, int GameRecord[MAXGAMELENGTH][BOUNDARYSIZE][BOUNDARYSIZE]){

    int original_turn = turn;
    int MoveList[HISTORYLENGTH];
    int move_count;
    int current_move;
    bool hasPassed[2] = {false, false};
    
    while (!(hasPassed[0] && hasPassed[1])) {
        turn = turn % 2 + 1;
        // generate all legal move and get legal move counts
        move_count = gen_legal_move(Board, turn, game_length, GameRecord, MoveList);

        // if there is no move, pass is the only choice
        if (move_count == 0){
            hasPassed[turn-1] = true;
        } else {
            hasPassed[turn-1] = false;
        }

        // when both player have passed, game end
        if (hasPassed[0] && hasPassed[1]) break;

        // before game end, pick a random move and go on
        current_move = rand_pick_move(move_count, MoveList);
        do_move(Board, turn, current_move);
        record(Board, GameRecord, game_length+1);
        game_length += 1;
    }
    // game end, get the final result of the game and reuturn
    int result = final_score(Board) - _komi;
    if(result == 0) 
        return 0;
    else if(result > 0 ^ original_turn == BLACK) 
        return -1;
    else 
        return 1;

}
double rand_simulate_score(int current_board[BOUNDARYSIZE][BOUNDARYSIZE], int current_turn, int current_depth){

    vector<int> rand_array;
    for(int i=0; i < BOARDSIZE*BOARDSIZE; i++) 
        rand_array.push_back(10 * (i / 9) + i % 9 + 11);
    double result;
    int move[2] = {-1, -1};
    bool move_check;

    while(!(move[0] == 0 && move[1] == 0) && current_depth < HISTORYLENGTH){
        random_shuffle(rand_array.begin(), rand_array.end());
        move_check = false;
        for (int k = 0; k < rand_array.size(); k++){
            if(move[current_turn] == rand_array[k]) continue;
            if ((update_board_check(current_board, rand_array[k]/10, rand_array[k]%10, current_turn)) == 1){
                move_check = true;
                move[current_turn - 1] = rand_array[k];
                break;
            }
        }
        if (!move_check)
            move[current_turn-1] = 0;
        current_turn = (current_turn%2)+1;
        current_depth ++;
    }
    result = final_score(current_board);
    return result;


}
bool CheckIfEat(int previousBoard[BOUNDARYSIZE][BOUNDARYSIZE], int currentBoard[BOUNDARYSIZE][BOUNDARYSIZE], int turn) {
    int previousCount = 0;
    int currentCount = 0;
    int color = (turn == BLACK) ? WHITE : BLACK;
    for (int i = 0; i < BOUNDARYSIZE; i++) {
        for (int j = 0; j < BOUNDARYSIZE; j++) {
            if (previousBoard[i][j] == color) {
                previousCount++;
            }
        }
    }
    for (int i = 0; i < BOUNDARYSIZE; i++) {
        for (int j = 0; j < BOUNDARYSIZE; j++) {
            if (currentBoard[i][j] == color) {
                currentCount++;
            }
        }
    }
    return (previousCount > currentCount) ? false : true;
}

/*
 * This function selects one move with Monte-Carlo pure search.
 * */
int MCS_pure_move(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int move_count, int MoveList[HISTORYLENGTH], clock_t end_t, int turn, int game_length, int GameRecord[MAXGAMELENGTH][BOUNDARYSIZE][BOUNDARYSIZE]){
    // no move, return 
    if (move_count == 0) return 0;
    int random_points[move_count];

    // initialize
    for (int i = 0; i < move_count; i++) {
        random_points[i] = 0;
    }
    // during thinking timeslot, go as much random calculate trial as possible
    while (clock() < end_t) {
        // every move
        for (int i = 0; i < move_count; i++) {
            // set the hole board
            for (int m = 0; m < BOUNDARYSIZE; m++) {
                for (int n = 0; n < BOUNDARYSIZE; n++) {
                    GameRecord[game_length][m][n] = Board[m][n];
                }
            }
            do_move(GameRecord[game_length+1], turn, MoveList[i]);
            random_points[i] += rand_simulate(GameRecord[game_length+1], turn, game_length+1, GameRecord);
        }
    } 

    // after counting points, choose one have best performance
    const int N = sizeof(random_points) / sizeof(int);
    int max_position = distance(random_points, max_element(random_points, random_points + N));

    return MoveList[max_position];
}

double getUCB(MCSNODE* node, int position) {
    return (node -> child[position].score)/(node -> child[position].visit) + UCBCONSTANT * sqrt(log(node -> visit)/node -> child[position].visit);
}

int MCTS_move(int Board[BOUNDARYSIZE][BOUNDARYSIZE], steady_clock::time_point end_t, int turn, int game_length,int GameRecord[MAXGAMELENGTH][BOUNDARYSIZE][BOUNDARYSIZE]){
    /*
    In every MCTS_move, we have to expand a MCTS tree from the root node,
    and find out the optimal next step we wanted
    */
    
    /* For principal variation */

    ofstream myfile;
    if (_newgame) {
         myfile.open("Principal_variation.txt");
         _newgame = false;
    } else {
        myfile.open( "Principal_variation.txt", ios::out | ios::app );
    }
    string principal_variation;
    move_counter += 1;
    principal_variation += "========================" + to_string(move_counter) + "=======================\n";

    /* initialize root node */
    MCSNODE root;

    root.parent = NULL;
    root.score = 0.0;
    root.square_score = 0.0;
    root.visit = 0;
    root.move = 0;
    root.turn = turn;
    for(int x = 0; x < BOUNDARYSIZE; x++)
        for(int y = 0; y < BOUNDARYSIZE; y++)
            root.Board[x][y] = Board[x][y];
    /* current pointer to root node */   
    MCSNODE* current;
    current = (&root);
    
    /* for choosing the max UCB node */ 
    int max_position;
    double max_UCB,tmp_UCB;
    /* parameters for progressive pruning */ 
    double mu[BOARDSIZE*BOARDSIZE],sigma[BOARDSIZE*BOARDSIZE],Mmr[BOARDSIZE*BOARDSIZE];
    double gamar_d = 2.0;
    double sigma_e = 0.2;
    double max_Mml,tmp_Mml;

    /* parameters for random simulation */ 
    int depth, current_depth, current_turn;
    bool move_check;
    int move[2];
    int current_board[BOUNDARYSIZE][BOUNDARYSIZE];
    vector<int> rand_array;
    for(int i=0; i < BOARDSIZE*BOARDSIZE; i++) 
        rand_array.push_back(10 * (i / 9) + i % 9 + 11);
    
    /* parameters for expansion */
    int move_count;
    int MoveList[HISTORYLENGTH];
    bool first = true;
    
    /* parameters for back probagation */
    double visit, score, square_score, total_visit, total_score, total_square_score, result;

    
    while(steady_clock::now() < end_t){
        depth = game_length;
        current = &root;

        /* selection: find the node with maximum UCB score*/
        while(current -> child.size() != 0){
            first = false;
            max_position = rand() % current -> child.size();
            max_UCB = getUCB(current, max_position);
            for(int i = 0; i < current -> child.size(); i++){
                tmp_UCB = getUCB(current, i);
                if(tmp_UCB > max_UCB){
                    max_UCB = tmp_UCB;
                    max_position = i;
                }
                /* calculate parameters for progressive pruning */
                mu[i] = (current -> child[i].score) / (current -> child[i].visit);
                sigma[i] = (current -> child[i].square_score - (2 * mu[i] * current -> child[i].score)) / (current -> child[i].visit) + mu[i]*mu[i];
                Mmr[i] = mu[i] * (mu[i] + gamar_d * sigma[i]);
                tmp_Mml = mu[i] * (mu[i] - gamar_d * sigma[i]);
                if (i == 0) max_Mml = tmp_Mml;
                else if (tmp_Mml > max_Mml)
                    max_Mml = tmp_Mml;
            }
            /* pruning */
            int index = 0;
            for(int i = 0; i < current -> child.size(); i++, index++){
                if(mu[i] * sigma[i] < sigma_e) continue;
                if(Mmr[index] < max_Mml){
                    current -> child.erase(current -> child.begin(), current -> child.begin()+i);
                    i--;
                }
            }
            current = &(current -> child[max_position]);
            record(current -> Board, GameRecord, depth);
            depth ++;
            principal_variation += "move" + to_string(depth-game_length) + ": " + to_string(current -> move) + " score: " + to_string(max_UCB) + " simulation time: " + to_string(current -> parent ->  visit) + "\n";
        }
        
        /* expansion: expand for the new node */
        move_count = gen_legal_move(current -> Board, current -> turn, depth, GameRecord, MoveList);
        if (move_count == 0) break;
        for (int i = 0; i < move_count; i++){
            if (MoveList[i] >= 100 && first) {
                return MoveList[i];
            }
            MCSNODE child_node;
            child_node.parent = current;
            child_node.score = 0.0;
            child_node.square_score = 0.0;
            child_node.visit = 0;
            child_node.move = MoveList[i];
            for(int x = 0; x < BOUNDARYSIZE; x++)
                for(int y = 0; y < BOUNDARYSIZE; y++)
                    child_node.Board[x][y] = current -> Board[x][y];
            do_move(child_node.Board, current -> turn, child_node.move);
            child_node.turn = (current -> turn % 2) + 1;
            current -> child.push_back(child_node);
        }
    
        /* simulation: rollout and get the score */
        total_visit = 0.0; total_score = 0.0; total_square_score = 0.0;
        for(int i = 0; i < current -> child.size(); i++){
                visit = 0.0; score = 0.0; square_score = 0.0;
                int count = 0;
                while (count < TRIALS) {
                    for(int x = 0; x < BOUNDARYSIZE; x++)
                        for(int y = 0; y < BOUNDARYSIZE; y++)
                            current_board[x][y] = current -> child[i].Board[x][y];
                    current_turn = current -> child[i].turn;
                    current_depth = depth;
                    count ++;                   
                    result = rand_simulate_score(current_board, current_turn, current_depth);
                    
                    if (result > 17) {
                        score = score + 1.0;
                        square_score = square_score + 1.0;
                    } else if (result > 10) {
                        score = score + 0.8;
                        square_score = square_score + 0.64;
                    } else if (result > 4) {
                        score = score + 0.6;
                        square_score = square_score + 0.36;
                    } else if (result > -3) {
                        score = score + 0.4;
                        square_score = square_score + 0.16;
                    } else {
                        score = score + 0.2;
                        square_score = square_score + 0.04;
                    }
                    visit = visit + 1;
                }
                if (current -> child[i].turn == BLACK) {
                    current -> child[i].score = score;
                    current -> child[i].square_score = square_score;
                } else {
                    current -> child[i].score = visit - score;
                    current -> child[i].square_score = visit - 2*score + square_score;
                }
                current -> child[i].visit = visit;
                total_score = total_score + score;
                total_square_score = total_square_score + square_score;
                total_visit = total_visit + visit;
        }

        /* back propagation: update the visit and score for the parent node */
        while(current != NULL) {
            if(current -> turn == BLACK){
                current -> score = current -> score + total_score;
                // current -> square_score = current -> square_score + total_square_score;
            } else {
                current -> score = current -> score + total_visit - total_score;
                // current -> square_score = current -> square_score + total_visit - 2*total_score + total_square_score;
            }
            current -> visit = current -> visit + total_visit;
            if(current -> parent == NULL){
                break;
            }
            current = current -> parent;
        }
    }

    if (root.child.size() == 0) return 0;
    max_UCB = root.child[0].score / root.child[0].visit;
    int return_move = root.child[0].move;
    for(int i = 0; i < root.child.size(); i++){
        tmp_UCB = root.child[i].score / root.child[i].visit;
        if(tmp_UCB > max_UCB){
            return_move = root.child[i].move;
            max_UCB = tmp_UCB;
        }
    }
    principal_variation += "================================================\n";
    principal_variation += "Move choosed: " + to_string(return_move) + "\n";
    principal_variation += "================================================\n";
    myfile << principal_variation;
    myfile.close();


    return return_move;
}

/*
 * This function update the Board with put 'turn' at (x,y)
 * where x = (move % 100) / 10 and y = move % 10.
 * Note this function will not check 'move' is legal or not.
 * */
void do_move(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int turn, int move) {
    int move_x = (move % 100) / 10;
    int move_y = move % 10;
    if (move<100)
        Board[move_x][move_y] = turn;
    else 
        update_board(Board, move_x, move_y, turn);
}
/* 
 * This function records the current game baord with current
 * game length "game_length"
 * */
void record(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int GameRecord[MAXGAMELENGTH][BOUNDARYSIZE][BOUNDARYSIZE], int game_length) {
        for (int i = 0 ; i < BOUNDARYSIZE; ++i) {
            for (int j = 0 ; j < BOUNDARYSIZE; ++j) {
            GameRecord[game_length][i][j] = Board[i][j];
            }
        }
}
/* 
 * This function randomly generate one legal move (x, y) with return value x*10+y,
 * if there is no legal move the function will return 0.
 * */
int genmove(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int turn, int time_limit, int game_length, int GameRecord[MAXGAMELENGTH][BOUNDARYSIZE][BOUNDARYSIZE]) {
    // clock_t start_t, end_t;
    // record start time
    // start_t = clock();
    // // calculate the time bound
    // end_t = start_t + CLOCKS_PER_SEC * time_limit;
    // int MoveList[HISTORYLENGTH];



    steady_clock::time_point start_t = steady_clock::now();
    steady_clock::time_point end_t = start_t + seconds(time_limit - 2);

    int return_move = 0;
    // int move_count = gen_legal_move(Board, turn, game_length, GameRecord, MoveList);

    
    for(int i=0;i<BOUNDARYSIZE;i++){
        Board[0][i] = Board[BOUNDARYSIZE-1][i] = Board[i][0] = Board[i][BOUNDARYSIZE-1] = BOUNDARY;
    }
    return_move = MCTS_move(Board, end_t, turn, game_length, GameRecord);
    // cout << "move chosen:" << return_move << endl;
    // return_move = rand_pick_move(num_legal_moves, MoveList);
    // return_move = MCS_pure_move(Board, move_count, MoveList, end_t, turn, game_length, GameRecord);

    do_move(Board, turn, return_move);
    //cerr << (clock() - start_t)/CLOCKS_PER_SEC << endl;
    return return_move % 100;
}
/*
 * This function counts the number of points remains in the board by Black's view
 * */
double final_score(int Board[BOUNDARYSIZE][BOUNDARYSIZE]) {
    int black, white;
    black = white = 0;
    int is_black, is_white;
    for (int i = 1 ; i <= BOARDSIZE; ++i) {
    for (int j = 1; j <= BOARDSIZE; ++j) {
        switch(Board[i][j]) {
        case EMPTY:
            is_black = is_white = 0;
            for(int d = 0 ; d < MAXDIRECTION; ++d) {
            if (Board[i+DirectionX[d]][j+DirectionY[d]] == BLACK) is_black = 1;
            if (Board[i+DirectionX[d]][j+DirectionY[d]] == WHITE) is_white = 1;
            }
            if (is_black + is_white == 1) {
            black += is_black;
            white += is_white;
            }
            break;
        case WHITE:
            white++;
            break;
        case BLACK:
            black++;
            break;
        }
    }
    }
    return black - white;
}
/* 
 * Folloscoreg are commands for Go Text Protocol (GTP)
 *
 * */
const char *KnownCommands[]={
    "protocol_version",
    "name",
    "version",
    "known_command",
    "list_commands",
    "quit",
    "boardsize",
    "clear_board",
    "komi",
    "play",
    "genmove",
    "undo",
    "quit",
    "showboard",
    "final_score"
};

void gtp_final_score(int Board[BOUNDARYSIZE][BOUNDARYSIZE]) {
    double result;
    result = final_score(Board);
    result -= _komi;
    cout << "= ";
    if (result > 0.0) { // Black score
    cout << "B+" << result << endl << endl<< endl;;
    }
    if (result < 0.0) { // White score
    cout << "W+" << -result << endl << endl<< endl;;
    }
    else { // draw
    cout << "0" << endl << endl<< endl;;
    }
}
void gtp_undo(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int game_length, int GameRecord[MAXGAMELENGTH][BOUNDARYSIZE][BOUNDARYSIZE]) {
    if (game_length!=0) {
    for (int i = 1; i <= BOARDSIZE; ++i) {
        for (int j = 1; j <= BOARDSIZE; ++j) {
        Board[i][j] = GameRecord[game_length][i][j];
        }
    }
    }
    cout << "= " << endl << endl;
}
void gtp_showboard(int Board[BOUNDARYSIZE][BOUNDARYSIZE]) {
    for (int i = 1; i <=BOARDSIZE; ++i) {
    cout << "#";
    cout <<10-i;
    for (int j = 1; j <=BOARDSIZE; ++j) {
        switch(Board[i][j]) {
        case EMPTY: cout << " .";break;
        case BLACK: cout << " X";break;
        case WHITE: cout << " O";break;
        }
    }
    cout << endl;
    }
    cout << "#  ";
    for (int i = 1; i <=BOARDSIZE; ++i) 
    cout << LabelX[i] <<" ";
    cout << endl;
    cout << endl;

}
void gtp_protocol_version() {
    cout <<"= 2"<<endl<< endl;
}
void gtp_name() {
    cout <<"= TCG-randomGo99" << endl<< endl;
}
void gtp_version() {
    cout << "= 1.02" << endl << endl;
}
void gtp_list_commands(){
    cout <<"= ";
    for (int i = 0 ; i < NUMGTPCOMMANDS; ++i) {
    cout <<KnownCommands[i] << endl;
    }
    cout << endl;
}
void gtp_known_command(const char Input[]) {
    for (int i = 0 ; i < NUMGTPCOMMANDS; ++i) {
    if (strcmp(Input, KnownCommands[i])==0) {
        cout << "= true" << endl<< endl;
        return;
    }
    }
    cout << "= false" << endl<< endl;
}
void gtp_boardsize(int size) {
    if (size!=9) {
    cout << "? unacceptable size" << endl<< endl;
    }
    else {
    _board_size = size;
   cout << "= "<<endl<<endl;
    }
}
void gtp_clear_board(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int NumCapture[]) {
    reset(Board);
    NumCapture[BLACK] = NumCapture[WHITE] = 0;
    cout << "= "<<endl<<endl;
}
void gtp_komi(double komi) {
    _komi = komi;
    cout << "= "<<endl<<endl;
}
void gtp_play(char Color[], char Move[], int Board[BOUNDARYSIZE][BOUNDARYSIZE], int game_length, int GameRecord[MAXGAMELENGTH][BOUNDARYSIZE][BOUNDARYSIZE]) {
    int turn, move_i, move_j;
    if (Color[0] =='b' || Color[0] == 'B')
    turn = BLACK;
    else
    turn = WHITE;
    if (strcmp(Move, "PASS") == 0 || strcmp(Move, "pass")==0) {
    record(Board, GameRecord, game_length+1);
    }
    else {
    // [ABCDEFGHJ][1-9], there is no I in the index.
    Move[0] = toupper(Move[0]);
    move_j = Move[0]-'A'+1;
    if (move_j == 10) move_j = 9;
    move_i = 10-(Move[1]-'0');
    update_board(Board, move_i, move_j, turn);
    record(Board, GameRecord, game_length+1);
    }
    cout << "= "<<endl<<endl;
}
void gtp_genmove(int Board[BOUNDARYSIZE][BOUNDARYSIZE], char Color[], int time_limit, int game_length, int GameRecord[MAXGAMELENGTH][BOUNDARYSIZE][BOUNDARYSIZE]){
    int turn = (Color[0]=='b'||Color[0]=='B')?BLACK:WHITE;
    int move = genmove(Board, turn, time_limit, game_length, GameRecord);
    int move_i, move_j;
    record(Board, GameRecord, game_length+1);
    if (move==0)
        cout << "= PASS" << endl<< endl<< endl;
    else {
        move_i = (move%100)/10;
        move_j = (move%10);
        cout << "= " << LabelX[move_j]<<10-move_i<<endl<< endl;
    }
}
/*
 * This main function is used of the gtp protocol
 * */
void gtp_main(int display) {
    char Input[COMMANDLENGTH]="";
    char Command[COMMANDLENGTH]="";
    char Parameter[COMMANDLENGTH]="";
    char Move[4]="";
    char Color[6]="";
    int ivalue;
    double dvalue;
    int Board[BOUNDARYSIZE][BOUNDARYSIZE]={{0}};
    int NumCapture[3]={0};// 1:Black, 2: White
    int time_limit = DEFAULTTIME;
    int GameRecord[MAXGAMELENGTH][BOUNDARYSIZE][BOUNDARYSIZE]={{{0}}};
    int game_length = 0;
    if (display==1) {
    gtp_list_commands();
    gtp_showboard(Board);
    }
    while (gets(Input) != 0) {
    sscanf(Input, "%s", Command);
    if (Command[0]== '#')
        continue;

    if (strcmp(Command, "protocol_version")==0) {
        gtp_protocol_version();
    }
    else if (strcmp(Command, "name")==0) {
        gtp_name();
    }
    else if (strcmp(Command, "version")==0) {
        gtp_version();
    }
    else if (strcmp(Command, "list_commands")==0) {
        gtp_list_commands();
    }
    else if (strcmp(Command, "known_command")==0) {
        sscanf(Input, "known_command %s", Parameter);
        gtp_known_command(Parameter);
    }
    else if (strcmp(Command, "boardsize")==0) {
        sscanf(Input, "boardsize %d", &ivalue);
        gtp_boardsize(ivalue);
    }
    else if (strcmp(Command, "clear_board")==0) {
        gtp_clear_board(Board, NumCapture);
        game_length = 0;
    }
    else if (strcmp(Command, "komi")==0) {
        sscanf(Input, "komi %lf", &dvalue);
        gtp_komi(dvalue);
    }
    else if (strcmp(Command, "play")==0) {
        sscanf(Input, "play %s %s", Color, Move);
        gtp_play(Color, Move, Board, game_length, GameRecord);
        game_length++;
        if (display==1) {
        gtp_showboard(Board);
        }
    }
    else if (strcmp(Command, "genmove")==0) {
        sscanf(Input, "genmove %s", Color);
        gtp_genmove(Board, Color, time_limit, game_length, GameRecord);
        game_length++;
        if (display==1) {
        gtp_showboard(Board);
        }
    }
    else if (strcmp(Command, "quit")==0) {
        break;
    }
    else if (strcmp(Command, "showboard")==0) {
        gtp_showboard(Board);
    }
    else if (strcmp(Command, "undo")==0) {
        game_length--;
        gtp_undo(Board, game_length, GameRecord);
        if (display==1) {
        gtp_showboard(Board);
        }
    }
    else if (strcmp(Command, "final_score")==0) {
        if (display==1) {
        gtp_showboard(Board);
        }
        gtp_final_score(Board);
    }
    }
}
int main(int argc, char* argv[]) {
//    int type = GTPVERSION;// 1: local version, 2: gtp version
    int type = GTPVERSION;// 1: local version, 2: gtp version
    int display = 0; // 1: display, 2 nodisplay
    srand(time(NULL));
    if (argc > 1) {
    if (strcmp(argv[1], "-display")==0) {
        display = 1;
    }
    if (strcmp(argv[1], "-nodisplay")==0) {
        display = 0;
    }
    }
    gtp_main(display);
    return 0;
}
