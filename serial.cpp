#include <iostream>
#include <vector>
using namespace std;
class Solution {
public:
    void gameOfLife(vector<vector<int>>& board) {
        int n = board.size();
        if(n == 0)
            return;
        int m = board[0].size();
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < m; j++) {
                if(board[i][j] > 0) {
                    for(int r = i-1; r <= i+1; r++) {
                        for(int c = j-1; c <= j+1; c++) {
                            if(r == i && c == j)
                                continue;
                            update(board, r, c);
                        }
                    }
                }
            }
        }
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < m; j++) {
                if(board[i][j] == -3 || board[i][j] == 3 || board[i][j] == 4)
                    board[i][j] = 1;
                else
                    board[i][j] = 0;
            }
        }
    }
private:
    void update(vector<vector<int>>& board, int i, int j) {
        int n = board.size();
        int m = board[0].size();
        if(i >= 0 && i < n  && j >= 0 && j < m) {
            if(board[i][j] > 0)
                board[i][j]++;
            else
                board[i][j]--;
        }
    }
};

int main() {
    int size;
    int gens;
    int num_ghost;
    cin >> size >> gens >> num_ghost;
    vector<vector<int>> board(size, vector<int>(size));
    char cell;
    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            cin >> cell;
            board[i][j] = cell == '#' ? 1 : 0;
        }
    }
    Solution sol;
    for(int i = 0; i < gens; i++)
        sol.gameOfLife(board);
    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            if(board[i][j] == 1)
                cout << "#";
            else
                cout << '.';
        }
        cout << endl;
    }
}